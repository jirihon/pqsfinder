/**
 * Implementation of PQS search algorithm.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/01/17
 * Package: pqsfinder
 */

#include <Rcpp.h>
#include <string>
#include <algorithm>
#include <boost/regex.hpp>
#ifdef _GLIBCXX_DEBUG
#include <google/profiler.h>
#endif

using namespace Rcpp;
using namespace std;


/*
 * Important notes to C++ iterator semantics that are extensively used in this code:
 *
 * string::end method returns an iterator pointing to the __past-the-end__ character
 * of the string! The past-the-end character is a theoretical character that would
 * follow the last character in the string. It __shall not be dereferenced.__
 *
 * The reason for that is because the ranges used by functions of the standard library
 * __do not include__ the element pointed by their closing iterator.
 */


// Implementation constants
const int RUN_CNT = 4;
const int CHECK_INT_PERIOD = 2e7;


/**
 * Cache for low complexity regions. It is usefull just for dealing with almost
 * G-complete sequence regions, which does not usually occur in human genome.
 */
class cache {
public:
  class entry {
  public:
    vector<int> density;
    int score;
    int len;
    entry(int max_len) : density(max_len), score(0), len(0) {}
  };
  typedef std::map<string, cache::entry>::iterator iterator;
  typedef std::map<string, cache::entry>::value_type value_type;
  static const int use_treshold = 1000;

  std::map<string, cache::entry> table;
  int max_len;

  cache(const int max_len) : max_len(max_len) {}

  inline entry *get(const std::string::const_iterator &s, const std::string::const_iterator &e) {
    cache::iterator it = this->table.find(std::string(s, e));
    if (it == this->table.end())
      return NULL;
    else
      return &((*it).second);
  }
  inline void put(const std::string::const_iterator &s, const std::string::const_iterator &e,
                  const cache::entry &entry) {
    std::string seq = std::string(s, e);
    cache::iterator it = this->table.find(seq);
    if (it == this->table.end())
      this->table.insert(cache::value_type(seq, entry));
    else if (entry.density[0] > (*it).second.density[0])
      (*it).second = entry;
  }
};


class results {
public:
  vector<int> start;
  vector<int> len;
  vector<int> score;
  int *density;

  int min_score;

  results(const int seq_len, const int min_score) : min_score(min_score) {
    this->density = (int *)calloc(seq_len, sizeof(int));
    if (this->density == NULL)
      stop("Unable to allocate memory for results density vector.");
  }
  ~results() {
    if (this->density != NULL)
      free(this->density);
  }
  inline void save(const int start, const int len, const int score) {
    if (score >= this->min_score) {
      this->start.push_back(start + 1);
      this->len.push_back(len);
      this->score.push_back(score);
    }
  }
  inline void print(const string::const_iterator &ref) const {
    Rcout << "Results" << endl;
    for (unsigned i = 0; i < this->start.size(); i++) {
      Rcout << "PQS[" << i << "]: " << this->start[i] << " "
            << string(ref + this->start[i], ref + this->start[i] + this->len[i])
            << " " << this->score[i] << endl;
    }
  }
};


class scoring {
public:
  int g_bonus;
  int bulge_penalty;
  Function *user_fn;

  scoring() {
    user_fn = NULL;
  }
  ~scoring() {
    if (user_fn != NULL)
      delete user_fn;
  }
  void set_user_fn(SEXP user_fn) {
    this->user_fn = new Function(user_fn);
  }
};

typedef struct flags {
  bool use_cache;
  bool use_re;
  bool debug;
} flags_t;

typedef struct opts {
  int max_len;
  int min_score;
  int run_min_len;
  int run_max_len;
  int loop_min_len;
  int loop_max_len;
} opts_t;

// Structure to store the currently best pqs during search
typedef struct pqs {
  string::const_iterator s;
  string::const_iterator e;
  int score;
} pqs_t;

// Class representing one run from quadruplex
class run_match {
public:
  string::const_iterator first;
  string::const_iterator second;
  int length() const {
    return second - first;
  };
};


/**
 * Print quadruplex summary
 *
 * @param m Quadruplex runs
 * @param score Score
 * @param ref Reference point, typically start of sequence
 * @param cnt Counter
 */
inline void print_pqs(const run_match m[], int score, const string::const_iterator ref, const int cnt)
{
  Rcout << cnt << ": " << m[0].first - ref << "[" << string(m[0].first, m[0].second) << "]";
  for (int i = 1; i < RUN_CNT; i++)
    Rcout << string(m[i-1].second, m[i].first) << "[" << string(m[i].first, m[i].second) << "]";
  Rcout << " " << score << endl;
}


/**
 * Check lenghts of runs in quadruplex
 *
 * @param score Quadruplex score
 * @param m Quadruplex runs
 * @param sc Scoring table
 */
inline void check_run_lengths(int &score, const run_match m[])
{
  int w1,w2,w3,w4;
  w1 = m[0].length();
  w2 = m[1].length();
  w3 = m[2].length();
  w4 = m[3].length();
  /*
   * Allowed length combinations:
   * r r r r
   * R r r r
   * r R r r
   * r r R r
   * r r r R
   */
  if ((w1 == w2 && w1 == w3 && w1 == w4) ||
      (w1 >  w2 && w2 == w3 && w2 == w4) ||
      (w1 <  w2 && w1 == w3 && w1 == w4) ||
      (w1 == w2 && w2 <  w3 && w1 == w4) ||
      (w1 == w2 && w1 == w3 && w1 <  w4))
    score++;
}

/**
 * Count number of G's in G-run
 *
 * NOTE: This simple counting implementation allows multiple bulges
 * in one g-run.
 *
 * @param m G-run match
 * @return Number of guanines in G-run
 */
inline int count_g_num(const run_match &m) {
  int cnt = 0;
  for (string::const_iterator s = m.first; s != m.second; s++) {
    if (*s == 'G') cnt++;
  }
  return cnt;
}

/**
 * Check content of runs
 *
 * @param score Quadruplex score
 * @param m Quadruples runs
 * @param sc Scoring table
 */
inline void check_run_content(int &score, const run_match m[], const scoring &sc)
{
  int g1,g2,g3,g4,w1,w2,w3,w4;
  w1 = m[0].length();
  w2 = m[1].length();
  w3 = m[2].length();
  w4 = m[3].length();
  g1 = count_g_num(m[0]);
  g2 = count_g_num(m[1]);
  g3 = count_g_num(m[2]);
  g4 = count_g_num(m[3]);

  /*
   * Allowed combinations:
   * g g g g
   * R g g g
   * g R g g
   * g g R g
   * g g g R
   *
   * Legend: g means that wi = gi
   *         R repesents one longer run
   */
  if (g1 == w1 && g2 == w2 && w3 == g3 && w4 == g4 && w1 == w2 && w2 == w3 && w3 == w4)
    score += g1 * sc.g_bonus; // canonical g-quadruplex
  else if ((g2 == w2 && g1 >= g2 && g2 == g3 && g2 == g4))
    score += g2 * sc.g_bonus - sc.bulge_penalty;
  else if ((g1 == w1 && g1 <= g2 && g1 == g3 && g1 == g4) ||
           (g1 == w1 && g1 == g2 && g1 <= g3 && g1 == g4) ||
           (g1 == w1 && g1 == g2 && g1 == g3 && g1 <= g4))
    score += g1 * sc.g_bonus - sc.bulge_penalty;
  else
    score = 0;
}


/**
 * Check loop lengths
 *
 * @param score Quadruplex score
 * @param m Quadruples runs
 */
inline void check_loop_lengths(int &score, const run_match m[])
{
  int l1, l2, l3;
  l1 = m[1].first - m[0].second;
  l2 = m[2].first - m[1].second;
  l3 = m[3].first - m[2].second;

  int mean = (l1 + l2 + l3)/3;

  int d1, d2, d3;
  d1 = (l1 - mean)*(l1 - mean);
  d2 = (l2 - mean)*(l2 - mean);
  d3 = (l3 - mean)*(l3 - mean);

  int s = sqrt((d1 + d2 + d3)/3);

  score = max(score - mean - s, 0);
}


/**
 * Check user scoring function
 *
 * @param score Quadruplex score
 * @param m Quadruples runs
 * @param sc Scoring table
 * @example user function in R
   my_fn <- function(subject, score, start, width, loop_1, run_2, loop_2, run_3, loop_3, run_4) {
     print(subject)
     cat(score, start, width, loop_1, run_2, loop_2, run_3, loop_3, run_4, "\n")
     return(score)
   }
 */
inline void check_user_fn(int &score, const run_match m[], const scoring &sc, SEXP subject, const string::const_iterator ref)
{
  int start, width, loop_1, run_2, loop_2, run_3, loop_3, run_4;

  start = m[0].first - ref + 1;
  width = m[3].second - m[0].first;
  loop_1 = m[0].second - ref + 1;
  run_2 = m[1].first - ref + 1;
  loop_2 = m[1].second - ref + 1;
  run_3 = m[2].first - ref + 1;
  loop_3 = m[2].second - ref + 1;
  run_4 = m[3].first - ref + 1;

  score = as<int>((*sc.user_fn)(subject, score, start, width, loop_1, run_2, loop_2, run_3, loop_3, run_4));
}


/**
 * Check GC skewness
 *
 * @param score Quadruplex score
 * @param m Quadruples runs
 * @param sc Scoring table
 */
inline void check_gc_skew(int &score, run_match m[])
{
  int gc_skew = 0;

  // TODO: Implement exactly the same GC skew calculation as published in journal
  for (string::const_iterator it = m[0].first; it < m[3].second; ++it) {
    if (*it == 'G')
      ++gc_skew;
    else if (*it == 'C')
      --gc_skew;
  }
  int run_sum_len = m[0].length() + m[1].length() + m[2].length() + m[3].length();
  score = score - (run_sum_len - max(gc_skew, 0));
}

/**
 * Perform run search on particular sequence region
 *
 * @param s Start of region
 * @param e End of region (iterator pointing past-the-end)
 * @param m Match info structure
 * @param run_re_c Run regular expression
 * @param opts Algorithm options and limits
 * @param flags Algorithm flags
 * @return True on success, false otherwise
 */
inline bool find_run(
    string::const_iterator s,
    string::const_iterator e,
    run_match &m,
    const boost::regex &run_re_c,
    const opts_t &opts,
    const flags_t &flags)
{
  static boost::smatch boost_m;
  bool status = false;

  if (flags.use_re) {
    status = boost::regex_search(s, e, boost_m, run_re_c, boost::match_default);
    if (status) {
      m.first = boost_m[0].first;
      m.second = boost_m[0].second;
    }
  } else {
    while (*s != 'G' && s < e) ++s;
    e = min(s + opts.run_max_len, e);
    --e; // <e> points to past-the-end character and as such should not be dereferenced
    while (*e != 'G' && e > s) --e;
    status = (s < e);
    if (status) {
      m.first = s;
      m.second = e + 1; // correction to point on past-the-end character
    }
  }
  return status;
}

/**
 * Recursively idetify 4 consecutive runs making quadruplex
 *
 * @param i Odinal number of quadruplex run
 * @param start Start position for the current run
 * @param end Limit end position for the current run
 * @param m Array of run matches
 * @param run_re_c Compiled run regular expression
 */
void find_all_runs(
    SEXP subject,
    int i,
    string::const_iterator start,
    string::const_iterator end,
    run_match m[],
    const boost::regex &run_re_c,
    const opts_t &opts,
    const flags_t &flags,
    const scoring &sc,
    const string::const_iterator &ref,
    const size_t len,
    string::const_iterator &pqs_start,
    pqs_t &pqs_best,
    cache &ctable,
    cache::entry &pqs_cache,
    int &pqs_cnt,
    results &res)
{
  string::const_iterator s, e;
  int score;
  cache::entry *cache_hit;

  if (i > 0 && start - m[i-1].second < opts.loop_min_len)
    start = min(m[i-1].second + opts.loop_min_len, end); // skip too short loop

  for (s = start; s < end; s++)
  {
    if (i == 0)
    {// Specific code for the first run matching
      if (flags.use_cache && pqs_cache.density[0] > cache::use_treshold)
      {
        Rcout << "Cache get " << string(s, min(s + opts.max_len, end)) << endl;
        cache_hit = ctable.get(s, min(s + opts.max_len, end));

        if (cache_hit != NULL) {
          if (flags.debug)
            Rcout << "Cache hit: " << s - ref  << " " << string(s, s+cache_hit->len)
                  << " " << cache_hit->score << endl;

          for (int k = 0; k < opts.max_len; ++k) {
            res.density[s - ref + k] += cache_hit->density[k];
          }

          if (pqs_best.score && s >= pqs_best.e)
          {// Export PQS because no further overlapping pqs can be found
            res.save(pqs_best.s - ref, pqs_best.e - pqs_best.s, pqs_best.score);
            pqs_best.score = 0;
          }
          if (cache_hit->score > pqs_best.score) {
            pqs_best.score = cache_hit->score;
            pqs_best.s = s;
            pqs_best.e = s + cache_hit->len;
          }
          continue;
        }
      }
      // Reset score of best PQS starting at current position and density info
      pqs_cache.score = 0;
      for (int k = 0; k < opts.max_len; ++k)
        pqs_cache.density[k] = 0;
    }
    for (e = end; s < e && find_run(s, e, m[i], run_re_c, opts, flags); e--)
    {
      if (m[i].length() < opts.run_min_len)
        break; // skip too short G-run, try next position

      // Update search bounds
      s = m[i].first;
      e = m[i].second;

      if (m[i].length() > opts.run_max_len) {
        e = s + opts.run_max_len; // skip too long G-runs
        continue;
      }
      if (i > 0 && s - m[i-1].second > opts.loop_max_len)
        return; // skip too long loops

      if (i == 0)
        // Enforce G4 total length limit to be relative to the first G-run start
        find_all_runs(subject, i+1, e, min(s + opts.max_len, end), m, run_re_c, opts, flags, sc,
                   ref, len, s, pqs_best, ctable, pqs_cache, pqs_cnt, res);
      else if (i < 3)
        find_all_runs(subject, i+1, e, end, m, run_re_c, opts, flags, sc,
                   ref, len, pqs_start, pqs_best, ctable, pqs_cache, pqs_cnt, res);
      else {
        /* Check user interrupt after reasonable amount of PQS identified to react
         * on important user signals. I.e. he might want to abort the computation. */
        if (++pqs_cnt == CHECK_INT_PERIOD)
        {
          pqs_cnt = 0;
          checkUserInterrupt();
          Rcout << "Search status: " << ceilf((m[0].first - ref)/(float)len*100) << " %\r";
        }

        if (pqs_best.score && pqs_start >= pqs_best.e)
        {// Export PQS because no further overlapping pqs can be found
          res.save(pqs_best.s - ref, pqs_best.e - pqs_best.s, pqs_best.score);
          pqs_best.score = 0;
        }

        score = 0;
        check_run_lengths(score, m);
        if (score)
          check_run_content(score, m, sc);
        if (score)
          check_loop_lengths(score, m);
        if (score && sc.user_fn != NULL)
          check_user_fn(score, m, sc, subject, ref);

        if (score) {
          // Current PQS satisfied all constraints.

          for (int k = 0; k < e - pqs_start; ++k)
            ++pqs_cache.density[k];

          if (score > pqs_cache.score) {
            // Update properties of caching candidate
            pqs_cache.score = score;
            pqs_cache.len = e - pqs_start;
          }
          if (score > pqs_best.score) {
            // Update properties of best PQS to be exported
            pqs_best.score = score;
            pqs_best.s = pqs_start;
            pqs_best.e = e;
          }
          if (flags.debug)
            print_pqs(m, score, ref, pqs_cache.density[0]);
        }
      }
    }
    if (i == 0) {
      if (flags.use_cache && pqs_cache.density[0] > cache::use_treshold)
        ctable.put(s, min(s + opts.max_len, end), pqs_cache);

      // Add locally accumulated density to global density array
      for (int k = 0; k < opts.max_len; ++k)
        res.density[s - ref + k] += pqs_cache.density[k];
    }
  }
}


/**
 * Perform quadruplex search on given DNA sequence.
 *
 * @param seq DNA sequence
 * @param run_re Run regular expression
 * @param sc Scoring options
 * @param opt Algorihtm options
 * @param flags Algorithm flags
 * @param res Results
 */
void pqs_search(
    SEXP subject,
    const string &seq,
    const string &run_re,
    const scoring &sc,
    const opts_t &opts,
    const flags_t &flags,
    results &res)
{
  boost::regex run_re_c(run_re);
  run_match m[RUN_CNT];

  cache ctable(opts.max_len);
  cache::entry pqs_cache(opts.max_len);

  pqs_t pqs_best;
  pqs_best.score = 0;

  string::const_iterator pqs_start;
  int pqs_cnt = 0;

  // Global sequence length is the only limit for the first G-run
  find_all_runs(subject, 0, seq.begin(), seq.end(), m, run_re_c, opts, flags, sc,
                seq.begin(), seq.length(), pqs_start, pqs_best, ctable, pqs_cache, pqs_cnt, res);

  if (pqs_best.score)
    res.save(pqs_best.s - seq.begin(), pqs_best.e - pqs_best.s, pqs_best.score);
}


//' Identificate potential quadruplex forming sequences.
//'
//' Function for identification of all potential intramolecular quadruplex
//' patterns (PQS) in DNA sequence.
//'
//' @param subject DNAString object.
//' @param max_len Maximal lenth of PQS.
//' @param min_score Minimal PQS score.
//' @param run_min_len Minimal length of quadruplex run.
//' @param run_max_len Maximal length of quadruplex run.
//' @param loop_min_len Minimal length of quadruplex loop.
//' @param loop_max_len Maxmimal length of quadruplex loop.
//' @param g_bonus Score bonus for one complete G tetrade.
//' @param bulge_penalty Penalization for a bulge in quadruplex run.
//' @param run_re Regular expression specifying one run of quadruplex.
//' @param user_fn Custom quadruplex scoring function. It takes the following 10
//' arguments: \code{subject} - Input DNAString object, \code{score} - implicit PQS score,
//' \code{start} - PQS start position, \code{width} - PQS width, \code{loop_1} - start pos. of loop #1,
//' \code{run_2} - start pos. of run #2, \code{loop_2} - start pos. of loop #2,
//' \code{run_3} - start pos. of run #3, \code{loop_3} - start pos. of loop #3,
//' \code{run_4} - start pos. of run #4. Return value of the function should be new score
//' represented as asingle integer value.
//' @param use_cache Use cache for low complexity regions?
//' @param use_re Use regular expression engine to validate quadruplex run?
//' @param use_prof Enables profiling.
//' @param debug Enables detailed debugging output. Turn it on if you want
//' to see all possible quadruplexes found at each positions and not just
//' the best one.
//' @return \code{\link{PQSViews}} object
//'
//' @examples
//' pv <- pqsfinder(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
//'
// [[Rcpp::export]]
SEXP pqsfinder(
    SEXP subject,
    int max_len = 70,
    int min_score = 0,
    int run_min_len = 3,
    int run_max_len = 11,
    int loop_min_len = 0,
    int loop_max_len = 30,
    int g_bonus = 20,
    int bulge_penalty = 10,
    std::string run_re = "G{1,5}.{0,5}G{1,5}",
    SEXP user_fn = R_NilValue,
    bool use_cache = 1,
    bool use_re = 0,
    bool use_prof = 0,
    bool debug = 0)
{
  Function as_character("as.character");
  Function get_class("class");

  CharacterVector subject_class = as_character(get_class(subject));

  if (subject_class[0] != "DNAString")
    stop("Subject must be DNAString object.");

  string seq = as<string>(as_character(subject));

  if (run_re != "G{1,5}.{0,5}G{1,5}")
    // User specified its own regexp, force to use regexp engine
    use_re = true;

  flags_t flags;
  flags.use_cache = use_cache;
  flags.use_re = use_re;
  flags.debug = debug;

  opts_t opts;
  opts.max_len = max_len;
  opts.min_score = min_score;
  opts.loop_max_len = loop_max_len;
  opts.loop_min_len = loop_min_len;
  opts.run_max_len = run_max_len;
  opts.run_min_len = run_min_len;

  scoring sc;
  sc.g_bonus = g_bonus;
  sc.bulge_penalty = bulge_penalty;

  results res(seq.length(), opts.min_score);

  if (flags.debug) {
    Rcout << "G-run regexp: " << run_re << endl;
    Rcout << "Use cache: " << flags.use_cache << endl;
    Rcout << "Use regexp engine: " << flags.use_re << endl;
    Rcout << "Input sequence length: " << seq.length() << endl;
    Rcout << "Use user fn: " << (user_fn != R_NilValue) << endl;
  }

  if (user_fn != R_NilValue) {
    sc.set_user_fn(user_fn);
    if (flags.debug) {
      Rcout << "User function: " << endl;
      Function show("show");
      show(user_fn);
    }
  }

  #ifdef _GLIBCXX_DEBUG
  if (use_prof)
    ProfilerStart("samples.log");
  #endif

  pqs_search(subject, seq, run_re, sc, opts, flags, res);

  #ifdef _GLIBCXX_DEBUG
  if (use_prof)
    ProfilerStop();
  #endif

  NumericVector res_start(res.start.begin(), res.start.end());
  NumericVector res_width(res.len.begin(), res.len.end());
  NumericVector res_score(res.score.begin(), res.score.end());

  NumericVector res_density(seq.length());
  for (unsigned i = 0; i < seq.length(); ++i)
    res_density[i] = res.density[i];

  Function pqsviews("PQSViews");
  return pqsviews(subject, res_start, res_width, res_score, res_density);
}
