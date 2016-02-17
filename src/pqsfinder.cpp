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
 * NOTE: Use Rcpp::checkUserInterrupt() every second
 * Conclusion: Cache is usefull just for dealing with almost G-complete sequences,
 * which does not usually occur in human genome
 */

/*
 * Important notes to C++ iterators semantic
 *
 * string::end
 * Returns an iterator pointing to the __past-the-end__ character of the string.
 *
 * The past-the-end character is a theoretical character that would follow
 * the last character in the string. It __shall not be dereferenced.__
 *
 * Because the ranges used by functions of the standard library __do not include__
 * the element pointed by their closing iterator, this function is often used
 * in combination with string::begin to specify a range including all
 * the characters in the string.
 */

typedef struct results {
  vector<int> start;
  vector<int> len;
  vector<int> score;
  int *density;
} results_t;

void pqs_export(int start, int len, int score, results_t &res) {
  res.start.push_back(start + 1);
  res.len.push_back(len);
  res.score.push_back(score);
}

typedef struct scoring {
  int g_bonus;
  int bulge_penalty;
} scoring_t;

typedef struct cache_entry {
  string seq;
  int score;
  int len;
  int cnt;
} cache_entry_t;

typedef struct pqs {
  string::const_iterator s;
  string::const_iterator e;
  int score;
} pqs_t;

class run_match {
public:
  string::const_iterator first;
  string::const_iterator second;
  int length() {
    return second - first;
  };
};

typedef struct flags {
  bool use_cache;
  bool use_re;
  bool debug;
} flags_t;

typedef struct opts {
  int max_len;
  int run_min_len;
  int run_max_len;
  int loop_min_len;
  int loop_max_len;
} opts_t;

const int CACHE_SIZE = 1024*1024;
const int RUN_CNT = 4;
const int CHECK_INT_PERIOD = 100000;
const int USE_CACHE_TRESHOLD = 1000;

cache_entry_t cache_table[CACHE_SIZE];

inline unsigned cache_hash(string::const_iterator s, const string::const_iterator &e)
{
  unsigned hash = 0;
  while (s < e) {
    hash = 31 * hash + *s;
    ++s;
  }
  return hash % CACHE_SIZE;
}

inline cache_entry_t *cache_get(const string::const_iterator &s, const string::const_iterator &e)
{
  unsigned hash = cache_hash(s, e);
  if (cache_table[hash].len == 0 || (cache_table[hash].seq != string(s, e)))
    return NULL;
  else
    return &cache_table[hash];
}

inline void cache_put(const string::const_iterator &s, const string::const_iterator &e, const cache_entry_t &entry)
{
  unsigned hash = cache_hash(s, e);
  if (entry.cnt > cache_table[hash].cnt) {
    // replace cached entry
    cache_table[hash] = entry;
    cache_table[hash].seq = string(s, e);
  }
}

/**
 * Print quadruplex summary
 *
 * @param m Quadruplex runs
 * @param score Score
 * @param ref Reference point, typically start of sequence
 * @param cnt Counter
 */
void print_pqs(run_match m[], int score, string::const_iterator ref, int cnt)
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
inline void check_run_lengths(int &score, run_match m[])
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
 * This simple counting allows multiple bulges in one g-run.
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
inline void check_run_content(int &score, run_match m[], const scoring_t &sc)
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
inline bool run_search(
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
void match_runs(
    int i,
    string::const_iterator start,
    string::const_iterator end,
    run_match m[],
    const boost::regex &run_re_c,
    const opts_t &opts,
    const flags_t &flags,
    const scoring_t &sc,
    const string::const_iterator &ref,
    string::const_iterator &pqs_start,
    pqs_t &pqs_best,
    cache_entry_t &pqs_cache,
    int &pqs_cnt,
    results_t &res)
{
  string::const_iterator s, e;
  int score;
  cache_entry_t *cache_hit;

  if (i > 0 && start - m[i-1].second < opts.loop_min_len)
    start = min(m[i-1].second + opts.loop_min_len, end); // skip too short loop

  for (s = start; s < end; s++)
  {
    if (i == 0)
    {// Specific code for the first run matching
      if (flags.use_cache && res.density[max((int)(s - ref - 1), 0)] > USE_CACHE_TRESHOLD)
      {
        cache_hit = cache_get(s, min(s + opts.max_len, end));
        if (cache_hit != NULL) {
          if (flags.debug)
            Rcout << "Cache hit: " << s - ref  << " " << string(s, s+cache_hit->len) << " " << cache_hit->score << endl;

          pqs_export(s - ref, cache_hit->len, cache_hit->score, res);
          res.density[s - ref] = cache_hit->cnt;
          continue;
        }
      }
      pqs_cache.score = 0;
    }
    for (e = end; s < e && run_search(s, e, m[i], run_re_c, opts, flags); e--)
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
        match_runs(i+1, e, min(s + opts.max_len, end), m, run_re_c, opts, flags, sc,
                   ref, s, pqs_best, pqs_cache, pqs_cnt, res);
      else if (i < 3)
        match_runs(i+1, e, end, m, run_re_c, opts, flags, sc,
                   ref, pqs_start, pqs_best, pqs_cache, pqs_cnt, res);
      else {
        /* Check user interrupt after reasonable amount of PQS identified to react
         * on important user signals. I.e. he might want to abort the computation. */
        if (++pqs_cnt % CHECK_INT_PERIOD == 0)
        {
          checkUserInterrupt();
          if (flags.debug)
            Rcout << m[0].first - ref << endl;
        }

        if (pqs_best.score && pqs_start >= pqs_best.e)
        {// Export PQS because no further overlapping pqs can be found
          pqs_export(pqs_best.s - ref, pqs_best.e - pqs_best.s, pqs_best.score, res);
          pqs_best.score = 0;
        }

        score = 0;
        check_run_lengths(score, m);
        if (score)
          check_run_content(score, m, sc);

        if (score) {
          ++res.density[pqs_start - ref];

          if (score > pqs_cache.score) {
            pqs_cache.score = score;
            pqs_cache.cnt = res.density[pqs_start - ref];
            pqs_cache.len = e - pqs_start;
          }

          if (score > pqs_best.score) {
            pqs_best.score = score;
            pqs_best.s = pqs_start;
            pqs_best.e = e;
          }
          if (flags.debug)
            print_pqs(m, score, ref, res.density[pqs_start - ref]);
        }
      }
    }
    if (i == 0 && flags.use_cache && res.density[s - ref] > USE_CACHE_TRESHOLD) {
      cache_put(s, min(s + opts.max_len, end), pqs_cache);
    }
  }
}

void print_res(results_t &res, const string::const_iterator &ref) {
  Rcout << "Results" << endl;
  for (unsigned i = 0; i < res.start.size(); i++) {
    Rcout << "PQS[" << i << "]: " << res.start[i] << " "
          << string(ref + res.start[i], ref + res.start[i] + res.len[i]) << " " << res.score[i] << endl;
  }
}

//'
// [[Rcpp::export]]
SEXP pqsfinder(
    SEXP subject,
    std::string run_re = "G{1,5}.{0,5}G{1,5}",
    int max_len = 70,
    int run_min_len = 3,
    int run_max_len = 11,
    int loop_min_len = 0,
    int loop_max_len = 30,
    int g_bonus = 2,
    int bulge_penalty = 3,
    bool use_cache = 1,
    bool use_re = 0,
    bool use_prof = 0,
    bool debug = 0)
{
#ifdef _GLIBCXX_DEBUG
  if (use_prof)
    ProfilerStart("samples.log");
#endif

  Function as_character("as.character");
  string seq = as<string>(as_character(subject));

  Rcout << "G-run regexp: " << run_re << endl;
  boost::regex run_re_c(run_re);

  if (run_re != "G{1,5}.{0,5}G{1,5}")
  // Override user flag
    use_re = true;

  Rcout << "Use cache: " << use_cache << endl;
  Rcout << "Use regexp: " << use_re << endl;
  Rcout << "Debug: " << debug << endl;
  Rcout << "Input length: " << seq.length() << endl;

  run_match m[RUN_CNT];

  flags_t flags;
  flags.use_cache = use_cache;
  flags.use_re = use_re;
  flags.debug = debug;

  opts_t opts;
  opts.loop_max_len = loop_max_len;
  opts.loop_min_len = loop_min_len;
  opts.max_len = max_len;
  opts.run_max_len = run_max_len;
  opts.run_min_len = run_min_len;

  pqs_t pqs_best;
  pqs_best.score = 0;

  cache_entry_t pqs_cache;
  pqs_cache.len = 0;
  // Clear cache table
  for (int i = 0; i < CACHE_SIZE; ++i)
    cache_table[i] = pqs_cache;

  string::const_iterator pqs_start;

  scoring_t sc;
  sc.g_bonus = g_bonus;
  sc.bulge_penalty = bulge_penalty;

  int pqs_cnt = 0;

  results_t res;
  res.density = (int *)calloc(seq.length(), sizeof(int));
  if (res.density == NULL)
    Rcpp::stop("Unable to allocate enough memory for results.");

  // Global sequence length is the only limit for the first G-run
  match_runs(0, seq.begin(), seq.end(), m, run_re_c, opts, flags, sc, seq.begin(), pqs_start, pqs_best, pqs_cache, pqs_cnt, res);
  if (pqs_best.score)
    pqs_export(pqs_best.s - seq.begin(), pqs_best.e - pqs_best.s, pqs_best.score, res);

  // print_res(res, seq.begin());

  NumericVector res_start(res.start.begin(), res.start.end());
  NumericVector res_width(res.len.begin(), res.len.end());
  NumericVector res_score(res.score.begin(), res.score.end());

  NumericVector res_density(seq.length());
  for (unsigned i = 0; i < seq.length(); ++i)
    res_density[i] = res.density[i];

  free(res.density);

#ifdef _GLIBCXX_DEBUG
  if (use_prof)
    ProfilerStop();
#endif

  Function pqsviews("PQSViews");
  return pqsviews(subject, res_start, res_width, res_score, res_density);
}
