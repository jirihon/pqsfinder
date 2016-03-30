/**
 * Implementation of PQS search algorithm.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/01/17
 * Package: pqsfinder
 */

#include <Rcpp.h>
#include <string>
#include <cstring>
#include <cstdlib>
#include <algorithm>
#include <boost/regex.hpp>
#ifdef _GLIBCXX_DEBUG
#include <google/profiler.h>
#endif
#include "results.h"
#include "pqs_storage.h"
#include "pqs_cache.h"

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
static const int RUN_CNT = 4;
static const int CHECK_INT_PERIOD = 2e7;


class scoring {
public:
  int tetrad_bonus;
  int bulge_penalty;
  int mismatch_penalty;
  int loop_mean_factor;
  int loop_sd_factor;
  int max_bulges;
  int max_mimatches;
  int max_defects;
  Function *custom_scoring_fn;

  scoring() {
    custom_scoring_fn = NULL;
  }
  ~scoring() {
    if (custom_scoring_fn != NULL)
      delete custom_scoring_fn;
  }
  void set_custom_scoring_fn(SEXP custom_scoring_fn) {
    this->custom_scoring_fn = new Function(custom_scoring_fn);
  }
};

typedef struct flags {
  bool use_cache;
  bool use_re;
  bool use_prof;
  bool verbose;
  bool debug;
  bool use_default_scoring;
} flags_t;

typedef struct opts {
  int max_len;
  int min_score;
  int run_min_len;
  int run_max_len;
  int loop_min_len;
  int loop_max_len;
} opts_t;

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
  Rcout << m[0].first - ref + 1 << " " << cnt << " "  << "[" << string(m[0].first, m[0].second) << "]";
  for (int i = 1; i < RUN_CNT; i++)
    Rcout << string(m[i-1].second, m[i].first) << "[" << string(m[i].first, m[i].second) << "]";
  Rcout << " " << score << endl;
}


/**
 * Count number of G's in G-run.
 * This way of counting leads to only one bulge allowed in G-run.
 *
 * @param m G-run match
 * @return Number of guanines in G-run
 */
inline int count_g_num(const run_match &m) {
  string::const_iterator s = m.first, e = m.second;
  int cnt = 0;
  while (*s == 'G' && s < e) { ++s; ++cnt; }
  --e; // <e> points to past-the-end character and as such should not be dereferenced
  while (*e == 'G' && e > s) { --e; ++cnt; }
  return cnt;
}


/**
 * R interface to test count_g_num function
 *
 * @param seq G-run sequence
 * @return Number of canonical Gs
 */
// ![[Rcpp::export]]
void count_g(std::string seq) {
  run_match m;
  m.first = seq.begin();
  m.second = seq.end();
  int cnt = count_g_num(m);
  Rcout << cnt << endl;
}


/**
 * Check content of runs
 *
 * @param score Quadruplex score
 * @param m Quadruples runs
 * @param sc Scoring table
 */
inline void score_run_content(int &score, const run_match m[], const scoring &sc)
{
  int w[RUN_CNT], g[RUN_CNT];
  int pi = -1, mismatches = 0, bulges = 0, perfects = 0;

  w[0] = m[0].length();
  w[1] = m[1].length();
  w[2] = m[2].length();
  w[3] = m[3].length();

  /*
   * Allowed length combinations:
   * r r r r
   * R r r r
   * r R r r
   * r r R r
   * r r r R
   */
  if (sc.max_bulges == 1 && !(
       (w[0] == w[1] && w[0] == w[2] && w[0] == w[3]) ||
       (w[0] >  w[1] && w[1] == w[2] && w[1] == w[3]) ||
       (w[0] <  w[1] && w[0] == w[2] && w[0] == w[3]) ||
       (w[0] == w[1] && w[1] <  w[2] && w[0] == w[3]) ||
       (w[0] == w[1] && w[0] == w[2] && w[0] <  w[3])
      )) {
    score = 0;
    return;
  }

  for (int i = 0; i < RUN_CNT; ++i) {
    g[i] = count_g_num(m[i]);
    if (g[i] == w[i] && (pi == -1 || w[i] < w[pi]))
      pi = i;
  }
  if (pi < 0)
  {// at least one run has to be perfect
    score = 0;
    return;
  }
  for (int i = 0; i < RUN_CNT; ++i) {
    if (w[i] == w[pi] && g[i] == g[pi])
      ++perfects;
    else if (w[i] == w[pi] && g[i] == g[pi] - 1)
      ++mismatches;
    else if (w[i] > w[pi] && g[i] >= g[pi] && *(m[i].first) == 'G' && *(m[i].second - 1) == 'G')
      ++bulges;
    else {
      score = 0;
      return;
    }
  }
  if (mismatches <= sc.max_mimatches && bulges <= sc.max_bulges && mismatches + bulges <= sc.max_defects)
    score = (w[pi] - 1) * sc.tetrad_bonus - mismatches * sc.mismatch_penalty - bulges * sc.bulge_penalty;
  else
    score = 0;
}


/**
 * Check loop lengths
 *
 * @param score Quadruplex score
 * @param m Quadruples runs
 */
inline void score_loop_lengths(int &score, const run_match m[], const scoring &sc)
{
  if (score == 0)
    return;

  int l1, l2, l3;
  l1 = m[1].first - m[0].second;
  l2 = m[2].first - m[1].second;
  l3 = m[3].first - m[2].second;

  int mean = (l1 + l2 + l3)/3;

  int d1, d2, d3;
  d1 = (l1 - mean)*(l1 - mean);
  d2 = (l2 - mean)*(l2 - mean);
  d3 = (l3 - mean)*(l3 - mean);

  int sd = sqrt((d1 + d2 + d3)/2.0);

  score = max(score - sc.loop_mean_factor * mean - sc.loop_sd_factor * sd, 0);
}


/**
 * Check user scoring function
 *
 * @param score Quadruplex score
 * @param m Quadruples runs
 * @param sc Scoring table
 * @example user function in R
   my_fn <- function(
      subject, score, start, width, loop_1, run_2, loop_2,
      run_3, loop_3, run_4)
   {
      len <- loop_1 - start
      if (len == loop_2 - run_2 && len == loop_3 - run_3 &&
          len == start + width - run_4)
      return(200)
   }
 */
inline void check_custom_scoring_fn(
    int &score, const run_match m[], const scoring &sc, SEXP subject,
    const string::const_iterator ref)
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

  score = as<int>((*sc.custom_scoring_fn)(
    subject, score, start, width,
    loop_1, run_2, loop_2, run_3, loop_3, run_4));
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
    const string::const_iterator start,
    const string::const_iterator end,
    run_match &m,
    const boost::regex &run_re_c,
    const opts_t &opts,
    const flags_t &flags)
{
  static boost::smatch boost_m;
  bool status = false;

  string::const_iterator s = start, e;

  if (flags.use_re) {
    try {
      status = boost::regex_search(start, end, boost_m, run_re_c, boost::match_default);
    } catch (bad_alloc &ba) {
      stop(string("Regexp engine failed with exception: ") + ba.what());
    }
    if (status) {
      m.first = boost_m[0].first;
      m.second = boost_m[0].second;
    }
  } else {
    while (*s != 'G' && s < end) ++s;
    e = min(s + opts.run_max_len, end);
    --e; // <e> points to past-the-end character and as such should not be dereferenced
    while (*e != 'G' && e > s) --e;

    status = (s < e); // this means that the run contains at least two guanines.
    ++e; // correction to point on past-the-end character
    if ((e - s) + 1 < opts.run_min_len)
    // definitely too short to be run
      status = false;

    if (status) {
      m.first = s;//max(s - 1, start); // if it is possible to extend one mismatch left, do it.
      m.second = e;//min(min(e + 1, s + opts.run_max_len), end); // if it is possible to extend one mismatch right, do it.
    }
  }
  return status;
}


/**
 * Recursively idetify 4 consecutive runs making quadruplex
 *
 * @param subject DNAString object
 * @param strand Strand specification
 * @param i Odinal number of quadruplex run
 * @param start Start position for the current run
 * @param end Limit end position for the current run
 * @param m Array of run matches
 * @param run_re_c Compiled run regular expression
 */
void find_all_runs(
    SEXP subject,
    const string &strand,
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
    pqs_storage &pqs_storage,
    pqs_cache &ctable,
    pqs_cache::entry &pqs_cache,
    int &pqs_cnt,
    results &res)
{
  string::const_iterator s, e, min_e;
  int score;
  pqs_cache::entry *cache_hit;

  if (i > 0 && start - m[i-1].second < opts.loop_min_len)
    start = min(m[i-1].second + opts.loop_min_len, end); // skip too short loop

  for (s = start; s < end; s++)
  {
    if (i == 0)
    {// Specific code for the first run matching
      if (flags.use_cache && pqs_cache.density[0] > pqs_cache::use_treshold)
      {
        cache_hit = ctable.get(s, min(s + opts.max_len, end));

        if (cache_hit != NULL) {
          if (flags.debug)
            Rcout << "Cache hit: " << s - ref  << " " << string(s, s+cache_hit->len)
                  << " " << cache_hit->score << endl;

          res.save_density(s, ref, strand, cache_hit->density, opts.max_len);

          if (pqs_storage.best.score && s >= pqs_storage.best.e)
          {// Export PQS because no further overlapping pqs can be found
            pqs_storage.export_pqs(res, ref, strand);
          }
          pqs_storage.insert(cache_hit->score, s, s + cache_hit->len);
          continue;
        }
      }
      // Reset score of best PQS starting at current position and density info
      pqs_cache.score = 0;
      for (int k = 0; k < opts.max_len; ++k)
        pqs_cache.density[k] = 0;
    }
    min_e = s + opts.run_min_len;

    for (e = end; e >= min_e && find_run(s, e, m[i], run_re_c, opts, flags); e--)
    {
      if (m[i].length() < opts.run_min_len)
        break; // skip too short G-run, try next position

      // Update search bounds
      s = m[i].first;
      e = m[i].second;

      if (m[i].length() > opts.run_max_len) {
        e = s + opts.run_max_len; // skip too long G-runs, can be improved in <find_run>
        continue;
      }
      if (i > 0 && s - m[i-1].second > opts.loop_max_len)
        return; // skip too long loops

      if (i == 0)
        // Enforce G4 total length limit to be relative to the first G-run start
        find_all_runs(
          subject, strand, i+1, e, min(s + opts.max_len, end), m, run_re_c,
          opts, flags, sc, ref, len, s, pqs_storage, ctable, pqs_cache,
          pqs_cnt, res
        );
      else if (i < 3)
        find_all_runs(
          subject, strand, i+1, e, end, m, run_re_c, opts, flags, sc, ref, len,
          pqs_start, pqs_storage, ctable, pqs_cache, pqs_cnt, res
        );
      else {
        /* Check user interrupt after reasonable amount of PQS identified to react
         * on important user signals. I.e. he might want to abort the computation. */
        if (++pqs_cnt == CHECK_INT_PERIOD)
        {
          pqs_cnt = 0;
          checkUserInterrupt();
          if (!flags.verbose)
            Rcout << "Search status: " << ceilf((m[0].first - ref)/(float)len*100) << " %\r" << flush;
        }

        if (pqs_storage.best.score && pqs_start >= pqs_storage.best.e)
        {// Export PQS because no further overlapping pqs can be found
          pqs_storage.export_pqs(res, ref, strand);
        }

        score = 0;
        if (flags.use_default_scoring) {
          score_run_content(score, m, sc);
          score_loop_lengths(score, m, sc);
        }
        if ((score || !flags.use_default_scoring) && sc.custom_scoring_fn != NULL) {
          check_custom_scoring_fn(score, m, sc, subject, ref);
        }

        if (score && score >= opts.min_score) {
          // Current PQS satisfied all constraints.
          pqs_storage.insert(score, pqs_start, e);

          for (int k = 0; k < e - pqs_start; ++k)
            ++pqs_cache.density[k];

          if (score > pqs_cache.score ||
              (score == pqs_cache.score && pqs_cache.len < e - pqs_start)) {
            // Update properties of caching candidate
            pqs_cache.score = score;
            pqs_cache.len = e - pqs_start;
          }
          if (flags.verbose)
            print_pqs(m, score, ref, pqs_cache.density[0]);
        }
      }
    }
    if (i == 0) {
      if (flags.use_cache && pqs_cache.density[0] > pqs_cache::use_treshold)
        ctable.put(s, min(s + opts.max_len, end), pqs_cache);

      // Add locally accumulated density to global density array
      res.save_density(s, ref, strand, pqs_cache.density, opts.max_len);
    }
  }
}


/**
 * Perform quadruplex search on given DNA sequence.
 *
 * @param subject DNAString object
 * @param seq DNA sequence
 * @param strand Strand specification
 * @param run_re_c Run regular expression
 * @param ctable PQS cache table
 * @param sc Scoring options
 * @param opts Algorihtm options
 * @param flags Algorithm flags
 * @param res Results object
 */
void pqs_search(
    SEXP subject,
    const string &seq,
    const string strand,
    const boost::regex &run_re_c,
    pqs_cache &ctable,
    const scoring &sc,
    const opts_t &opts,
    const flags_t &flags,
    results &res)
{
  run_match m[RUN_CNT];
  pqs_cache::entry pqs_cache(opts.max_len);
  string::const_iterator pqs_start;
  int pqs_cnt = 0;
  pqs_storage pqs_storage;

  // Global sequence length is the only limit for the first G-run
  find_all_runs(
    subject, strand, 0, seq.begin(), seq.end(), m, run_re_c, opts, flags, sc,
    seq.begin(), seq.length(), pqs_start, pqs_storage, ctable,
    pqs_cache, pqs_cnt, res
  );
  if (pqs_storage.best.score)
    pqs_storage.export_pqs(res, seq.begin(), strand);
}


//' Identificate potential quadruplex forming sequences.
//'
//' Function for identification of all potential intramolecular quadruplex
//' patterns (PQS) in DNA sequence.
//'
//' @param subject DNAString object.
//' @param strand Strand specification. Allowed values are "+", "-" or "*",
//'   where the last one represents both strands. Implicitly, the input
//'   DNAString object is assumed to encode the "+" strand.
//' @param max_len Maximal lenth of PQS.
//' @param min_score Minimal PQS score.
//' @param run_min_len Minimal length of quadruplex run.
//' @param run_max_len Maximal length of quadruplex run.
//' @param loop_min_len Minimal length of quadruplex loop.
//' @param loop_max_len Maxmimal length of quadruplex loop.
//' @param max_bulges Maximal number of runs with bulge.
//' @param max_mismatches Maximal number of runs with mismatch.
//' @param max_defects Maximum number of defects in total (\code{max_bulges +
//'   max_mismatches}).
//' @param tetrad_bonus Score bonus for one complete G tetrade.
//' @param bulge_penalty Penalization for a bulge in quadruplex run.
//' @param mismatch_penalty Penalization for a mismatch in tetrad.
//' @param loop_mean_factor Penalization factor of loop lengths mean.
//' @param loop_sd_factor Penalization factor of loop lengths standard
//'   deviation.
//' @param run_re Regular expression specifying one run of quadruplex.
//' @param custom_scoring_fn Custom quadruplex scoring function. It takes the
//'   following 10 arguments: \code{subject} - Input DNAString object,
//'   \code{score} - implicit PQS score, \code{start} - PQS start position,
//'   \code{width} - PQS width, \code{loop_1} - start pos. of loop #1,
//'   \code{run_2} - start pos. of run #2, \code{loop_2} - start pos. of loop
//'   #2, \code{run_3} - start pos. of run #3, \code{loop_3} - start pos. of
//'   loop #3, \code{run_4} - start pos. of run #4. Return value of the function
//'   has to be new score represented as a single integer value. Please note
//'   that if \code{use_default_scoring} is enabled, the custom scoring function
//'   is evaluated AFTER the default scoring system but ONLY IF the default
//'   scoring system resulted in non-zero score (for performance reasons). On
//'   the other hand, when \code{use_default_scoring} is disabled, custom
//'   scoring function is evaluated on every PQS.
//' @param use_default_scoring Enables default internal scoring system. This
//'   option is particularly useful in case you intend to radically change the
//'   default behavior and specify your own scoring function. By disabling the
//'   default scoring you will get a full control above the underlying detection
//'   algorithm.
//' @param verbose Enables detailed output. Turn it on if you want to see all
//'   possible PQS found at each positions and not just the best one. It is
//'   highly recommended to use this option for debugging custom quadruplex
//'   scoring function. Each PQS is reported on separate row in the following
//'   format: \code{start cnt pqs_sequence score}, where \code{start} is the PQS
//'   starting position, \code{pqs_sequence} shows the PQS sequence structure
//'   with each run surrounded by square brackets and \code{score} is the score
//'   assigned to the particular PQS by all applied scoring functions.
//' @return \code{\link{PQSViews}} object
//'
//' @examples
//' pv <- pqsfinder(DNAString("CCCCCCGGGTGGGTGGGTGGGAAAA"))
//'
// [[Rcpp::export]]
SEXP pqsfinder(
    SEXP subject,
    std::string strand = "*",
    int max_len = 50,
    int min_score = 42,
    int run_min_len = 3,
    int run_max_len = 11,
    int loop_min_len = 0,
    int loop_max_len = 30,
    int max_bulges = 3,
    int max_mismatches = 2,
    int max_defects = 3,
    int tetrad_bonus = 45,
    int bulge_penalty = 20,
    int mismatch_penalty = 31,
    float loop_mean_factor = 1,
    float loop_sd_factor = 1,
    std::string run_re = "G{1,5}.{0,5}G{1,5}",
    SEXP custom_scoring_fn = R_NilValue,
    bool use_default_scoring = true,
    bool verbose = false)
{
  if (max_len < 1)
    stop("Maximal length of PQS has to be positive value.");
  if (min_score < 1)
    stop("Minimal PQS score has to be positive value.");

  if (run_min_len < 0)
    stop("Minimal PQS run length has to be non-negative value.");
  if (run_max_len < 0)
    stop("Maximal PQS run length has to be non-negative value.");
  if (run_min_len > run_max_len)
    stop("Minimal PQS run length can't be greater than maximal PQS run length.");

  if (loop_min_len < 0)
    stop("Minimal PQS loop length has to be non-negative value.");
  if (loop_max_len < 0)
    stop("Maximal PQS loop length has to be non-negative value.");
  if (loop_min_len > loop_max_len)
    stop("Minimal PQS loop length can't be greater than maximal PQS loop length.");

  if (max_bulges < 0 || max_bulges > 3)
    stop("Maximum number of runs with bulges has to be from the range 0-3.");
  if (max_mismatches < 0 || max_mismatches > 3)
    stop("Maximum number of runs with mismatches has to be from the range 0-3.");
  if (max_defects < 0 || max_defects > 3)
    stop("Maximum number of runs with defects (bulge or mismatch) has to be from the range 0-3.");
  if (strand != "+" && strand != "-" && strand != "*")
    stop("Strand specification must be +, - or *.");

  Function as_character("as.character");
  Function get_class("class");

  CharacterVector subject_class = as_character(get_class(subject));

  if (subject_class[0] != "DNAString")
    stop("Subject must be DNAString object.");

  flags_t flags;
  flags.use_cache = true;
  flags.use_re = false;
  flags.use_prof = false;
  flags.debug = false;
  flags.verbose = verbose;
  flags.use_default_scoring = use_default_scoring;

  if (run_re != "G{1,5}.{0,5}G{1,5}")
    // User specified its own regexp, force to use regexp engine
    flags.use_re = true;

  opts_t opts;
  opts.max_len = max_len;
  opts.min_score = min_score;
  opts.loop_max_len = loop_max_len;
  opts.loop_min_len = loop_min_len;
  opts.run_max_len = run_max_len;
  opts.run_min_len = run_min_len;

  scoring sc;
  sc.tetrad_bonus = tetrad_bonus;
  sc.bulge_penalty = bulge_penalty;
  sc.mismatch_penalty = mismatch_penalty;
  sc.loop_mean_factor = loop_mean_factor;
  sc.loop_sd_factor = loop_sd_factor;
  sc.max_bulges = max_bulges;
  sc.max_mimatches = max_mismatches;
  sc.max_defects = max_defects;

  string seq = as<string>(as_character(subject));
  Function reverseComplement("reverseComplement");
  SEXP subject_rc = reverseComplement(subject);
  string seq_rc = as<string>(as_character(subject_rc));

  results res(seq.length(), opts.min_score);
  pqs_cache ctable(opts.max_len);
  boost::regex run_re_c(run_re);

  if (flags.debug) {
    Rcout << "G-run regexp: " << run_re << endl;
    Rcout << "Use cache: " << flags.use_cache << endl;
    Rcout << "Use regexp engine: " << flags.use_re << endl;
    Rcout << "Input sequence length: " << seq.length() << endl;
    Rcout << "Use user fn: " << (custom_scoring_fn != R_NilValue) << endl;
  }

  if (custom_scoring_fn != R_NilValue) {
    sc.set_custom_scoring_fn(custom_scoring_fn);
    if (flags.debug) {
      Rcout << "User function: " << endl;
      Function show("show");
      show(custom_scoring_fn);
    }
  }

  #ifdef _GLIBCXX_DEBUG
  if (flags.use_prof)
    ProfilerStart("profiling.log");
  #endif

  if (strand == "+" || strand == "*") {
    Rcout << "Searching on sense strand..." << endl;
    pqs_search(subject, seq, "+", run_re_c, ctable, sc, opts, flags, res);
  }
  if (strand == "-" || strand == "*") {
    Rcout << "Searching on antisense strand..." << endl;
    pqs_search(subject_rc, seq_rc, "-", run_re_c, ctable, sc, opts, flags, res);
  }

  #ifdef _GLIBCXX_DEBUG
  if (flags.use_prof)
    ProfilerStop();
  #endif

  NumericVector res_start(res.start.begin(), res.start.end());
  NumericVector res_width(res.len.begin(), res.len.end());
  NumericVector res_score(res.score.begin(), res.score.end());
  CharacterVector res_strand(res.strand.begin(), res.strand.end());

  NumericVector res_density(seq.length());
  for (unsigned i = 0; i < seq.length(); ++i)
    res_density[i] = res.density[i];

  Function pqsviews("PQSViews");
  return pqsviews(subject, res_start, res_width, res_strand, res_score, res_density);
}
