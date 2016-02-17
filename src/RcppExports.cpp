// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// pqsfinder
SEXP pqsfinder(SEXP subject, std::string run_re, int max_len, int run_min_len, int run_max_len, int loop_min_len, int loop_max_len, int g_bonus, int bulge_penalty, bool use_cache, bool use_re, bool use_prof, bool debug);
RcppExport SEXP pqsfinder_pqsfinder(SEXP subjectSEXP, SEXP run_reSEXP, SEXP max_lenSEXP, SEXP run_min_lenSEXP, SEXP run_max_lenSEXP, SEXP loop_min_lenSEXP, SEXP loop_max_lenSEXP, SEXP g_bonusSEXP, SEXP bulge_penaltySEXP, SEXP use_cacheSEXP, SEXP use_reSEXP, SEXP use_profSEXP, SEXP debugSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< SEXP >::type subject(subjectSEXP);
    Rcpp::traits::input_parameter< std::string >::type run_re(run_reSEXP);
    Rcpp::traits::input_parameter< int >::type max_len(max_lenSEXP);
    Rcpp::traits::input_parameter< int >::type run_min_len(run_min_lenSEXP);
    Rcpp::traits::input_parameter< int >::type run_max_len(run_max_lenSEXP);
    Rcpp::traits::input_parameter< int >::type loop_min_len(loop_min_lenSEXP);
    Rcpp::traits::input_parameter< int >::type loop_max_len(loop_max_lenSEXP);
    Rcpp::traits::input_parameter< int >::type g_bonus(g_bonusSEXP);
    Rcpp::traits::input_parameter< int >::type bulge_penalty(bulge_penaltySEXP);
    Rcpp::traits::input_parameter< bool >::type use_cache(use_cacheSEXP);
    Rcpp::traits::input_parameter< bool >::type use_re(use_reSEXP);
    Rcpp::traits::input_parameter< bool >::type use_prof(use_profSEXP);
    Rcpp::traits::input_parameter< bool >::type debug(debugSEXP);
    __result = Rcpp::wrap(pqsfinder(subject, run_re, max_len, run_min_len, run_max_len, loop_min_len, loop_max_len, g_bonus, bulge_penalty, use_cache, use_re, use_prof, debug));
    return __result;
END_RCPP
}