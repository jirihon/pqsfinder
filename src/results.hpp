/**
 * Storage class for results.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/03/30
 * Package: pqsfinder
 */

#ifndef RESULTS_HEADER
#define RESULTS_HEADER

#include <Rcpp.h>
#include <cstdlib>

using namespace Rcpp;
using namespace std;

class results {
public:
  vector<int> start;
  vector<int> len;
  vector<int> score;
  vector<string> strand;
  int *density;

  const int min_score;
  const int seq_len;

  results(const int seq_len, const int min_score) :
    min_score(min_score), seq_len(seq_len)
  {
    this->density = (int *)calloc(seq_len, sizeof(int));
    if (this->density == NULL)
      stop("Unable to allocate memory for results density vector.");
  }
  ~results() {
    if (this->density != NULL)
      free(this->density);
  }
  inline void save_pqs(
      const int score, const string::const_iterator &s,
      const string::const_iterator &e, const string::const_iterator &ref,
      const string &strand)
  {
    if (score >= this->min_score) {
      if (strand == "+")
        this->start.push_back(s - ref + 1); // R indexing starts at 1
      else
        this->start.push_back(this->seq_len - (e - ref) + 1);

      this->len.push_back(e - s);
      this->score.push_back(score);
      this->strand.push_back(strand);
    }
  }
  inline void save_density(
      const string::const_iterator &s, const string::const_iterator &ref,
      const string &strand, const int *density, const int max_len)
  {
    int offset;
    if (strand == "+") {
      offset = s - ref;
      for (int k = 0; k < max_len; ++k)
        this->density[offset + k] += density[k];
    }
    else {
      offset = (this->seq_len - 1) - (s - ref);
      for (int k = 0; k < max_len; ++k)
        this->density[offset - k] += density[k];
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

#endif // RESULTS_HEADER
