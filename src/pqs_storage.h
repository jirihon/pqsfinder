/**
 * Storage class aggregating best overlapping PQS.
 *
 * Author: Jiri Hon <jiri.hon@gmail.com>
 * Date: 2016/03/30
 * Package: pqsfinder
 */

#ifndef PQS_STORAGE_HEADER
#define PQS_STORAGE_HEADER

#include <Rcpp.h>
#include "results.h"

using namespace Rcpp;
using namespace std;

class pqs_storage {
private:
  class range {
  public:
    string::const_iterator s;
    string::const_iterator e;
    range(string::const_iterator s, string::const_iterator e) :
      s(s), e(e) {};
    range() {};
  };
  typedef map< int, list<range> > storage_t;
  storage_t st;

public:
  typedef struct pqs {
    int score;
    string::const_iterator s;
    string::const_iterator e;
  } pqs_t;
  pqs_t best;

  pqs_storage() {
    this->best.score = 0;
  }
  inline void insert(
      int score, string::const_iterator s,
      string::const_iterator e)
  {
    if (score > this->best.score ||
        (score == this->best.score &&
         this->best.s <= s && e <= this->best.e)) {
      this->best.score = score;
      this->best.s = s;
      this->best.e = e;
    }
    storage_t::iterator it = st.find(score);
    if (it != st.end()) {
      list<range> &list = it->second;
      if (list.empty()) {
        list.push_back(range(s, e));
      }
      else {
        range &last = list.back();
        if (last.s <= s && e <= last.e) {
          // Replace existing by shorter one
          last.s = s;
          last.e = e;
        }
        else if (last.e <= s) {
          // Insert new non-overlapping pqs
          list.push_back(range(s, e));
        }
      }
    }
    else {
      st.insert(storage_t::value_type(score, list<range>(1, range(s, e))));
    }
  }
  inline void export_pqs(results &res, const string::const_iterator &ref,
                   const string &strand) {
    this->best.score = 0; // reset

    range best_pqs;
    storage_t::iterator it;
    storage_t::iterator temp;

    while (!st.empty()) {
      it = --st.end(); // decrement to point on last list
      best_pqs = it->second.back();

      res.save_pqs(it->first, best_pqs.s, best_pqs.e, ref, strand);

      while (true) {
        // remove all overlapping PQS with lower score
        // Rcout << "Removing from list " << it->first << endl;
        list<range> &list = it->second;
        while (!list.empty() &&
               ((list.back().s <= best_pqs.s && best_pqs.s < list.back().e) ||
               (best_pqs.s <= list.back().s && list.back().s < best_pqs.e)))  {
          list.pop_back();
        }
        if (it == st.begin()) {
          // the end of iteration
          if (list.empty())
            st.erase(it);
          break;
        }
        else if (list.empty()) {
          // delete empty list from storage
          temp = it; // erase operation invalidates iterator
          --it;
          st.erase(temp);
        }
        else {
          --it;
        }
      }
    }
  }
  inline void print() {
    for (storage_t::const_iterator it = st.begin(); it != st.end(); ++it) {
      Rcout << it->first << ":";
      for (list<range>::const_iterator lit = it->second.begin(); lit != it->second.end(); ++lit) {
        Rcout << " " << string(lit->s, lit->e);
      }
      Rcout << endl;
    }
  }
};

#endif // PQS_STORAGE_HEADER
