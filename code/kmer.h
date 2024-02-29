#ifndef KMER_H
#define KMER_H

#include "loadreads.h"
#include "GeneralSet.h"
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <list>
#include <boost/unordered_map.hpp>

using namespace std;
extern std::map<std::string,vector<std::string> >K_R_map;
extern std::map<std::string,vector<std::string> >R_K_map;
bool contains_non_gatc (string kmer);
bool Find_peak(int kmer_length);
void get_kmer(int kmer_length);
#endif
