#ifndef RI_H
#define RI_H

#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<sstream>
#include <map>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <string.h>
#include "loadreads.h"
#include "GeneralSet.h"
#include"kmer.h"

using namespace std;

//这里map用以存储kmer
int get_seed(int kmer_length);
string construct_contig(int kmer_length, vector<int> kmer_position,int max_read, string &first_kmer2, string &second_kmer2,vector<bool> &had_used_read);
string Extend_R (int kmer_length,int seed_read,vector<bool> &had_used_read);
string final_contig(int kmer_length, string contig);
#endif
