#ifndef LI_H
#define LI_H

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
string construct_contig_reverse(int kmer_length, vector<int> kmer_position,int max_read, string &first_kmer2, string &second_kmer2,vector<bool> &had_used_read);
string Extend_L (int kmer_length,int seed_read,vector<bool> &had_used_read);
string final_contig_Left(int kmer_length, string contig);
#endif