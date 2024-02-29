// This file was modified from a file named common.h in Binpack.
// The original copyright info is Copyright (c) 2013, The Broad Institute, Inc.
// Distributed under the  Distributed under Binpack Software LICENSE.



#include "GeneralSet.h"
#include <cstring>
#include <cstdlib>
#include <errno.h>
#include <map>

using namespace std;

//k-mer
int g_kmer_length=25;
std::map<std::string,int>kmer_map;
std::map<string, vector<pair<int, int> > >k_r_p;
int first_p=1;
int second_p=10;
int dep_search=0;
float g_ratio=0.01;
int g_window_length=1000;
int g_re_gene=150;
float g_com_kmer_ratio=0.3;
int avr_kmer_dep = 0;
//reads

vector<string> data_str;
vector<int >data_dep;
bool g_is_paired_end = true; 
int g_fr_strand = 2;
bool g_double_stranded_mode = false; 


//contig
extern std::vector<bool> had_used_read;
int break_const=1666;
float break_ratio=0.6;
//file
string g_reads_file = "";
string g_left_file = "";
string g_right_file = "";
string out_dir = ""; 
string g_file_type = "";
//others
bool g_help = false;


