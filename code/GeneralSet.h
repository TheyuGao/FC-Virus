#ifndef GENERALSET_H
#define GENERALSET_H



#include <vector>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <map>

using namespace std;

//k-mer
extern int first_p;
extern int second_p;
extern int dep_search;
extern int g_kmer_length;
extern float g_ratio;
extern int g_re_gene;
extern float g_com_kmer_ratio;
extern int g_window_length;
extern int avr_kmer_dep;



extern std::map<std::string,int>kmer_map;
extern std::map<string, vector<pair<int, int> > >k_r_p;

//reads

extern vector<string> data_str;
extern vector<int >data_dep;
extern bool g_is_paired_end; //false;
extern int g_fr_strand;
extern bool g_double_stranded_mode; //false;




//contig

extern std::vector<bool> had_used_read;
extern int break_const;
extern float break_ratio;
//file
extern string g_reads_file;
extern string g_left_file;
extern string g_right_file;
extern string out_dir; 
extern string g_file_type;

//others
extern bool g_help;


#define OPT_DOUBLE_STRANDED_MODE		303
#define OPT_FR_STRAND		304
#define OPT_LEFT		308
#define OPT_RIGHT	309
#define OPT_SINGLEFILE			311



#endif

