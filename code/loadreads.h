#ifndef LOADREADS_H
#define LOADREADS_H


#include <vector>
#include <map>
#include <string>


using namespace std;


typedef vector<unsigned long long> read_int_type;
bool check_read(string read);
void load_reads();
void load_reads_fa(string file, map<string,int>& input_data, bool rev);
void load_reads_fq(string file, map<string,int>& input_data, bool rev);

#endif