#include "loadreads.h"
#include "kmer.h"
#include "GeneralSet.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include<algorithm>
using namespace std;

const int MAX_STR = 1024;

bool check_read(string read)
{

    for (unsigned int i = 0; i < read.size(); ++i)
    {

        unsigned char c = read[i];

        if (c != 'A' && c != 'T' && c != 'G' && c != 'C' && c != 'a' && c != 't' && c != 'g' && c != 'c'&& c != 'N'&& c != 'n')
        {
            return (false);
        }
    }

    return (true);
}



string revcomp (const string kmer) {

    string revstring;

    for(int i = kmer.size() -1; i >= 0;i--) {

        char c = kmer[i];
        char revchar;

        switch (c) {

            case 'g':
                revchar = 'c';
                break;

            case 'G':
                revchar = 'C';
                break;

            case 'a':
                revchar = 't';
                break;

            case 'A':
                revchar = 'T';
                break;

            case 't':
                revchar = 'a';
                break;

            case 'T':
                revchar = 'A';
                break;

            case 'c':
                revchar = 'g';
                break;

            case 'C':
                revchar = 'G';
                break;

            default:
                revchar = 'N';
        }

        revstring += revchar;

    }

    return (revstring);

}

void load_reads() {
    map<string,int >input_data;
 if (!g_is_paired_end) {
    if(g_file_type=="fa"){
        load_reads_fa( g_reads_file,input_data,false);
    }
    else {
        load_reads_fq(g_reads_file,input_data,false);
    }
    } else {
        if (g_double_stranded_mode) {
            if(g_file_type=="fa"){
                load_reads_fa( g_left_file,input_data,false);
                load_reads_fa( g_right_file,input_data,false);
            }
            else {
                load_reads_fq( g_left_file,input_data,false);
                load_reads_fq( g_right_file,input_data,false);
            }

        } else {
            if (g_fr_strand == 2) {//--1-->  <--2--
            if(g_file_type=="fa"){
                load_reads_fa( g_left_file,input_data,false);
                load_reads_fa( g_right_file,input_data,true);
            }
            else {
                load_reads_fq( g_left_file,input_data,false);
                load_reads_fq( g_right_file,input_data,true);
            }
            }
            if (g_fr_strand == 1) {//<--1-- --2-->
            if(g_file_type=="fa"){
                load_reads_fa( g_left_file,input_data,true);
                load_reads_fa( g_right_file,input_data,false);
            }
            else {
                load_reads_fq( g_left_file,input_data,true);
                load_reads_fq( g_right_file,input_data,false);
            }
            }
            if (g_fr_strand == 3) {//--1--> --2-->  or <--1-- <--2--
            if(g_file_type=="fa"){
                load_reads_fa( g_left_file,input_data,false);
                load_reads_fa( g_right_file,input_data,false);
            }
            else {
                load_reads_fq( g_left_file,input_data,false);
                load_reads_fq( g_right_file,input_data,false);
            }
            }
        }
    }
for(map<string,int >::iterator it =input_data.begin();it!=input_data.end();it++){
    data_str.push_back(it->first);
    data_dep.push_back(it->second);
}
    cout<<input_data.size()<<" reads have found!"<<endl;

}



void load_reads_fa(string file, map<string,int >& input_data, bool rev) {

    time_t s_time = time(NULL);

    fstream in;
    in.open(file.c_str(), fstream::in);

    if (!in.is_open()) {

        cout << "Error! Can't open file " << file << endl;
        exit(1);

    }

    cout << "Loading reads from file " << file << " ..." << endl;

    char temp[MAX_STR];
    string read;
    in.getline(temp, MAX_STR);
    int read_num=0;

    while(getline(in,read))
    {
        if(read[0]=='>')
        
        {
            getline(in,read);
            if(read[0]=='>') { continue; }
            if (rev)
                read = revcomp(read);
            if(input_data.find(read)==input_data.end()) {
                input_data.insert(make_pair(read, 1));
            }
            else
                input_data[read]=input_data[read]+1;
            read_num++;
        }
        

    }



    cout <<"Total load " << read_num << " reads!" << endl;
    in.close();

    time_t e_time = time(NULL);
    cout << "Success! (total cost time: " << (e_time - s_time) << " s)" << endl;

}


void load_reads_fq(string file, map<string,int >& input_data, bool rev) {

	time_t s_time = time(NULL);
	
	fstream in;
	in.open(file.c_str(), fstream::in);

	if (!in.is_open()) {

		cout << "Error! Can't open file " << file << endl;
		exit(1);

	}

	cout << "Loading reads from file " << file << " ..." << endl;

	char temp[MAX_STR];
	string read;
	in.getline(temp, MAX_STR);

    while(getline(in,read))
{
    if(read[0]=='@')
    {
        getline(in,read);
        if(read[0]=='@') { continue; }
        if (rev)
            read = revcomp(read);
        if(input_data.find(read)==input_data.end()) {
            input_data.insert(make_pair(read, 1));
        }
        else
            input_data[read]=input_data[read]+1;
    }
}

	in.close();

	time_t e_time = time(NULL);
	cout << "Success! (total cost time: " << (e_time - s_time) << " s)" << endl;

}


