#include "GeneralSet.h"
#include "loadreads.h"
#include "kmer.h"
#include "ex_r.h"
#include "ex_l.h"
#include <iostream>
#include <ctime>
#include <errno.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <time.h>
#include <getopt.h>
#include <map>

using namespace std;


struct option opts[] = {
        {"kmer_length",   required_argument,   0,   'k'},
        {"out_dir",        required_argument,   0,   'o'},
        {"file_type",          required_argument,         0,   't'},
        {"help",          no_argument,         0,   'h'},
        {"double_stranded_mode", no_argument,  0,   OPT_DOUBLE_STRANDED_MODE},
        {"fr",     required_argument,   0,   OPT_FR_STRAND},
        {"left",          required_argument,   0,   OPT_LEFT},
        {"right",         required_argument,   0,   OPT_RIGHT},
        {"singlefile",    required_argument,   0,   OPT_SINGLEFILE},
        {0,0,0,0}

};


string usage() {

    stringstream usage_info;
    usage_info
            << endl
            << "===============================================================================" << endl
            << " IsoTree Usage " << endl
            << "===============================================================================" << endl
            << " ** Options: **" <<endl
            << "  -k <int>: length of kmer, default 25. " << endl
            << "  -o <string>: output directory. " << endl
            << "  -p : paired-end reads. " << endl
            << "  -t <string>: type of file, fa or fq. " << endl
            << "  -h : help information. " << endl
            << " If pair end reads: " << endl
            << "  --left <string>: left reads file name (.fasta). " << endl
            << "  --right <string>: right reads file name (.fasta). " << endl
            << " If single end reads: " << endl
            << "  -singlefile <string>: reads file name (.fasta). " << endl
            << "  --double_stranded_mode: indicate the pair-end read is double stranded mode" << endl
            << "  --fr <int>: only used for pair-end reads. 1: --1--> <--2--  2: <--1-- --2-->  3: --1--> --2--> or <--1-- <--2--, default 1. " << endl
            << "===============================================================================" << endl
            << endl;

    return usage_info.str();

}


int parse_options(int argc, char* argv[]) {

    int option_index = 0;
    int next_option;
    do {
        next_option = getopt_long(argc, argv, "k:o:t:h", opts, &option_index);
        switch (next_option) {
            case -1:
                break;
            case 'k':
                g_kmer_length = atoi(optarg);
                break;
            case 'o':
                out_dir = optarg;
                break;
            case 'h':
                g_help = true;
                break;
            case 't':
                g_file_type = optarg;
                break;
            case OPT_DOUBLE_STRANDED_MODE:
                g_double_stranded_mode = true;
                break;
            case OPT_FR_STRAND:
                g_fr_strand = atoi(optarg);
                break;
            case OPT_LEFT:
                g_left_file = optarg;
                break;
            case OPT_RIGHT:
                g_right_file = optarg;
                break;
            case OPT_SINGLEFILE:
                g_reads_file = optarg;
                break;
            default:
                exit(1);
        }

    } while (next_option != -1);

    if (g_help) {
        cout << usage();
        exit (1);
    }



    if (g_kmer_length > 32) {
        cout << "Error: the kmer length should shorter than 32." << endl;
        exit(1);
    }

    if (g_reads_file.length() > 1)
        g_is_paired_end = false;

    if (g_fr_strand != 1 && g_fr_strand != 2 && g_fr_strand != 3) {
        cout << "Error: --fr can only be 1, 2 or 3" << endl;
        exit(1);
    }

    return 0;

}



int main(int argc, char* argv[]){
    time_t s_time = time(NULL);
    int parse_ret = parse_options(argc,argv);
    if (parse_ret)
        return parse_ret;

    load_reads();

    get_kmer(g_kmer_length);
   
     int seed=get_seed(g_kmer_length);


 string contig="";
 vector<bool> had_used_read(data_str.size(),false);
 string file_name =out_dir+"/FC_strains_contig.fa";
     fstream GAP_file;
     GAP_file.open(file_name.c_str(), fstream::out);
     GAP_file<<"+FC_contig"<<endl;
     string contig_r=Extend_R(g_kmer_length,seed,had_used_read);
    string contig_l=Extend_L(g_kmer_length,seed,had_used_read);
       string contig_exr=final_contig(g_kmer_length,contig_r);
    if(dep_search==2){
        contig=contig_exr;
        GAP_file<<contig<<endl;
        cout<<"contig:"<<contig<<endl;
        time_t a_time = time(NULL);
        cout << "Success! (elapsed time: " << (a_time - s_time) << " s)" << endl;
        return 1;
    }
    string Contig_Ex_L = final_contig_Left(g_kmer_length,contig_l);
    contig=Contig_Ex_L+contig_exr;
    cout<<"contig："<<contig<<endl;
    GAP_file<<">contig"<<endl;
    time_t e_time = time(NULL);
    cout << "Success! (elapsed time: " << (e_time - s_time) << " s)" << endl;
    GAP_file<<contig<<endl;
    return 1;
}