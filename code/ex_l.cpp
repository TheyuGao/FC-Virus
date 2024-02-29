#include "GeneralSet.h"
#include "kmer.h"
#include "loadreads.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>
#include <vector>

using namespace std;


string construct_contig_reverse(int kmer_length, vector<int> kmer_position,
                                int max_read, string &first_kmer2,string &second_kmer2, vector<bool> &had_used_read) {
  int insert_size=0;

  string first_kmer;
  string second_kmer;
  map<string, int> insert_map;
  int break_count;
  string temp_contig = "";
  for (int i = 0; i < kmer_position.size() - 1; i++) {
    first_kmer = data_str[max_read].substr(kmer_position[i], kmer_length);
    second_kmer = data_str[max_read].substr(kmer_position[i + 1], kmer_length);
    vector<pair<int, int>> first_kmer_read = k_r_p[first_kmer];

    break_count = 0;
    for (int k = 0; k < first_kmer_read.size(); k++) {
      int search_read = first_kmer_read[k].first;
      int search_read_p = first_kmer_read[k].second;
      for (int j = search_read_p - 1; j > 0; j--) {
        if (break_count >= 1000) {
          break;
        }
        string kmer = data_str[search_read].substr(j, kmer_length);
        if (kmer == second_kmer) {
          break_count += data_dep[search_read];
          insert_size = search_read_p - j;
          string insert_str = data_str[search_read].substr(j, insert_size);
          if (insert_map.count(insert_str) == 0)
            {insert_map.insert({insert_str, 1});}
          else
            {insert_map[insert_str] += data_dep[search_read];}

          break;
        }
      }
    }

    int max_num = 0;
    string max_str = "";
    for (map<string, int>::iterator kmer_it = insert_map.begin();
         kmer_it != insert_map.end(); kmer_it++) {

      if (kmer_it->second > max_num) {
        max_num = kmer_it->second;
        max_str = kmer_it->first;
      }
    }

    temp_contig = max_str + temp_contig;

    insert_map.clear();
    insert_size = 0;
  }
  first_kmer2=first_kmer;
  second_kmer2=second_kmer;
  return temp_contig;
}

string Extend_L(int kmer_length, int seed_read, vector<bool> &had_used_read) {
  string test_contig = "";
  int max_read = seed_read;
  string first_kmer, second_kmer;
  string last_first_kmer;
  string last_second_kmer;
  bool ex_finish = false;
  vector<int> kmer_position;
  int star_position = 0;
  int last_max_read = 0;
  int cut_num = 0;
  string contig = "";
star_position = data_str[max_read].length() - kmer_length;
  while (max_read != -1) {
    kmer_position.clear();
    for (int j = star_position; j >= 0; j--) {
      const string &kmer = data_str[max_read].substr(j, kmer_length);
      if (kmer_map[kmer] >= first_p && kmer_map[kmer] <= second_p) {
        kmer_position.push_back(j);
      }
    }


    int test_lengtha=test_contig.length();

    string candi_contig =
        construct_contig_reverse(kmer_length, kmer_position, max_read, first_kmer,
                         second_kmer, had_used_read);

    if(test_contig.length()==0){cut_num=candi_contig.length();}
    cout<<"<-----left Constructing......"<<endl;
    test_contig = candi_contig+test_contig;

    int temp_length = test_contig.length() - test_lengtha;

    max_read = -1;
    vector<pair<int, int>> first_kmer_read = k_r_p[first_kmer];
    for (int i = 0; i < k_r_p[first_kmer].size(); i++) {
      if (had_used_read[first_kmer_read[i].first] == true)
        continue;
      int search_read = k_r_p[first_kmer][i].first;
      int search_read_p = k_r_p[first_kmer][i].second;
      int start_p;
      for (int j = search_read_p - 1; j >= 0; j--) {
        string kmer = data_str[search_read].substr(j, kmer_length);
        if (kmer == second_kmer) {
          start_p = j;
          int kmer_num = 0;
          for (int k = j - 1; k >= 0; k--) {
            if (kmer_map[data_str[first_kmer_read[i].first].substr(
                    k, kmer_length)] >= first_p &&
                kmer_map[data_str[first_kmer_read[i].first].substr(
                    k, kmer_length)] <= second_p) {
              star_position = k;
              kmer_num++;
            }
          }
          if (kmer_num > 0) {
            max_read = first_kmer_read[i].first;
            star_position = start_p;
            break;
          }
        }
      }
      if (max_read != -1)
        break;
    }
  }
  cout<<cut_num<<endl;
  test_contig.erase(test_contig.length() - cut_num, cut_num);
  return test_contig;
}

string final_contig_Left(int kmer_length, string contig) {
  int min_fre = 2;
  for (int i = 0; i <= contig.length() - kmer_length; i++) {
    string k_mer = contig.substr(i, kmer_length - 1);
    if (kmer_map.find("A"+k_mer ) != kmer_map.end()) {
      kmer_map["A"+k_mer] =0;

    }
    
  }

  while (true) {

    string temp_contig =
        contig.substr(0, kmer_length - 1);
    string bases = "ATGC";
    char max_base = 'A';
    int max_cov = 0;

    for (char base : bases) {
      string temp_base =base + temp_contig;
      if (kmer_map.find(temp_base) != kmer_map.end()) {
        int cov = kmer_map[temp_base];
        if (cov > max_cov) {
          max_cov = cov;
          max_base = base;
        }
      }
    }
    if (max_cov <= min_fre) {
      break;
    }
    contig = max_base+contig;
    for (char base : bases) {
      string temp_base = temp_contig + base;
      if(kmer_map.find(temp_base) != kmer_map.end()){
        if(kmer_map[temp_base]>avr_kmer_dep*3){
          kmer_map[temp_base]=avr_kmer_dep*3;
        }
        else {
        kmer_map[temp_base] = 0;
        }
      } 
    }
  }

  return contig;
}