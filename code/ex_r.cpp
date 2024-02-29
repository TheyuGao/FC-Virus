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

int get_seed(int kmer_length) {
  int max_read = 0;
  int max_Read_dep = 0;
  int temp_dep = 0;
  for (int i = 0; i < data_str.size(); i++) {
    temp_dep = 0;
    for (int j = 0; j <= data_str[i].length() - kmer_length; j++) {
      const string K_mer = data_str[i].substr(j, kmer_length);
      if (kmer_map[K_mer] >= first_p) {

         temp_dep += kmer_map[K_mer];
      }
    }

    if (temp_dep > max_Read_dep) {
      max_Read_dep = temp_dep;
      max_read = i;
    }
  }
  cout << data_str[max_read] << "    is seed read" << endl;
  return max_read;
}

string construct_contig(int kmer_length, vector<int> kmer_position,
                        int max_read, string &first_kmer2, string &second_kmer2,
                        vector<bool> &had_used_read) {
  map<string, int> insert_map;
  int break_count;
  string temp_contig = "";
  int ii = 0;
  int jj = 1;

  while (ii < kmer_position.size()-1 && jj < kmer_position.size()) {

    string first_kmer =
        data_str[max_read].substr(kmer_position[ii], kmer_length);
    string second_kmer =
        data_str[max_read].substr(kmer_position[jj], kmer_length);

    vector<pair<int, int>> first_kmer_read = k_r_p[first_kmer];
    int break_down = 0;
    for (int kk = 0; kk < first_kmer_read.size(); kk++) {
      break_down += data_dep[first_kmer_read[kk].first];
      if (break_count >= break_const) {
        break;
      }
    }


    break_down = break_down * break_ratio;

    break_count = 0;
    insert_map.clear();
    for (int k = 0; k < first_kmer_read.size(); k++) {
      int search_read = first_kmer_read[k].first;
      int search_read_p = first_kmer_read[k].second;
      std::string::size_type second_kmer_p = string::npos;
      if (search_read_p + kmer_length < data_str[search_read].length()) {
        second_kmer_p =
            data_str[search_read].substr(search_read_p + 1).find(second_kmer);
        if (second_kmer_p != string::npos) {
          break_count += data_dep[search_read];
          string insert_str =
              data_str[search_read].substr(search_read_p, second_kmer_p + 1);
          if (insert_map.count(insert_str) == 0)
            insert_map.insert({insert_str, 1});
          else
            insert_map[insert_str] += data_dep[search_read];
          if (break_count > break_down) {
            break;
          }
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
    temp_contig += max_str;
    first_kmer2 = first_kmer;
    second_kmer2 = second_kmer;
    ii++; 
    jj++;
  }
  return temp_contig;
}

string Extend_R(int kmer_length, int seed_read, vector<bool> &had_used_read) {
  string test_contig = "";
  int max_read = seed_read;
  string first_kmer, second_kmer;
  string last_first_kmer;
  string last_second_kmer;
  bool ex_finish = false;
  vector<int> kmer_position;
  int star_position = 0;
  int last_max_read = 0;

  string contig = "";

  while (max_read != -1) {
    last_max_read = max_read;
    kmer_position.clear();
    for (int j = star_position; j <= data_str[max_read].length() - kmer_length;
         j++) {
      const string &kmer = data_str[max_read].substr(j, kmer_length);
      if (kmer_map[kmer] >= first_p && kmer_map[kmer] <= second_p) {
        kmer_position.push_back(j);
      }
    }


    int test_lengtha=test_contig.length();

    string candi_contig =
        construct_contig(kmer_length, kmer_position, max_read, first_kmer,
                         second_kmer, had_used_read);

    
    cout<<"Constructing......"<<endl;
    test_contig += candi_contig;

    int temp_length = test_contig.length() - test_lengtha;

    max_read = -1;
    vector<pair<int, int>> first_kmer_read = k_r_p[first_kmer];

        for (int i = 0; i < first_kmer_read.size(); i++) {
      if (had_used_read[first_kmer_read[i].first] == true)
        continue;
      int search_read = k_r_p[first_kmer][i].first;
      int search_read_p = k_r_p[first_kmer][i].second;
      int start_p;
      for (int j = search_read_p;
           j < data_str[search_read].length() - kmer_length; j++) {
        string kmer = data_str[search_read].substr(j, kmer_length);
        if (kmer == second_kmer) {
          start_p = j;
          int kmer_num = 0;
          for (int k = j + 1; k <= data_str[search_read].length() - kmer_length;
               k++) {
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

  test_contig = test_contig + second_kmer;

  return test_contig;
}
string final_contig(int kmer_length, string contig) {
  int min_fre = 0;
  for (int i = 0; i <= contig.length() - kmer_length; i++) {
    string k_mer = contig.substr(i, kmer_length - 1);
    if (kmer_map.find(k_mer + "A") != kmer_map.end()) {
      kmer_map[k_mer + "A"] =0;

    }
    
  }


  if (dep_search==2) {
    min_fre=30;
  }
  else {
    min_fre=2;
  }



  while (true) {

    string temp_contig =
        contig.substr(contig.length() - kmer_length + 1, kmer_length - 1);
    string bases = "ATGC";
    char max_base = 'A';
    int max_cov = 0;

    for (char base : bases) {
      string temp_base = temp_contig + base;
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

    contig += max_base;

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