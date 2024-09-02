
#include "kmer.h"
#include "GeneralSet.h"
#include "loadreads.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <math.h>
#include <sstream>
#include <stdlib.h>

using namespace std;

bool contains_non_gatc(string kmer) {

  for (unsigned int i = 0; i < kmer.size(); ++i) {

    unsigned char c = kmer[i];

    if (c != 'A' && c != 'T' && c != 'G' && c != 'C' && c != 'a' && c != 't' &&
        c != 'g' && c != 'c') {
      return (false);
    }
  }

  return (true);
}

bool Find_peak(int kmer_length) {

  int max_kmerfre = 0;
  int dep_judge = 0;

  for (map<string, int>::iterator it = kmer_map.begin(); it != kmer_map.end();
       it++) {
    if (max_kmerfre < it->second) {
      max_kmerfre = it->second;
    }
  }
  cout << max_kmerfre << "is max_kmer_fre" << endl;

  int count_kmer[max_kmerfre] = {0};

  for (map<string, int>::iterator it = kmer_map.begin(); it != kmer_map.end();
       it++) {

    count_kmer[it->second]++;
  }

  int U = 0;
  int Seq_Dep = 0;

  int start = max_kmerfre * g_ratio;
  int end = max_kmerfre * (1 - g_ratio);
  if (start < 2) {
    start = 2;
  }
  if (end == max_kmerfre) {
    end -= 1;
  }
  for (int i = start; i <= end; i++) {
    if (count_kmer[i] > 0) {
      Seq_Dep += count_kmer[i];
      U++;
    }
  }
  Seq_Dep = Seq_Dep / U;
  cout << "Seq dep is" << Seq_Dep << endl;

  int peak[max_kmerfre] = {0};

  int step_length = max_kmerfre / g_window_length;
  if (max_kmerfre < g_window_length) {
    step_length = 1;
  }
  int win_length = step_length * 10;

  int mark = 0;
  for (int i = 0; i < max_kmerfre - win_length; i += step_length) {
    int none_zero = 0;
    for (int j = 0; j < (win_length); j++) {

      if (count_kmer[i + j] != 0) {
        none_zero++;
        peak[i] += count_kmer[i + j];
      }
    }
    if (none_zero > 0) {
      peak[i] = peak[i] / none_zero;
      mark = i;
    }
  }

  cout<<"1"<<endl;

  for (int i = mark; i >= 0; i -= step_length) {
    dep_judge = 0;
    int j = i - step_length;

    while (j >= 0 && peak[j] >= peak[j + step_length]) {
      j -= step_length;
    }
    if (j > 0) {
      j -= step_length;
      while (j >= 0 && peak[j] >= peak[j + step_length]) {
        j -= step_length;
      }
      while (j >= Seq_Dep && peak[j] <= peak[j + step_length]) {
        j -= step_length;
      }
      if (j >= Seq_Dep) {
        j -= step_length;
        while (j >= Seq_Dep && peak[j] <= peak[j + step_length]) {
          j -= step_length;
        }
      }
      cout<<"2"<<endl;


      //debug
      if(2*j-max_kmerfre>0){
        j=j-(max_kmerfre-j);
      }

      int judge = 0;
      int kmer_dep = Seq_Dep * U;

      // 1
      first_p = j;
      second_p = max_kmerfre;
      dep_search = 0;
      int find_peak_right_num = 0;
      // 2

      for (int ii = j; ii <= max_kmerfre-1; ii++) {
        //debug
        if(count_kmer[ii]>500){
          cout<<"异常区间"<<ii<<endl;
          cout<<count_kmer[ii]<<endl;
        }
        //debug



        dep_judge += count_kmer[ii];
        if (count_kmer[ii] > 0) {
          find_peak_right_num++;
        }
      }
      cout<<"3"<<endl;

      for (int ii = 0; ii < data_str.size(); ii++) {
        avr_kmer_dep += (data_str[ii].length() - kmer_length+ 1) * data_dep[ii];
      }
      avr_kmer_dep = avr_kmer_dep / kmer_map.size();
      cout << " avr_kmer_dep" << avr_kmer_dep << endl;

      cout<<"4"<<endl;


     ////这个值需要具体重新算一下

     cout<<"计算最后一个峰区间起始位置"<<j<<endl;
      if (j > avr_kmer_dep*g_re_gene)
        dep_search = 2;



      cout<<"选取到的k-mer数量"<<dep_judge<<endl;
      cout<<"总共数量"<<kmer_dep<<endl;
      while (dep_judge < kmer_dep * g_com_kmer_ratio && first_p > 2 &&
             find_peak_right_num < g_window_length) {
        first_p = first_p - 1;
        dep_judge += count_kmer[first_p];
        if (count_kmer[first_p] > 0) {
          find_peak_right_num++;
        }
      }
      cout<<"最终取得的起始位置"<<first_p<<endl;
      cout<<"5"<<endl;
      // 算平均k-mer厚度，j是否大于平均k-mer厚度的二倍


      cout<<"取得的dep_search状态"<<dep_search<<endl;



      return true;
    }
  }
  return false;
}

void get_kmer(int kmer_length) {

  int data_size = data_str.size();
  cout << "kmer map contructing ..." << endl;
  time_t start_time = time(NULL);

  if (data_str.empty()) {
    cout << "data empty!!" << endl;
    return;
  }

  for (int num = 0; num < data_str.size(); num++) {
    // cout<<"kmer"<<endl;
    const string &read = data_str[num];
    {
      if (read.length() < kmer_length)
        continue;
    }
    for (int j = 0; j <= read.length() - kmer_length; j++) {
      const string &kmer = read.substr(j, kmer_length);

      if (!contains_non_gatc(kmer)) {
        continue;
      }

      if (kmer_map.find(kmer) == kmer_map.end()) {
        kmer_map.insert(make_pair(kmer, data_dep[num]));
      } else {
        kmer_map[kmer] += data_dep[num];
      }
    }
  }
  time_t kmer_map_time = time(NULL);
  cout << " kmer_map" << (kmer_map_time - start_time) << " s)" << endl;

  if (!Find_peak(kmer_length)) {
    cout << "can not find peak!" << endl;
    return;
  }
  time_t find_peak_time = time(NULL);
  cout << " kmer_map" << (find_peak_time - start_time) << " s)" << endl;

  for (int i = 0; i < data_str.size(); i++) {
    const string &read = data_str[i];

    if (read.length() <= kmer_length)
      continue;

    for (int j = 0; j <= read.length() - kmer_length; j++) {
      const string &kmer = read.substr(j, kmer_length);
      if (kmer_map[kmer] >= first_p && kmer_map[kmer] <= second_p) {

        k_r_p[kmer].push_back(make_pair(i, j));
      }
    }
  }
  cout << "k_r_p success found! size is " << k_r_p.size() << endl;
  time_t end_time = time(NULL);
  cout << "Kmer_map has been constructed, total " << kmer_map.size()
       << " kmers! time total: " << (end_time - start_time) << " s)" << endl;
}
