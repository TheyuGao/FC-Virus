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


// string construct_contig_reverse(int kmer_length, vector<int> kmer_position,
//                         int max_read, string &first_kmer2, string &second_kmer2,
//                         vector<bool> &had_used_read) {
//   map<string, int> insert_map;
//   int break_count;
//   string temp_contig = "";
//   int ii = 0;
//   int jj = 1;

//   while (ii < kmer_position.size()-1 && jj < kmer_position.size()) {

//     string first_kmer =
//         data_str[max_read].substr(kmer_position[ii], kmer_length);
//     string second_kmer =
//         data_str[max_read].substr(kmer_position[jj], kmer_length);

//     vector<pair<int, int>> first_kmer_read = k_r_p[first_kmer];
//     int break_down = 0;
//     for (int kk = 0; kk < first_kmer_read.size(); kk++) {
//       break_down += data_dep[first_kmer_read[kk].first];
//       if (break_count >= break_const) {
//         break;
//       }
//     }

//     // 改成全局变量的输入，还有1666
//     break_down = break_down * break_ratio;

//     break_count = 0;
//     insert_map.clear();
//     for (int k = 0; k < first_kmer_read.size(); k++) {
//       int search_read = first_kmer_read[k].first;
//       int search_read_p = first_kmer_read[k].second;
//       std::string::size_type second_kmer_p = string::npos;
//       if (search_read>0) {
//         second_kmer_p =
//             data_str[search_read].find(second_kmer);
//         if (second_kmer_p != string::npos) {
//           break_count += data_dep[search_read];
//           string insert_str =
//               data_str[search_read].substr(second_kmer_p,search_read_p-search_read_p);
//           if (insert_map.count(insert_str) == 0)
//             insert_map.insert({insert_str, 1});
//           else
//             insert_map[insert_str] += data_dep[search_read];
//           if (break_count > break_down) {
//             break;
//           }
//         }
//       }
//     }

//     int max_num = 0;
//     string max_str = "";
//     for (map<string, int>::iterator kmer_it = insert_map.begin();
//          kmer_it != insert_map.end(); kmer_it++) {

//       if (kmer_it->second > max_num) {
//         max_num = kmer_it->second;
//         max_str = kmer_it->first;
//       }
//     }
//     temp_contig = max_str + temp_contig;
//     first_kmer2 = first_kmer;
//     second_kmer2 = second_kmer;
//     ii++; 
//     jj++;
//   }
//   return temp_contig;
// }
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
    cout<<test_contig.length()<<endl;
    test_contig = candi_contig+test_contig;

    int temp_length = test_contig.length() - test_lengtha;

    max_read = -1;
    //改成和ex_r一样的新代码
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
// 需要删除参数
string final_contig_Left(int kmer_length, string contig) {
  int min_fre = 2;
  for (int i = 0; i <= contig.length() - kmer_length; i++) {
    string k_mer = contig.substr(i, kmer_length - 1);
    if (kmer_map.find("A"+k_mer ) != kmer_map.end()) {//for循环，contig相同或有一个差异的k-mer的厚度之和/(contig.length()-kmer_length+1),这个值作为平均厚度。
      kmer_map["A"+k_mer] =0;
 // 有重复的情况，如果kmer_map[kmer+A]>平均k-mer厚度*3（这个值需要确定），不减去，直接等于平均厚度，如果小的话直接赋值为0.TGC同理，需要具体情况分析分析。
    }
    
  }
//   for(int i=0;i<contig.length()-kmer_length;i++){
//     string k_min_mer=contig.substr(i,kmer_length-1);
//     string bases = "ATGC";
//     char max_base = 'A';
//     for (char base : bases) {
//       string temp_base = k_min_mer + base;
//       for(int j=0;j<data_str.size();j++){
//         if (data_str[j].find(temp_base)) {
//           min_fre+=data_dep[j];
//         }
//       }
//     }
//   }
//   cout<<"测试厚度之和为"<<min_fre<<endl;
//   min_fre=min_fre/contig.length();
//   cout<<"最终算得min_fre"<<min_fre<<endl;


  while (true) {
    // 获取当前 contig 的最后 kmer_length - 1 个字符
    string temp_contig =
        contig.substr(0, kmer_length - 1);
    string bases = "ATGC";
    char max_base = 'A';
    int max_cov = 0;

    // 遍历所有碱基，找到最大覆盖度及对应碱基
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

    // 如果最大覆盖度不大于最小频率阈值，则终止循环
    if (max_cov <= min_fre) {
      break;
    }

    // 将最大覆盖度对应的碱基追加到 contig 中，并更新 temp_contig
    contig = max_base+contig;

    // 将与 temp_contig 相关的所有 kmer 的频率设为 0
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