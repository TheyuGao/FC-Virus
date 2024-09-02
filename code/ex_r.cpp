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

// 重新试一下直接使用两个以上公共k-mer的read对结果是否具有影响？
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
        // temp_dep+=1;
      }
    }
    // if(temp_dep>10){
    //   max_read=i;
    //   cout << max_read << "is max read" << endl;
    //   return max_read;
    // }
    if (temp_dep > max_Read_dep) {
      max_Read_dep = temp_dep;
      max_read = i;
    }
  }
  cout << max_read << "is max read" << endl;
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

    // 改成全局变量的输入，还有1666
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
  //  if(break_count<break_down) {
  //   vector<pair<int, int>> first_kmer_read = k_r_p[first_kmer];
  //   int had_second_kmer_read=0;
  //   for (int h=0; h<first_kmer_read.size(); h++) {
  //   if (data_str[first_kmer_read[h].first].find(second_kmer)) {
  //     had_second_kmer_read++;
  //   }
  //   }
  //   if(had_used_read*4 < first_kmer_read.size())
  //   jj+1;continue;
  //   }

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

    
    cout<<test_contig.length()<<endl;
    test_contig += candi_contig;

    int temp_length = test_contig.length() - test_lengtha;

    max_read = -1;
    vector<pair<int, int>> first_kmer_read = k_r_p[first_kmer];

//     for (int i = 0; i < first_kmer_read.size(); i++) {
//       //最大的reads标注，所有用过的都进行标注
//       if (had_used_read[first_kmer_read[i].first] == true)
//         continue;
//       int search_read = k_r_p[first_kmer][i].first;
//       int search_read_p = k_r_p[first_kmer][i].second;
//       int start_p;
//       int second_kmer_p;
//       if (search_read_p + kmer_length < data_str[search_read].length()) {
//         second_kmer_p =
//             data_str[search_read].substr(search_read_p + 1).find(second_kmer);
//         if (second_kmer_p != string::npos) {
// //debug一下是否生效，没有生效的话输出second k-mer
//           for (int k = second_kmer_p+search_read_p  + 1;
//                k <= data_str[search_read].length() - kmer_length; k++) {
//             if (kmer_map[data_str[first_kmer_read[i].first].substr(
//                     k, kmer_length)] >= first_p &&
//                 kmer_map[data_str[first_kmer_read[i].first].substr(
//                     k, kmer_length)] <= second_p) {
//               star_position = k;
//               max_read = first_kmer_read[i].first;
//               break;
//             }
//           }
//           if (max_read != -1) {
//             break;
//           }
//         }
//       }
//     }
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
    if (kmer_map.find(k_mer + "A") != kmer_map.end()) {//for循环，contig相同或有一个差异的k-mer的厚度之和/(contig.length()-kmer_length+1),这个值作为平均厚度。
      kmer_map[k_mer + "A"] =0;
 // 有重复的情况，如果kmer_map[kmer+A]>平均k-mer厚度*3（这个值需要确定），不减去，直接等于平均厚度，如果小的话直接赋值为0.TGC同理，需要具体情况分析分析。
    }
    
  }
  // for(int i=0;i<contig.length()-kmer_length;i++){
  //   string k_min_mer=contig.substr(i,kmer_length-1);
  //   string bases = "ATGC";
  //   char max_base = 'A';
  //   for (char base : bases) {
  //     string temp_base = k_min_mer + base;
  //     for(int j=0;j<data_str.size();j++){
  //       if (data_str[j].find(temp_base)) {
  //         min_fre+=data_dep[j];
  //       }
  //     }
  //   }
  // }
  // cout<<"测试厚度之和为"<<min_fre<<endl;
  if (dep_search==2) {
    min_fre=30;
  }
  else {
    min_fre=2;
  }

  while (true) {
    // 获取当前 contig 的最后 kmer_length - 1 个字符
    string temp_contig =
        contig.substr(contig.length() - kmer_length + 1, kmer_length - 1);
    string bases = "ATGC";
    char max_base = 'A';
    int max_cov = 0;

    // 遍历所有碱基，找到最大覆盖度及对应碱基
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

    // 如果最大覆盖度不大于最小频率阈值，则终止循环
    if (max_cov <= min_fre) {
      break;
    }

    // 将最大覆盖度对应的碱基追加到 contig 中，并更新 temp_contig
    contig += max_base;

    // 将与 temp_contig 相关的所有 kmer 的频率设为 0
    for (char base : bases) {
      string temp_base = temp_contig + base;
      if(kmer_map.find(temp_base) != kmer_map.end()){
        if(kmer_map[temp_base]>avr_kmer_dep*100){
          kmer_map[temp_base]=avr_kmer_dep*100;
        }
        else {
        kmer_map[temp_base] = 0;
        }
      } 
    }
  }

  return contig;
}
