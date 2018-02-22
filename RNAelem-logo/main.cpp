//
//  main.cpp
//  RNAelem-logo
//
//  Created by Hiroshi Miyake on 2018/02/19.
//  Copyright © 2018 Hiroshi Miyake. All rights reserved.
//
#include<iostream>
#include<type_traits>
#include"util.hpp"
#include"logo.hpp"
using namespace iyak;

template<class T>vector<vector<T>> parse_table(string const& s){
  vector<vector<T>> x{};
  VI stack{};
  int j0=(int)s.find_first_of("[");
  int j1=(int)s.find_last_of("]");
  for(int j=j0+1;j<j1;++j){
    if('['==s[j])
      stack.push_back(j);
    else if(']'==s[j]){
      int i=stack.back();
      x.push_back(split<T>(s.substr(i+1,j-i-1),","));
      stack.pop_back();
    }
  }
  return x;
}

void foreach(string const& s,std::function<void(string const&,int)>callback){
  int j0=(int)s.find_first_of("["),j1=(int)s.find_last_of("]");
  if(npos==j0 and npos==j1){
    callback(s,0);
    return;
  }
  int from=j0+1,i=0;
  for(int j=from;j<j1;++j){
    if('['==s[j]){
      ++j;
      for(int nest=1;0<nest and j<j1;++j)nest+=('['==s[j])?1:(']'==s[j])?-1:0;
      --j;
    }else if(','==s[j]){
      callback(s.substr(from,j-from),i);
      from=j+1;
      ++i;
    }
  }
  callback(s.substr(from,j1-from),i);
}

void parse_input(RNAlogo& logo,istream& fi){
  string seq="",val="",color="",meta="";
  while(not fi.eof()){
    string s;
    getline(fi,s);
    auto p=split<string>(s,":");
    if(size(p)<2)continue;
    check(2==p.size(),"fail to parse");
    if("seq"==strip(p[0]))seq=strip(p[1]);
    else if("val"==strip(p[0]))val=strip(p[1]);
    else if("color"==strip(p[0]))color=strip(p[1]);
    else if("meta"==strip(p[0]))meta=strip(p[1]);
    else;
  }
  foreach(seq,[&](string const& ss,int){
    logo.data().columns.push_back({});
    foreach(ss,[&](string const& sss,int){
      logo.data().columns.back().blocks.emplace_back(sss);
    });
  });
  foreach(val,[&](string const& ss,int i){
    foreach(ss,[&](string const& sss,int j){
      logo.data().columns[i].blocks[j].val=iss_cast<double>(sss);
    });
  });
  foreach(color,[&](string const& ss,int i){
    foreach(ss,[&](string const& sss,int j){
      foreach(sss,[&](string const& ssss,int k){
        if(size(strip(ssss)))
          logo.data().columns[i].blocks[j].colors[k]=ssss;
      });
    });
  });
  foreach(meta,[&](string const& ss,int i){
    logo.data().columns[i].meta=ss;
  });
}

int main(int argc,char *argv[]){
  RNAlogo logo;
  logo.set_ostream(std::cout);
  try{
    if(1==argc){
      parse_input(logo,std::cin);
    }else{
      ifstream ifs(argv[1]);
      check(!!ifs,"cannot open:",argv[1]);
      parse_input(logo,ifs);
    }
    logo.plot();
  }catch(...){
    //logo.no_logo();
  }
  return 0;
}
