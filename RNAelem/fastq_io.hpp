//
//  fastq_io.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/14.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef fastq_io_h
#define fastq_io_h

#include<fstream>
#include<iostream>
#include<vector>
#include<random>
#include<string>
#include"util.hpp"
#include"bio_sequence.hpp"

namespace iyak {
  using std::getline;

  class FastqReader{
    string _fname;
    string _encoding;
    int _base;
    int _cnt=0;
    int _cnt_shf=0;
    string _seq_cat;
    string _qual_cat;
    string _id_cat;
    VI _i_seq,_l_seq;
    VI _i_qual,_l_qual;
    VI _i_id,_l_id;
    int _N;
    void static qual_int2str(VI const& qual,string& s,int base=33){
      s.resize(qual.size());
      for(int i=0;i<size(qual);++i)
        s[i]=char(base+qual[i]);
    }
    void static qual_str2int(string const& s,VI& qual,int base=33){
      qual.resize(s.size());
      for(int i=0;i<size(s);++i)
        qual[i]=int(s[i])-base;
    }
    void seq_region(int from,int len,VI& dest){
      dest.resize(len);
      for(int i=0;i<len;++i){
        char c=_seq_cat[from+i];
        dest[i]=
        'A'==c or 'a'==c?1:
        'C'==c or 'c'==c?2:
        'G'==c or 'g'==c?3:
        'U'==c or 'u'==c or 'T'==c or 't'==c?4:
        0;
      }
    }
    void qual_region(int from,int len,VI& dest){
      dest.resize(len);
      for(int i=0;i<len;++i)
        dest[i]=int(_qual_cat[i])-_base;
    }
  public:
    void set_fq_fname(string const& fname,string const& encoding="sanger"){
      _fname=fname;
      _encoding=encoding;
      _base=
      "sanger"==encoding?33:
      "solexa"==encoding?64:
      "illumina1.3"==encoding?64:
      "illumina1.5"==encoding?64:
      "illumina1.8"==encoding?33:
      -1;
      check(-1!=_base,"wrong encoding:",encoding);
      ifstream ifs;
      ifs.open(fname);
      check(!!ifs,"could not open:",fname);
      string id,seq,dummy,qual;
      _i_seq.clear();
      _l_seq.clear();
      _i_qual.clear();
      _l_qual.clear();
      _i_id.clear();
      _l_id.clear();
      for(int i=0,j=0,k=0;not ifs.eof();){
        getline(ifs,id,'\n');
        getline(ifs,seq,'\n');
        getline(ifs,dummy,'\n');
        getline(ifs,qual,'\n');
        if(ifs.eof())
          break;
        _id_cat.append(id);
        _seq_cat.append(seq);
        _qual_cat.append(qual);
        _i_seq.push_back(i);
        _l_seq.push_back(size(seq));
        _i_qual.push_back(j);
        _l_qual.push_back(size(qual));
        _i_id.push_back(k);
        _l_id.push_back(size(id));
        i+=size(seq);
        j+=size(qual);
        k+=size(id);
      }
      _N=size(_i_seq);
      _cnt=0;
      _cnt_shf=0;
    }
    void get_read(string& id,VI& seq,VI& qual,string&/*rss*/){
      id=_id_cat.substr(_i_id[_cnt],_l_id[_cnt]);
      seq_region(_i_seq[_cnt],_l_seq[_cnt],seq);
      qual_region(_i_qual[_cnt],_l_qual[_cnt],qual);
      ++_cnt;
    }
    void shuffle(){
      std::mt19937 m;
      m.seed(_cnt_shf);std::shuffle(_i_id.begin(),_i_id.end(),m);
      m.seed(_cnt_shf);std::shuffle(_l_id.begin(),_l_id.end(),m);
      m.seed(_cnt_shf);std::shuffle(_i_seq.begin(),_i_seq.end(),m);
      m.seed(_cnt_shf);std::shuffle(_l_seq.begin(),_l_seq.end(),m);
      m.seed(_cnt_shf);std::shuffle(_i_qual.begin(),_i_qual.end(),m);
      m.seed(_cnt_shf);std::shuffle(_l_qual.begin(),_l_qual.end(),m);
      ++_cnt_shf;
    }
    bool is_end(){return _cnt==_N;}
    void skip(){++_cnt;}
    void clear(){_cnt=0;}
    int cnt(){return _cnt;}
    int N(){return _N;}
  };

  class FastqBatchReader{
    FastqReader _qr;
    int _N_batch;
    int _cnt=0;
    int _cnt_epoc=0;
  public:
    void set_fq_name(string const& fname, string const& encoding="sanger"){
      _qr.set_fq_fname(fname,encoding);
      _cnt=0;
      _cnt_epoc=0;
    }
    void set_batch_size(int N){
      if(N<0) _N_batch=_qr.N();
      else _N_batch=N;
    }
    void get_read(string& id,VI& seq,VI& qual,string& rss){
      _qr.get_read(id,seq,qual,rss);
      ++_cnt;
    }
    bool is_end(){return _N_batch<=_cnt or _qr.is_end();}
    bool is_end_epoc(){return _qr.is_end();}
    void clear(){
      if(is_end_epoc()){
        _qr.shuffle();
        _qr.clear();
        ++_cnt_epoc;
      }
      _cnt=0;
    }
    void skip(){++_cnt;_qr.skip();}
    int cnt(){return _cnt;}
    int cnt_epoc(){return _cnt_epoc;}
    int N(){return _qr.N();}
  };
}

#endif /* fastq_io_h */
