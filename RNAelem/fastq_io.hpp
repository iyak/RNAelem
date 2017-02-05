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
#include<string>
#include"util.hpp"
#include"bio_sequence.hpp"

namespace iyak {

  using std::getline;

  class FastqReader {

    string _fname;
    string _encoding;
    int _base;
    int _cnt;

    int _N;

    ifstream _ifs;

    string _id;
    string _seq_s;
    string _qual_s;
    string _rss;

    void store_one_record() {
      string dummy;
      getline(_ifs, _id, '\n');
      getline(_ifs, _seq_s, '\n');
      getline(_ifs, dummy, '\n');
      getline(_ifs, _qual_s, '\n');
      if (debug&DBG_FIX_RSS)
        getline(_ifs, _rss, '\n');
    }

    public:

    void static qual_stoi(string const& s, VI& qual, int base=33) {
      qual.resize(s.size());
      for (int i=0; i < size(s); ++ i) {
        qual[i] = int(s[i]) - base;
      }
    }

    /* getter */
    int const& N() {return _N;}

    /* setter */
    void set_fq_fname(string fname, string encoding = "sanger") {
      _fname = fname;
      _encoding = encoding;

      _base =
      "sanger" == encoding? 33:
      "solexa" == encoding? 64:
      "illumina1.3" == encoding? 64:
      "illumina1.5" == encoding? 64:
      "illumina1.8" == encoding? 33:
      -1;
      check(-1 != _base, "wrong encoding:", encoding);

      if (_ifs.is_open()) _ifs.close();
      _ifs.open(fname);
      check(!!_ifs, "could not open:", fname);

      clear();

      /* count the number of reads */
      string id;
      VI seq;
      VI qual;
      string rss;
      for (_N=0; not is_end(); ++_N) {read_seq(id, seq, qual, rss);}

      clear();
     };

    ~FastqReader() {
      if (_ifs.is_open()) _ifs.close();
    }

    void read_seq(string& id, VI& seq, VI& qual, string& rss) {
      id = _id;
      seq_stoi(_seq_s, seq);
      qual_stoi(_qual_s, qual, _base);
      if (debug&DBG_FIX_RSS) {
        rss = _rss;
        check(size(seq) == size(rss), "mismatch len:", seq, rss);
      }

      store_one_record();
      ++ _cnt;
    }

    void clear() {
      _ifs.clear();
      _ifs.seekg(0);
      store_one_record();
      _cnt = 0;
    }

    bool is_end() {
      return !(_ifs.is_open()) or _ifs.eof();
    }

    int cnt() {return _cnt;}
  };
}

#endif /* fastq_io_h */
