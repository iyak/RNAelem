//
//  bio_sequence.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/09.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef bio_sequence_h
#define bio_sequence_h

#include<string>
#include"util.hpp"

namespace iyak {

  enum {nchar = 5, nchar2 = 7,};
  string const nacgu = "NACGU";

  int const bp[nchar][nchar]
  /*         {N,A,C,G,U} */
    ={/* N */{0,0,0,0,0},
      /* A */{0,0,0,0,5},
      /* C */{0,0,0,1,0},
      /* G */{0,0,2,0,3},
      /* U */{0,6,0,4,0}};

  void seq_str2int(string const& s, VI& seq) {
    seq.resize(size(s));
    for (int i=0; i < size(s); ++ i) {
      switch (s[i]) {
        case 'A': case 'a': seq[i] = 1; break;
        case 'C': case 'c': seq[i] = 2; break;
        case 'G': case 'g': seq[i] = 3; break;
        case 'T': case 't': case 'U': case 'u': seq[i] = 4; break;
        default: seq[i] = 0;
      };
    }
  }

  void seq_int2str(VI const& seq, string& s) {
    s.resize(size(seq));
    for (int i=0; i < size(seq); ++ i) {
      s[i] = nacgu[seq[i]];
    }
  }

  string seq_int2str(VI const& seq) {
    string s = "";
    for (int i=0; i < size(seq); ++ i) {
      s += nacgu[seq[i]];
    }
    return s;
  }

  string slice(VI const& s, int i, int j) {
    string t = "";
    for (int k=i; k < j; ++ k) {
      t += nacgu[s[k]];
    }
    return t;
  }

}

#endif /* bio_sequence_h */
