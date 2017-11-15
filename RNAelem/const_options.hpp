//
//  const_options.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/12/11.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef const_options_h
#define const_options_h

namespace iyak {
  enum {
    DBG_NONE = 0,
    DBG_NO_THETA = 1<<0,
    DBG_NO_ENE = 1<<1, // obsolete
    DBG_FIX_RSS = 1<<2,
    DBG_VERBOSE = 1<<3,
    DBG_NO_TURN = 1<<4,
    DBG_PROOF = 1<<5,
    DBG_ARRAY = 1<<6,
    DBG_KEEP_TMP = 1<<7,
    DBG_CORE_FILE = 1<<8,
    DBG_NO_LOGSUM = 1<<9,
  };
  enum {
    debug = DBG_NONE
  };
#define DBG_NO_MULTI 1
#define MAX_SEQLEN 9999
}

#endif /* const_options_h */
