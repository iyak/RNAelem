//
//  arrayjob_manager.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2017/01/21.
//  Copyright © 2017 Kiryu Lab. All rights reserved.
//

#ifndef arrayjob_manager_h
#define arrayjob_manager_h

#include<cstdlib>
#include<array>
#include"util.hpp"

namespace iyak {

  class ArrayJobManager {

    string _command;
    string _env_tid;
    string _array;
    string _binary;
    string _sync;
    string _cwd;
    string _env;
    string _other;
    unsigned int _set = 0;

    bool all_set() {return _set==(1<<8)-1;}

    void read_sge_option_stream(istream& is) {
      _set = 0;
      string l;
      while (getline(is, l)) {

        auto vs = split<string>(strip(l), ":");
        if (2 != size(vs)) continue;
        auto key = strip(vs[0]);
        auto val = strip(vs[1]);

        if ("command"==key) {
          _command = val;
          _set |= 1<<0;
        }
        else if ("task id"==key) {
          _env_tid = val;
          _set |= 1<<1;
        }
        else if ("array"==key) {
          _array = val;
          _set |= 1<<2;
        }
        else if ("binary"==key) {
          _binary = val;
          _set |= 1<<3;
        }
        else if ("sync"==key) {
          _sync = val;
          _set |= 1<<4;
        }
        else if ("cwd"==key) {
          _cwd = val;
          _set |= 1<<5;
        }
        else if ("environment"==key) {
          _env = val;
          _set |= 1<<6;
        }
        else if ("other"==key) {
          _other = val;
          _set |= 1<<7;
        }
        else {
          cry("not used:", key);
        }
      }
      check(all_set(), "grid_engine_opt broken.");
    }


    string exec(string const& cmd) {
      std::array<char, 512> buf;
      string res;
      sptr<FILE> pipe(popen((cmd + " 2>&1").c_str(), "r"), pclose);
      check(!!pipe, "fail popen shell");

      while(!feof(pipe.get()))
        if (fgets(buf.data(), 512, pipe.get()))
          res += buf.data();

      return res;
    }

  public:
    void set_sge_option(string const& fname) {
      if ("~DEFAULT~" == fname) {
        isstream iss (
#include"grid_engine_opt"
        );
        read_sge_option_stream(iss);
      }
      else {
        ifstream ifs (fname);
        check(!!ifs, "fail opne:", fname);
        read_sge_option_stream(ifs);
      }
    }

    int tid() {
      check(all_set(), "called tid() with options unset");
      if (debug&DBG_ARRAY) {
        return 5;
      } else {
        char const* tid_s = std::getenv(_env_tid.c_str());
        check(!!tid_s, "fail read env:", _env_tid);
        return iss_cast<int>(string(tid_s));
      }
    }

    void submit_array_job(string const& job,
                          int const n,
                          bool const show=false) {
      check(all_set());

      size_t i;

      string array(_array);
      while (npos != (i=array.find("$from"))) array.replace(i,5,"1");
      while (npos != (i=array.find("$to"))) array.replace(i,3,to_str(n));

      string total = paste1(
        _command,
        array,
        _binary,
        _sync,
        _cwd,
        _env,
        _other,
        "\"" + job + "\""
      );

      if (show) cry("submit:", total);
      string res = exec(total);
      if (show) cry(strip(res, "\n"));
    }

    PII assinged_range(int total, int n, int k) {
      int res = total - (total/n) * n;
      int from=0, to=0;
      for (int i=0; i<k; ++i) {
        from = to;
        to += total/n + (i<res? 1:0);
      }
      return {from, to};
    }
  };
}

#endif /* arrayjob_manager_h */
