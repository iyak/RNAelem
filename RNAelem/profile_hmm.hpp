//
//  profile_hmm.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2016/11/15.
//  Copyright © 2016 Kiryu Lab. All rights reserved.
//

#ifndef profile_hmm_h
#define profile_hmm_h

#include"util.hpp"
#include"bio_sequence.hpp"

namespace iyak {

  /* hidden state */
  struct intervalState {
    int id;
    int l;
    int r;
  };
  using IS = intervalState;
  bool operator== (IS const& a, IS const& b) {
    return a.id==b.id;
  }

  class ProfileHMM {

  public:

    using VIS = vector<IS>;
    using VVIS = vector<VIS>;

  private:

    int M;

    string _pattern;
    VI _node;
    VI _pair;
    VVI _edge_to;
    VVI _edge_from;
    VI _theta_id;
    VV _theta;
    VV _s;
    VVB _reachable;
    VVB _reachable_as_loop;
    VIS _state;
    VVIS _nodes_to_state;
    VIS _loop_state;
    VVIS _loop_right_trans;
    VVIS _loop_left_trans;
    VVIS _pair_trans;
    VVIS _loop_loop_states;

  public:

    /* getter */
    string const& pattern() {return _pattern;}
    VIS const& state() {return _state;}
    VIS const& loop_state() {return _loop_state;}
    VIS const& loop_right_trans(int s) {return _loop_right_trans[s];}
    VIS const& loop_left_trans(int s) {return _loop_left_trans[s];}
    VIS const& pair_trans(int s) {return _pair_trans[s];}
    VVIS const& loop_loop_states() {return _loop_loop_states;}
    int theta_id(int i) {return _theta_id[i];}

    int node(int i) {return _node[i];}
    bool reachable(int s, int s1) {return _reachable[s][s1];}

    void clear_emit_count(VV& e) {
      e.resize(theta().size(), V{});
      for (int i=0; i<iyak::size(theta()); ++i) {
        e[i].resize(theta()[i].size(), 0);
        for (auto& eij: e[i]) eij = 0;
      }
    }

    IS n2s(int h, int h1) {

      if (debug&DBG_PROOF)
        check(h < int(_nodes_to_state.size()) and
              h1 < int(_nodes_to_state[h].size()),
              "bad args to nodes to state:",h,h1);

      IS& is = (debug&DBG_PROOF)?
      _nodes_to_state.at(h).at(h1):
      _nodes_to_state[h][h1];
      if (debug&DBG_PROOF)
        check(is.id >= 0, "nodes to state failed:", h, h1);
      return is;
    }

    size_t size() {return _node.size();}
    VV& theta() {return _theta;}
    VV& s() {return _s;}
    void calc_theta(){
      _theta.clear();
      for(int i=0;i<iyak::size(_s);++i){
        _theta.push_back(V(iyak::size(_s[i]),-1));
        double tot=logsumexp(_s[i]);
        for(int j=0;j<iyak::size(_s[i]);++j)
          _theta[i][j]=_s[i][j]-tot;
      }
    }

    double theta(int const h, int const h1, int const i, int const j) {
      char const cl = _node[h];
      char const cr = _node[h1];
      if (')'==cr) {
        if (debug&DBG_PROOF)
          check(h==_pair[h1], h, _pair[h1]);
        return 0==bp[i][j] or debug&DBG_NO_THETA?
          oneL:
          _theta[_theta_id[h1]][bp[i][j]-1];
      }
      if (debug&DBG_PROOF)
        check((('z'==cl and 'z'==cr) or
               ('z'==cl and '*'==cr) or
               ('z'==cl and 'o'==cr) or
               ('*'==cl and '*'==cr) or
               ('*'==cl and 'o'==cr) or
               ('o'==cl and 'o'==cr)),
              "theta", cl, cr);
      return debug&DBG_NO_THETA?
        oneL:
        (0==i?oneL:_theta[_theta_id[h]][i-1])+
        (0==j?oneL:_theta[_theta_id[h1]][j-1]);
    }

    double theta(int const h, int const j) {
      return 0==j or debug&DBG_NO_THETA?
        oneL:
        _theta[_theta_id[h]][j-1];
    }

    /* setter */
    void add_emit_count(VV& e, int const h, int const h1,
                        int const i, int const j, double const w) {
      char const cl = _node[h];
      char const cr = _node[h1];

      if (')'==cr) {

        if (debug&DBG_PROOF)
          check(h==_pair[h1], h, _pair[h1]);
        if(0<bp[i][j]){
            double &c = e[_theta_id[h1]][bp[i][j]-1];
            c += w;
        }

      } else {

        if (debug&DBG_PROOF)
          check((('z'==cl and 'z'==cr) or
                 ('z'==cl and '*'==cr) or
                 ('z'==cl and 'o'==cr) or
                 ('*'==cl and '*'==cr) or
                 ('*'==cl and 'o'==cr) or
                 ('o'==cl and 'o'==cr)),

                "add_emit_count", cl, cr);

        if(0!=i){
          double &c = e[_theta_id[h]][i-1];
          c += w;
        }
        if(0!=j){
          double &d = e[_theta_id[h1]][j-1];
          d += w;
        }
      }
    }

    void add_emit_count(VV& e, int const h, int const j, double const w) {
      if(0!=j){
        double &c = e[_theta_id[h]][j-1];
        c += w;
      }
    }

    void build(string const& str) {
      _pattern = str;

      check(!str.empty(), "empty motif");

      set_node(str);
      M = (int)_node.size();

      set_pair();
      set_edge();
      set_s_theta();

      set_reachable();
      set_interval_state();
      set_interval_state_trans();

      set_states_pairs();

      save();
    }

    void set_node(string str) {
      _node = VI{'z'};
      _node.insert(_node.end(), str.begin(), str.end());
      _node.insert(_node.end(), {'o'}); 
    }

    /* recognize paired brackets */
    void set_pair() {
      _pair.assign(M, -1);

      VI stack {};
      for (int h=0; h < M; ++h) {
        switch (_node[h]) {
          case '(':
            stack.push_back(h);
            break;
          case ')':
            check(0 < stack.size(), "unmatched brackets");
            int hl = stack.back();
            stack.pop_back();
            _pair[hl] = h;
            _pair[h] = hl;
            break;
        }
      }
      check(0==stack.size(), "unmatched brackets");
    }

    /* define transition between nodes */
    void set_edge() {

      _edge_to.assign(M, VI());
      _edge_from.assign(M, VI());

      for (int h = 0; h<M; ++h) {
        int const& c = _node[h];

        if (0 < h) {
          int const& c1 = _node[h-1];

          if ('*' == c1) { /* no consecutive *s */
            _edge_to[h].push_back(h - 2);
            _edge_from[h - 2].push_back(h);
          }

          _edge_to[h].push_back(h - 1);
          _edge_from[h - 1].push_back(h);
        }

        /* self loop */
        if ('<' != c and '>' != c) {
          _edge_to[h].push_back(h);
          _edge_from[h].push_back(h);
        }
      }
    }

    /* set initial log-ordered profile over bases/base pairs */
    void set_s_theta() {

      _theta_id.assign(_node.size(),-1);
      _s.assign(1,V(nchar-1,0));

      for (int h=0; h < (int)_node.size(); ++h) {
        switch (_node[h]) {
          case ')':
            _theta_id[h]=(int)_s.size();
            _s.push_back(V(nchar2-1,0));
            break;
          case '.':
            _theta_id[h]=(int)_s.size();
            _s.push_back(V(nchar-1,0));
            break;
          case '*':case 'z':case 'o':
            _theta_id[h]=0;
            break;
          case '<':case '>':case '(':
            /* no-op */
            break;
          default:
            die("bad motif char:", _node[h]);
            break;
        };
      }
      calc_theta();
    }

    /* define reachable nodes, which will form IS.*/
    void set_reachable() {

      _reachable.assign(M, VB(M, 0));
      _reachable_as_loop.assign(M, VB(M, 0));

      for (int h = 0; h < M; ++h) {
        int c = _node[h];

        if (c==')') {
          for (int h1: _edge_to[_pair[h]]) {
            _reachable[h1][h] = 1;
          }
        }
        else if (c=='(') {
          // no op
        }
        else if ('>'==c) {
          for (int h1: _edge_to[_pair[h]]) {
            _reachable[h1][h] = 1;
            _reachable_as_loop[h1][h] = 1;
          }
        }
        else if ('<'==c) {
          // no op
        }
        else {
          for (int h1: _edge_to[h]) {
            _reachable[h1][h]         = 1;
            _reachable_as_loop[h1][h] = 1;
          }
        }

        _reachable[h][h] = 1;
        _reachable_as_loop[h][h] = 1;
      }

      compute_transitive_closure(_reachable);
      compute_transitive_closure(_reachable_as_loop);
    }

    /* Warshall's algorithm to fill the adjacency matrix */
    void compute_transitive_closure(VVB& mat) {
      int n = (int)mat.size();
      for (int k = 0; k < n; ++k) {
        for (int i = 0; i < n; ++i) {
          for (int j = 0; j < n; ++j) {
            if (mat[i][k] and mat[k][j]) { mat[i][j] = 1;}
          }
        }
      }
    }

    /* name interval states */
    void set_interval_state() {

      _state.clear();
      for (int hr=0; hr<M; ++hr)
        for (int hl=hr; 0<=hl; --hl)
          if (_reachable[hl][hr])
            _state.push_back(IS{int(_state.size()), hl, hr});

      _nodes_to_state.assign(M, VIS(M, IS{-1, -1, -1}));
      for (auto const& s: _state) {_nodes_to_state[s.l][s.r] = s;}

      _loop_state.clear();
      for (auto const& s: _state) {
        if(_reachable_as_loop[s.l][s.r]) _loop_state.push_back(s);
      }
    }

    /* name transitions between interval states */
    void set_interval_state_trans() {

      _loop_right_trans.assign(_state.size(), VIS());
      for (auto const& s: _state) {
        auto const cr = _node[s.r];
        if (cr=='z' or cr=='.' or
            cr=='*' or cr=='o') {
          for (int h: _edge_to[s.r]) {
            if (s.l <= h and _reachable[s.l][h]) {
              auto const& s1 = n2s(s.l, h);
              _loop_right_trans[s.id].push_back(s1);
            }
          }
        }
      }

      _loop_left_trans.assign(_state.size(), VIS());
      for (auto const& s: _state) {
        auto const cl = _node[s.l];
        if (cl=='z' or cl=='.' or
            cl=='*' or cl=='o') {
          for (int h: _edge_to[s.l]) {
            if (h <= s.r and _reachable[h][s.r]) {
              auto const& s1 = n2s(h, s.r);
              _loop_left_trans[s1.id].push_back(s);
            }
          }
        }
      }

      _pair_trans.assign(_state.size(), VIS());
      for (int hr=0; hr<M; ++ hr) {
        if (')'==_node[hr]) {
          int kl = _pair[hr];
          for (int hl: _edge_to[kl]) {
            auto const& s = n2s(hl, hr);
            for (int kr: _edge_to[hr]) {
              if(_reachable[kl][kr]) {
                _pair_trans[s.id].push_back(n2s(kl, kr));
              }
            }
          }
        }
      }

      for (auto const& s: _state) {
        if ('z'==_node[s.r] or
            'o'==_node[s.r] or
            '*'==_node[s.r]) {
          for  (int hl: _edge_from[s.l]) {
            if ('z'==_node[hl] or
                'o'==_node[hl] or
                '*'==_node[hl]) {
              for (int hr: _edge_to[s.r]) {
                if (_reachable[hl][hr]) {
                  _pair_trans[s.id].push_back(n2s(hl, hr));
                }
              }
            }
          }
        }
      }
    }

    void set_states_pairs() {
      _loop_loop_states.clear();
        for (auto const s2: loop_state()) {
          for (auto const s3: loop_state()) {
            if (s3.r < s2.l or
                !reachable(s2.r, s3.l) or
                !reachable(s2.l, s3.r)) continue;
            IS const s = n2s(s2.l, s3.r);
            IS const s1 = n2s(s2.r, s3.l);
            _loop_loop_states.push_back(VIS{s,s1,s2,s3});
          }
        }
    }

    void save() {
      for (auto const& si: _state) {
        say(si.id, ":", si.l, si.r);
        for (auto const& sij: _loop_right_trans[si.id]) {
          say("\tright:", sij.l, sij.r, sij.id);
        }
        for (auto const& sij: _loop_left_trans[si.id]) {
          say("\tleft:", sij.l, sij.r, sij.id);
        }
        for (auto const& sij: _pair_trans[si.id]) {
          say("\tpair:", sij.l, sij.r, sij.id);
        }
      }
    }
  };
}

#endif /* profile_hmm_h */
