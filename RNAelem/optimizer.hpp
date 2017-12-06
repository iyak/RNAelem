//
//  optimizer.hpp
//  RNAelem
//

#ifndef optimizer_h
#define optimizer_h

namespace iyak {

  class Naive{

    double _step = 1.e-1;
    double _converged = false;
    double _precision = 1.;

    public:

    void set_precision(double p) {_precision = p;}

    bool converged() {return _converged;}
    void reset() {_converged = false;}
    void update(double &a, double g) {
      a -= _step * g;
      if (_precision < fabs(g)) _converged = false;
    }
    template <typename T>
      void update(std::vector<T> &a, std::vector<T> &g) {
        _converged = true;
        for (int i = 0; i < a.size(); ++ i) {
          update(a[i], g[i]);
        }
      }
  };

  class ExMax {

    double _eps;
    double _fn;
    double _iter;

    int _maxit;
    V x0;
    V _x;
    V _ex;

    public:

    void set_eps(double e) {_eps = e;}
    void set_maxit(int m) {_maxit = m;}

    template <typename Block1, typename Blcok2>
      void maximize(const V& x0, Block1 e_step, Blcok2 m_step) {
        _iter = 0;
        _ex.assign(x0.size(), 0);
        _x = x0;
        //_fn = -inf;
        _fn=0.;
        double old_fn;
        do {
          old_fn = _fn;
          e_step(_x, _fn, _ex);
          m_step(_x, _ex);
          ++_iter;
        } while (_eps < _fn-old_fn and _iter < _maxit);
        say("em converged. iter:", _iter);
        say("final value:", _fn);
        return;
      }
  };

  class Adam {
    double _alpha=0.001; //stepsize
    double _beta1=0.9,_beta2=0.999; //exponential decay rates
    V _m{},_v{}; //moment vectors
    double _m0i=0,_v0i=0; //initial moment vector elements
    double _epsilon=1e-8; //stepsize auxiliary parameter
    V _rho{}; //regularization term
    V _x{}; //parameter vector
    V _xl{},_xu{}; //boudary of x vector
    VI _xb{}; //boundary type of x vector. 1:lower 2:upper 3:both
    VI _xr{}; //regularization of x. 0:None, 1:L1, 2:L2
    double _rgl_term; //sum of regularization terms
    int _t; //time
  public:
    Adam(){}
    void set_hp(double m0i,double v0i,double alpha,double beta1,
                double beta2,double epsilon){
      _alpha=alpha;
      _beta1=beta1;_beta2=beta2;
      _m0i=m0i;_v0i=v0i;
      _epsilon=epsilon;
    }
    void set_bounds(V& lower,V& upper,VI& xb){
      _xl=lower;_xu=upper;_xb=xb;
    }
    void set_regularization(VI const& xr,V const& rho){_xr=xr;_rho=rho;}
    bool converged(V const& gr,double y){
      return norm2(gr)<(y+1.)*1.e-5;
    }
    void before_update(double& y,V& gr){
      /* regularize */
      for(int i=0;i<size(_x);++i){
        if(1==_xr[i]){
          _rgl_term+=_rho[i]*abs(_x[i]);
          gr[i]+=_rho[i]*0<_x[i]?1:-1;
        }
        else if(2==_xr[i]){
          _rgl_term+=_rho[i]*_x[i]*_x[i]/2.;
          gr[i]+=_rho[i]*_x[i];
        }
      }
      y+=_rgl_term;
    }
    void after_update(double& y,V& gr){
      /* clip bounds */
      for(int i=0;i<size(_x);++i){
        if(1==_xb[i] or 3==_xb[i])
          if(_x[i]<_xl[i])
            _x[i]=_xl[i];
        if(2==_xb[i] or 3==_xb[i])
          if(_xu[i]<_x[i])
            _x[i]=_xu[i];
      }
    }
    template<class T>
    void minimize(
        T& f,//stochastic objective function with parameter x
        V const& x0, //initial parameter vector
        int max_iter=100 //max count of iteration
    ){
      _t=0;
      _m.assign(size(x0),_m0i);
      _v.assign(size(x0),_v0i);
      _x=x0;
      _xl.resize(size(x0),-inf);
      _xu.resize(size(x0),inf);
      _xr.resize(size(x0),0);
      V gr(size(x0),0);
      double y=0;
      double beta1t=_beta1,beta2t=_beta2;
      do{
        ++_t;
        f(_x,y,gr);
        before_update(y,gr);
        beta1t*=_beta1;
        beta2t*=_beta2;
        for(int i=0;i<size(_x);++i){
          _m[i]+=(1.-_beta1)*(gr[i]-_m[i]);
          _v[i]+=(1.-_beta2)*(gr[i]*gr[i]-_v[i]);
          double mhat=_m[i]/(1.-beta1t);
          double vhat=_v[i]/(1.-beta2t);
          _x[i]-=_alpha*mhat/(sqrt(vhat)+_epsilon);
        }
        after_update(y,gr);
        cry("iter:",_t,", y:",y,", |gr|:",norm2(gr),", p|x|:",_rgl_term);
      }while(not converged(gr,y) and _t<max_iter);
    }
    V& x(){return _x;}
    int itercount(){return _t-1;}
  };

  class Lbfgsb {
    public:
      using VC = vector<char>;
      enum Nbd {
        NBD_NONE = 0,
        NBD_LO   = 1,
        NBD_HI   = 3,
        NBD_BOTH = 2
      };
      int    c__1;
      int    c__11;
      int    _n;
      int    _m;
      V      _x;
      V      _l;
      V      _u;
      VI     _nbd;
      VI     _xr;
      V _rho;
      double _f;
      V      _g;
      double _factr;
      double _pgtol;
      int    _maxit;
      int    _iprint;
      VC     _task;
      VI     _lsave;
      VI     _isave;
      V      _dsave;
      V      _wa;
      VI     _iwa;
      int    _fncount;
      int    _grcount;
      int    _fdfcount;
      int    _fail;
      int    _iter;
      double _best_fn;
      V      _best_x;
      double _rgl_term;
      Lbfgsb():c__1(1),c__11(11),_n(0),_m(5),_f(HUGE_VAL),
      _factr(1.0e7),_pgtol(1.0e-3),_maxit(200),_iprint(0),_fncount(0),
      _grcount(0),_fdfcount(0),_fail(0),_best_fn(HUGE_VAL),_rgl_term(0){}
      double fn() const {return _best_fn;}
      double best_fn() const {return _best_fn;}
      const V& best_x() const { return _best_x;}
      const V& gr() const {return _g;}
      const V& x() const {return _x;}
      int iter() const { return _iter;}
      int fncount() const {return _fncount;}
      int grcount() const {return _grcount;}
      int fdfcount() const { return _fdfcount;}
      string message() const {return string(_task.begin(), _task.end());}
      void set_initial_point(const V& x) {_x = x; _n = (int)x.size();}
      void set_num_corrections(int m) {_m = m;}
      // type of bound. nbd[i] = {0=none, 1=l, 2=l&u, 3=u}
      void set_regularization(VI const& xr,V const& rho){_xr=xr;_rho=rho;}
      void set_bounds(const V& l, const V& u, 
          const VI& nbd) { _l = l; _u = u; _nbd = nbd;}
      void set_eps(double eps) { _pgtol = eps;}
      void set_factr(double factr) { _factr = factr;}
      void set_pgtol(double pgtol) { _pgtol = pgtol;}
      void set_maxit(int maxit) { _maxit = maxit;}
      //  iprint<0    no output is generated;
      //  iprint=0    print only one line at the last iteration;
      //  0<iprint<99 print also f and |proj g| every iprint iterations;
      //  iprint=99   print details of every iteration except n-vectors;
      //  iprint=100  print also the changes of active set and final x;
      //  iprint>100  print details of every iteration including x and g;
      //  When iprint > 0, the file iterate.dat will be created to summarize the iteration.
      void set_verbosity(int n) {_iprint = n;}
      void before_update(double& y,V& gr){
        /* regularize */
        _rgl_term=0.;
        for(int i=0;i<size(_x);++i){
          if(1==_xr[i]){
            _rgl_term+=_rho[i]*abs(_x[i]);
            gr[i]+=_rho[i]*0<_x[i]?1:-1;
          }
          else if(2==_xr[i]){
            _rgl_term+=_rho[i]*_x[i]*_x[i]/2.;
            gr[i]+=_rho[i]*_x[i];
          }
        }
        y+=_rgl_term;
      }
      template <typename Block>
        void minimize(const V& x0, Block& compute_fn_gr) {
          set_initial_point(x0);
          _task.resize(60);
          _lsave.resize(4);
          _isave.resize(44);
          _dsave.resize(29);
          _g.resize(_n);
          _wa.resize(2*_m*_n+4*_n+11*_m*_m+8*_m);
          _iwa.resize(3*_n);

          if ((int)_l.size() != _n || (int)_u.size() != _n || (int)_nbd.size() != _n) {
            _l.assign(_n, 0);
            _u.assign(_n, 0);
            _nbd.assign(_n, 0); // all variables are unbounded
          }

          _f = HUGE_VAL;
          _best_fn = _f;
          _best_x  = x0;
          fill(_task.begin(), _task.end(), 0);
          fill(_lsave.begin(), _lsave.end(), 0);
          fill(_dsave.begin(), _dsave.end(), 0);
          fill(_g.begin(), _g.end(), 0);
          /* this needs to be zeroed for snd in mainlb to be zeroed */
          fill(_wa.begin(), _wa.end(), 0);
          fill(_iwa.begin(), _iwa.end(), 0);

          _fail = 0;
          _iter = 0;
          _fdfcount = 0;
          strcpy(&_task[0], "START");
          while (true) {
            setulb(_n, _m, &_x[0], &_l[0], &_u[0], &_nbd[0], &_f, &_g[0], 
                _factr, &_pgtol, &_wa[0], &_iwa[0], &_task[0], 
                _iprint, &_lsave[0], &_isave[0], &_dsave[0]);
            if (strncmp(&_task[0], "FG", 2) == 0) {
              int st = compute_fn_gr(_x, _f, _g);
              before_update(_f,_g);
              ++_fdfcount;
              if (_f < _best_fn) { _best_fn = _f; _best_x = _x;}
              if (st) return;
              if (!std::isfinite(_f)) {
                print_message("L-BFGS-B needs finite values of 'fn'");
                break;
              }
            } else if (strncmp(&_task[0], "NEW_X", 5) == 0) {
              if (++_iter > _maxit) {
                _fail = 1;
                break;
              }
            } else if (strncmp(&_task[0], "WARN", 4) == 0) {
              _fail = 51;
              break;
            } else if (strncmp(&_task[0], "CONV", 4) == 0) {
              break;
            } else if (strncmp(&_task[0], "ERROR", 5) == 0) {
              _fail = 52;
              break;
            } else { /* some other condition that is not supposed to happen */
              _fail = 52;
              break;
            }
          }
          _fncount = _grcount = _isave[33];
          if (-2 < _iprint) {
            if (_iter < _maxit && _fail == 0) {
              print_message("lbfgsb converged. iter: %d\n", _iter);
            } else {
              print_message("stopped after %d iterations.\n", _iter);
            }
            print_message("final value: %f\n", _f);
          }
        }
      void print_message(const char* format, ...) {
        va_list ap;
        va_start(ap, format);
        std::vfprintf(stderr, format, ap);
        va_end(ap);
      }
      void timer(double* ttime) { *ttime = 0.0;}
      void setulb(int n, int m, double *x, double *l, double *u, int *nbd,
          double *f, double *g, double factr, double *pgtol,
          double *wa, int* iwa, char *task, int iprint,
          int *lsave, int *isave, double *dsave) {
        char csave[60] = {};
        int lsnd; 
        // int l1, l2, l3;
        int ld, lr, lt;
        int lz, lwa, lwn, lss, lws, lwt, lsy, lwy;
        /* Parameter adjustments */
        --wa;
        --isave;

        if (strncmp(task, "START", 5) == 0) {
          isave[1] = m * n;
          isave[2] = m * m;
          isave[3] = m * m << 2;
          isave[4] = 1;
          isave[5] = isave[4] + isave[1];
          isave[6] = isave[5] + isave[1];
          isave[7] = isave[6] + isave[2];
          isave[8] = isave[7] + isave[2];
          isave[9] = isave[8];
          isave[10] = isave[9] + isave[2];
          isave[11] = isave[10] + isave[3];
          isave[12] = isave[11] + isave[3];
          isave[13] = isave[12] + n;
          isave[14] = isave[13] + n;
          isave[15] = isave[14] + n;
          isave[16] = isave[15] + n;
        }
        // l1 = isave[1];
        // l2 = isave[2];
        // l3 = isave[3];
        lws = isave[4];
        lwy = isave[5];
        lsy = isave[6];
        lss = isave[7];
        lwt = isave[9];
        lwn = isave[10];
        lsnd = isave[11];
        lz = isave[12];
        lr = isave[13];
        ld = isave[14];
        lt = isave[15];
        lwa = isave[16];
        mainlb(n, m, x, l, u, nbd, f, g, factr, pgtol,
            &wa[lws], &wa[lwy], &wa[lsy],&wa[lss], &wa[lwt],&wa[lwn],
            &wa[lsnd], &wa[lz], &wa[lr], &wa[ld], &wa[lt], &wa[lwa],
            iwa, &iwa[n], &iwa[n << 1], task, iprint, csave, lsave, &isave[22], dsave);
      } 
      void mainlb(int n, int m, double *x, double *l, double *u, int *nbd, double *f, double *g,
          double factr, double *pgtol, double *ws, double * wy,
          double *sy, double *ss, double *wt, double *wn,
          double *snd, double *z, double *r, double *d,
          double *t, double *wa, int *indx, int *iwhere, int *indx2, char *task, int iprint,
          char *csave, int *lsave, int *isave, double *dsave) {
        int ws_offset=0, wy_offset=0, sy_offset=0, ss_offset=0, wt_offset=0;
        int wn_offset=0, snd_offset=0, i__1;
        double d__1, d__2;
        int head;
        double fold;
        int nact;
        double ddum;
        int info;
        double time;
        int nfgv;
        int ifun, iter, nint;
        char word[4]; /* allow for terminator */
        double time1, time2;
        int i, iback, k = 0; /* -Wall */
        double gdold;
        int nfree;
        int boxed;
        int itail;
        double theta;
        double dnorm;
        int nskip, iword;
        double xstep = 0.0, stpmx; /* xstep is printed before being used */
        double gd, dr, rr;
        int ileave;
        int itfile;
        double cachyt, epsmch;
        int updatd;
        double sbtime;
        int prjctd;
        int iupdat;
        int cnstnd;
        double sbgnrm;
        int nenter;
        double lnscht;
        int nintol;
        double dtd;
        int col;
        double tol;
        int wrk;
        double stp, cpu1, cpu2;
        /* Parameter adjustments */
        --indx2;
        --iwhere;
        --indx;
        --t;
        --d;
        --r;
        --z;
        --g;
        --nbd;
        --u;
        --l;
        --x;
        --wa;
        --lsave;
        --isave;
        --dsave;

        if (strncmp(task, "START", 5) == 0) {
          timer(&time1);
          /* Generate the current machine precision. */
          epsmch = DBL_EPSILON;

          fold = 0.;
          dnorm = 0.;
          cpu1 = 0.;
          gd = 0.;
          sbgnrm = 0.;
          stp = 0.;
          xstep = 0.;
          stpmx = 0.;
          gdold = 0.;
          dtd = 0.;
          /* Initialize counters and scalars when task='START'. */
          /* for the limited memory BFGS matrices: */
          col = 0;
          head = 1;
          theta = 1.;
          iupdat = 0;
          updatd = false;
          iback = 0;
          itail = 0;
          ifun = 0;
          iword = 0;
          nact = 0;
          ileave = 0;
          nenter = 0;
          /* for operation counts: */
          iter = 0;
          nfgv = 0;
          nint = 0;
          nintol = 0;
          nskip = 0;
          nfree = n;
          /* for stopping tolerance: */
          tol = factr * epsmch;
          /*	     for measuring running time: */
          cachyt = 0.;
          sbtime = 0.;
          lnscht = 0.;
          /* 'word' records the status of subspace solutions. */
          strcpy(word, "---");
          /* 'info' records the termination information. */
          info = 0;
          itfile = 0;
          /* Check the input arguments for errors. */
          errclb(n, m, factr, &l[1], &u[1], &nbd[1], task, &info, &k);
          if (strncmp(task, "ERROR", 5) == 0) {
            prn3lb(n, x+1, f, task, iprint, info, iter, nfgv, nintol, nskip, nact, sbgnrm,
                nint, word, iback, stp, xstep, k);
            return;
          }
          prn1lb(n, m, l+1, u+1, x+1, iprint, epsmch);
          /* Initialize iwhere & project x onto the feasible set. */
          active(n, &l[1], &u[1], &nbd[1], &x[1], &iwhere[1], iprint, &prjctd, &cnstnd, &boxed);
          /* The end of the initialization. */
        } else {
          /* restore local variables. */
          prjctd = lsave[1];
          cnstnd = lsave[2];
          boxed = lsave[3];
          updatd = lsave[4];

          nintol = isave[1];
          itfile = isave[3];
          iback = isave[4];
          nskip = isave[5];
          head = isave[6];
          col = isave[7];
          itail = isave[8];
          iter = isave[9];
          iupdat = isave[10];
          nint = isave[12];
          nfgv = isave[13];
          info = isave[14];
          ifun = isave[15];
          iword = isave[16];
          nfree = isave[17];
          nact = isave[18];
          ileave = isave[19];
          nenter = isave[20];

          theta = dsave[1];
          fold = dsave[2];
          tol = dsave[3];
          dnorm = dsave[4];
          epsmch = dsave[5];
          cpu1 = dsave[6];
          cachyt = dsave[7];
          sbtime = dsave[8];
          lnscht = dsave[9];
          time1 = dsave[10];
          gd = dsave[11];
          stpmx = dsave[12];
          sbgnrm = dsave[13];
          stp = dsave[14];
          gdold = dsave[15];
          dtd = dsave[16];
          /* After returning from the driver go to the point where execution */
          /* is to resume. */
          if (strncmp(task, "FG_LN", 5) == 0) goto L666;
          if (strncmp(task, "NEW_X", 5) == 0) goto L777;
          if (strncmp(task, "FG_ST", 5) == 0) goto L111;

          if (strncmp(task, "STOP", 4) == 0) {
            if (strncmp(task + 6, "CPU", 3) == 0) {
              /* restore the previous iterate. */
              dcopy(&n, &t[1], &c__1, &x[1], &c__1);
              dcopy(&n, &r[1], &c__1, &g[1], &c__1);
              *f = fold;
            }
            goto L999;
          }
        }
        /* Compute f0 and g0. */
        strcpy(task, "FG_START");
        /* return to the driver to calculate f and g; reenter at 111. */
        goto L1000;
L111:
        nfgv = 1;
        /* Compute the infinity norm of the (-) projected gradient. */
        projgr(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);

        if (iprint >= 1) print_message("iter: %d , f: %.5g , |gr|: %.5g , p|x|: %.5g\n",
            iter, *f, sbgnrm, _rgl_term);

        if (sbgnrm <= *pgtol) {
          /* terminate the algorithm. */
          strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
          goto L999;
        }
        /* -- the beginning of the loop -- */
L222:
        if (iprint >= 99) print_message("Iteration %5d\n", iter);
        iword = -1;

        if (! cnstnd && col > 0) {
          /* skip the search for GCP. */
          dcopy(&n, &x[1], &c__1, &z[1], &c__1);
          wrk = updatd;
          nint = 0;
          goto L333;
        }
        /* ccccccccccccccccccccccccccccccccccccccccccccccc */
        /*     Compute the Generalized Cauchy Point (GCP). */
        /* ccccccccccccccccccccccccccccccccccccccccccccccc */
        timer(&cpu1);
        cauchy(n, &x[1], &l[1], &u[1], &nbd[1], &g[1], &indx2[1], &iwhere[1], 
            &t[1], &d[1], &z[1], m, &wy[wy_offset], &ws[ws_offset], &sy[sy_offset], 
            &wt[wt_offset], &theta, &col, &head, &wa[1], &wa[(m<< 1) + 1], &wa[(m << 2) + 1], 
            &wa[m * 6 + 1], &nint, iprint, &sbgnrm, &info, &epsmch);
        if (info != 0) {
          /* singular triangular system detected; refresh the lbfgs memory. */
          if (iprint >= 1) print_message("%s\n%s\n", "Singular triangular system detected;",
              " refresh the lbfgs memory and restart the iteration.");
          info = 0;
          col = 0;
          head = 1;
          theta = 1.;
          iupdat = 0;
          updatd = false;
          timer(&cpu2);
          cachyt = cachyt + cpu2 - cpu1;
          goto L222;
        }
        timer(&cpu2);
        cachyt = cachyt + cpu2 - cpu1;
        nintol += nint;
        /* Count the entering and leaving variables for iter > 0; */
        /* find the index set of free and active variables at the GCP. */
        freev(n, &nfree, &indx[1], &nenter, &ileave, &indx2[1], &iwhere[1], &
            wrk, &updatd, &cnstnd, iprint, &iter);
        nact = n - nfree;
L333:
        /* If there are no free variables or B=theta*I, then */
        /* skip the subspace minimization. */
        if (nfree == 0 || col == 0) goto L555;

        /* ccccccccccccccccccccccccccccccccccccccccccccccccccccc */
        /*     Subspace minimization. */
        /* ccccccccccccccccccccccccccccccccccccccccccccccccccccc */
        timer(&cpu1);
        /* Form  the LEL^T factorization of the indefinite */
        /*	 matrix	   K = [-D -Y'ZZ'Y/theta     L_a'-R_z'	] */
        /*		       [L_a -R_z	   theta*S'AA'S ] */
        /*	 where	   E = [-I  0] */
        /*		       [ 0  I] */
        if (wrk) {
          formk(n, &nfree, &indx[1], &nenter, &ileave, &indx2[1], &iupdat,
              &updatd, &wn[wn_offset], &snd[snd_offset], m, &ws[ws_offset],
              &wy[wy_offset], &sy[sy_offset], &theta, &col, &head, &info);
        }
        if (info != 0) {
          /* nonpositive definiteness in Cholesky factorization; */
          /* refresh the lbfgs memory and restart the iteration. */
          if (iprint >= 0) print_message("%s\n%s\n",
              "Nonpositive definiteness in Cholesky factorization in formk;",
              " refresh the lbfgs memory and restart the iteration.");
          info = 0;
          col = 0;
          head = 1;
          theta = 1.;
          iupdat = 0;
          updatd = false;
          timer(&cpu2);
          sbtime = sbtime + cpu2 - cpu1;
          goto L222;
        }
        /* compute r=-Z'B(xcp-xk)-Z'g (using wa(2m+1)=W'(xcp-x) */
        /* from 'cauchy'). */
        cmprlb(n, m, &x[1], &g[1], &ws[ws_offset], &wy[wy_offset], &sy[sy_offset],
            &wt[wt_offset], &z[1], &r[1], &wa[1], &indx[1], &theta,
            &col, &head, &nfree, &cnstnd, &info);
        if (info != 0) {
          goto L444;
        }
        /* call the direct method. */
        subsm(n, m, &nfree, &indx[1], &l[1], &u[1], &nbd[1], &z[1], &r[1],
            &ws[ws_offset], &wy[wy_offset], &theta, &col, &head, &iword, &wa[1],
            &wn[wn_offset], iprint, &info);
L444:
        if (info != 0) {
          /* singular triangular system detected; */
          /* refresh the lbfgs memory and restart the iteration. */
          if (iprint >= 1) print_message("%s\n%s\n", "Singular triangular system detected;",
              " refresh the lbfgs memory and restart the iteration.");
          info = 0;
          col = 0;
          head = 1;
          theta = 1.;
          iupdat = 0;
          updatd = false;
          timer(&cpu2);
          sbtime = sbtime + cpu2 - cpu1;
          goto L222;
        }
        timer(&cpu2);
        sbtime = sbtime + cpu2 - cpu1;
L555:
        /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
        /*     Line search and optimality tests. */
        /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
        /* Generate the search direction d:=z-x. */
        i__1 = n;
        for (i = 1; i <= i__1; ++i) {
          d[i] = z[i] - x[i];
        }
        timer(&cpu1);
L666:
        lnsrlb(n, &l[1], &u[1], &nbd[1], &x[1], f, &fold, &gd, &gdold, &g[1],
            &d[1], &r[1], &t[1], &z[1], &stp, &dnorm, &dtd, &xstep,
            &stpmx, &iter, &ifun, &iback, &nfgv, &info, task, &boxed, &cnstnd,
            csave, &isave[22], &dsave[17]);
        if (info != 0 || iback >= 20) {
          /* restore the previous iterate. */
          dcopy(&n, &t[1], &c__1, &x[1], &c__1);
          dcopy(&n, &r[1], &c__1, &g[1], &c__1);
          *f = fold;
          if (col == 0) {
            /* abnormal termination. */
            if (info == 0) {
              info = -9;
              /* restore the actual number of f and g evaluations etc. */
              --nfgv;
              --ifun;
              --iback;
            }
            strcpy(task, "ERROR: ABNORMAL_TERMINATION_IN_LNSRCH");
            ++iter;
            goto L999;
          } else {
            /* refresh the lbfgs memory and restart the iteration. */
            if (iprint >= 1) print_message("%s\n%s\n", "Bad direction in the line search;",
                " refresh the lbfgs memory and restart the iteration.");
            if (info == 0) --nfgv;

            info = 0;
            col = 0;
            head = 1;
            theta = 1.;
            iupdat = 0;
            updatd = false;
            strcpy(task, "RESTART_FROM_LNSRCH");
            timer(&cpu2);
            lnscht = lnscht + cpu2 - cpu1;
            goto L222;
          }
        } else if (strncmp(task, "FG_LN", 5) == 0) {
          /* return to the driver for calculating f and g; reenter at 666. */
          goto L1000;
        } else {
          /* calculate and print out the quantities related to the new X. */
          timer(&cpu2);
          lnscht = lnscht + cpu2 - cpu1;
          ++iter;
          /* Compute the infinity norm of the projected (-)gradient. */
          projgr(n, &l[1], &u[1], &nbd[1], &x[1], &g[1], &sbgnrm);
          /* Print iteration information. */
          prn2lb(n, x+1, f, g+1, iprint, iter, nfgv, nact,
              sbgnrm, nint, word, iword, iback, stp, xstep);
          goto L1000;
        }
L777:
        /* Test for termination. */
        if (sbgnrm <= *pgtol) {
          /* terminate the algorithm. */
          strcpy(task, "CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL");
          goto L999;
        }
        /* Computing MAX */
        d__1 = std::abs(fold), d__2 = std::abs(*f), d__1 = std::max(d__1,d__2);
        ddum = std::max(d__1,1.);
        if (fold - *f <= tol * ddum) {
          /* terminate the algorithm. */
          strcpy(task, "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH");
          if (iback >= 10) info = -5;
          /* i.e., to issue a warning if iback>10 in the line search. */
          goto L999;
        }
        /* Compute d=newx-oldx, r=newg-oldg, rr=y'y and dr=y's. */
        i__1 = n;
        for (i = 1; i <= i__1; ++i) {
          r[i] = g[i] - r[i];
        }
        rr = ddot(&n, &r[1], &c__1, &r[1], &c__1);
        if (stp == 1.) {
          dr = gd - gdold;
          ddum = -gdold;
        } else {
          dr = (gd - gdold) * stp;
          dscal(&n, &stp, &d[1], &c__1);
          ddum = -gdold * stp;
        }
        if (dr <= epsmch * ddum) {
          /* skip the L-BFGS update. */
          ++nskip;
          updatd = false;
          if (iprint >= 1)
            print_message("ys=%10.3e  -gs=%10.3e, BFGS update SKIPPED\n", dr, ddum);
          goto L888;
        }
        /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
        /*     Update the L-BFGS matrix. */
        /* ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc */
        updatd = true;
        ++iupdat;
        /* Update matrices WS and WY and form the middle matrix in B. */
        matupd(n, m, &ws[ws_offset], &wy[wy_offset], &sy[sy_offset], 
            &ss[ss_offset], &d[1], &r[1], &itail, &iupdat, &col, &head,
            &theta, &rr, &dr, &stp, &dtd);
        /* Form the upper half of the pds T = theta*SS + L*D^(-1)*L'; */
        /* Store T in the upper triangular of the array wt; */
        /* Cholesky factorize T to J*J' with */
        /* J' stored in the upper triangular of wt. */
        formt(m, &wt[wt_offset], &sy[sy_offset], &ss[ss_offset], &col, &theta, &info);
        if (info != 0) {
          /* nonpositive definiteness in Cholesky factorization; */
          /* refresh the lbfgs memory and restart the iteration. */
          if (iprint >= 0)
            print_message("%s\n%s\n",
                "Nonpositive definiteness in Cholesky factorization in formk;",
                "   refresh the lbfgs memory and restart the iteration.");
          info = 0;
          col = 0;
          head = 1;
          theta = 1.;
          iupdat = 0;
          updatd = false;
          goto L222;
        }
        /* Now the inverse of the middle matrix in B is */
        /* [  D^(1/2)	 O ] [ -D^(1/2)	 D^(-1/2)*L' ] */
        /* [ -L*D^(-1/2)	 J ] [	0	 J'	     ] */
L888:
        /* --the end of the loop -- */
        goto L222;
L999:
        timer(&time2);
        time = time2 - time1;
        if (0<time)
          ; /* to make gcc happy */
L1000:
        /* Save local variables. */
        lsave[1] = prjctd;
        lsave[2] = cnstnd;
        lsave[3] = boxed;
        lsave[4] = updatd;
        isave[1] = nintol;
        isave[3] = itfile;
        isave[4] = iback;
        isave[5] = nskip;
        isave[6] = head;
        isave[7] = col;
        isave[8] = itail;
        isave[9] = iter;
        isave[10] = iupdat;
        isave[12] = nint;
        isave[13] = nfgv;
        isave[14] = info;
        isave[15] = ifun;
        isave[16] = iword;
        isave[17] = nfree;
        isave[18] = nact;
        isave[19] = ileave;
        isave[20] = nenter;
        dsave[1] = theta;
        dsave[2] = fold;
        dsave[3] = tol;
        dsave[4] = dnorm;
        dsave[5] = epsmch;
        dsave[6] = cpu1;
        dsave[7] = cachyt;
        dsave[8] = sbtime;
        dsave[9] = lnscht;
        dsave[10] = time1;
        dsave[11] = gd;
        dsave[12] = stpmx;
        dsave[13] = sbgnrm;
        dsave[14] = stp;
        dsave[15] = gdold;
        dsave[16] = dtd;
        prn3lb(n, x+1, f, task, iprint, info, iter, nfgv, nintol, nskip, nact, sbgnrm,
            nint, word, iback, stp, xstep, k);
      } 
      void active(int n, double *l, double *u, int *nbd, double *x, int *iwhere, int iprint,
          int *prjctd, int *cnstnd, int *boxed) {
        int nbdd, i;
        /* Parameter adjustments */
        --iwhere;
        --x;
        --nbd;
        --u;
        --l;

        /* Initialize nbdd, prjctd, cnstnd and boxed. */
        nbdd = 0;
        *prjctd = false;
        *cnstnd = false;
        *boxed = true;
        /* Project the initial x to the easible set if necessary. */
        for (i = 1; i <= n; ++i) {
          if (nbd[i] > 0) {
            if (nbd[i] <= 2 && x[i] <= l[i]) {
              if (x[i] < l[i]) {
                *prjctd = true;
                x[i] = l[i];
              }
              ++nbdd;
            } else if (nbd[i] >= 2 && x[i] >= u[i]) {
              if (x[i] > u[i]) {
                *prjctd = true;
                x[i] = u[i];
              }
              ++nbdd;
            }
          }
        }
        /* Initialize iwhere and assign values to cnstnd and boxed. */
        for (i = 1; i <= n; ++i) {
          if (nbd[i] != 2) {
            *boxed = false;
          }
          if (nbd[i] == 0) {
            /* this variable is always free */
            iwhere[i] = -1;
            /* otherwise set x(i)=mid(x(i), u(i), l(i)). */
          } else {
            *cnstnd = true;
            if (nbd[i] == 2 && u[i] - l[i] <= 0.) {
              /* this variable is always fixed */
              iwhere[i] = 3;
            } else {
              iwhere[i] = 0;
            }
          }
        }
        if (iprint >= 0) {
          if (*prjctd)
            print_message("The initial X is infeasible.  Restart with its projection.\n");
          if (!*cnstnd) print_message("This problem is unconstrained.\n");
        }
        if (iprint > 0) print_message("At X0, %d variables are exactly at the bounds\n", nbdd);
      }
      void bmv(int m, double *sy, double *wt, int *col, double *v, double *p, int *info) {
        int sy_dim1, sy_offset, wt_dim1, wt_offset, Col;
        int i, k;
        int i2;
        double sum;
        /* Parameter adjustments */
        wt_dim1 = m;
        wt_offset = 1 + wt_dim1 * 1;
        wt -= wt_offset;
        sy_dim1 = m;
        sy_offset = 1 + sy_dim1 * 1;
        sy -= sy_offset;
        --p;
        --v;

        if (*col == 0) return;

        /*	PART I: solve [	 D^(1/2)      O ] [ p1 ] = [ v1 ]
         *		      [ -L*D^(-1/2)   J ] [ p2 ]   [ v2 ].
         *	solve Jp2=v2+LD^(-1)v1.
         */
        Col = *col;
        p[*col + 1] = v[*col + 1];
        for (i = 2; i <= Col; ++i) {
          i2 = *col + i;
          sum = 0.;
          for (k = 1; k <= i - 1; ++k) {
            sum += sy[i + k * sy_dim1] * v[k] / sy[k + k * sy_dim1];
          }
          p[i2] = v[i2] + sum;
        }
        /* Solve the triangular system */
        dtrsl(&wt[wt_offset], &m, col, &p[*col + 1], &c__11, info);
        if (*info != 0) return;

        /* solve D^(1/2)p1=v1. */
        for (i = 1; i <= Col; ++i) {
          p[i] = v[i] / sqrt(sy[i + i * sy_dim1]);
        }

        /*	PART II: solve [ -D^(1/2)   D^(-1/2)*L'	 ] [ p1 ] = [ p1 ]
         *		       [  0	    J'		 ] [ p2 ]   [ p2 ].
         *	solve J^Tp2=p2.
         */
        dtrsl(&wt[wt_offset], &m, col, &p[*col + 1], &c__1, info);
        if (*info != 0) return;

        /* compute p1=-D^(-1/2)(p1-D^(-1/2)L'p2) */
        /*		 =-D^(-1/2)p1 + D^(-1)L'p2. */
        for (i = 1; i <= Col; ++i) p[i] = -p[i] / sqrt(sy[i + i * sy_dim1]);

        for (i = 1; i <= Col; ++i) {
          sum = 0.;
          for (k = i + 1; k <= Col; ++k) {
            sum += sy[k + i * sy_dim1] * p[*col + k] / sy[i + i * sy_dim1];
          }
          p[i] += sum;
        }
      }
      void cauchy(int n, double *x, double *l, double *u, int *nbd,
          double *g, int *iorder, int * iwhere, double *t,
          double *d, double *xcp, int m, double *wy, double *ws, double *sy, double *wt,
          double *theta, int *col, int *head, double *p,
          double *c, double *wbp, double *v, int *nint,
          int iprint, double *sbgnrm, int *info, double * epsmch) {
        int wy_dim1, wy_offset, ws_dim1, ws_offset, sy_dim1, sy_offset, wt_dim1, wt_offset, i__2;
        double d__1;
        double bkmin, dibp, dibp2, zibp, neggi, tsum;
        double f1, f2, f2_org__, dt, tj, tj0, tl= 0.0, tu=0.0, dtm, wmc, wmp, wmw;
        int i, j, ibp, iter, bnded, nfree, nleft, nbreak, ibkmin, pointr;
        int xlower, xupper, col2;
        /* Parameter adjustments */
        --xcp;
        --d;
        --t;
        --iwhere;
        --iorder;
        --g;
        --nbd;
        --u;
        --l;
        --x;
        --v;
        --wbp;
        --c;
        --p;
        wt_dim1 = m;    wt_offset = 1 + wt_dim1 * 1;    wt -= wt_offset;
        sy_dim1 = m;    sy_offset = 1 + sy_dim1 * 1;    sy -= sy_offset;
        ws_dim1 = n;    ws_offset = 1 + ws_dim1 * 1;    ws -= ws_offset;
        wy_dim1 = n;    wy_offset = 1 + wy_dim1 * 1;    wy -= wy_offset;

        /* Check the status of the variables, reset iwhere(i) if necessary;
         * compute the Cauchy direction d and the breakpoints t; initialize
         * the derivative f1 and the vector p = W'd (for theta = 1).
         */
        if (*sbgnrm <= 0.) {
          if (iprint >= 0) print_message("Subgnorm = 0.  GCP = X.\n");
          dcopy(&n, &x[1], &c__1, &xcp[1], &c__1);
          return;
        }
        bnded = true;
        nfree = n + 1;
        nbreak = 0;
        ibkmin = 0;
        bkmin = 0.;
        col2 = *col << 1;
        f1 = 0.;
        if (iprint >= 99) print_message("\n---CAUCHY entered---\n\n");

        /* We set p to zero and build it up as we determine d. */
        for (i = 1; i <= col2; ++i)
          p[i] = 0.;

        /* In the following loop we determine for each variable its bound */
        /* status and its breakpoint, and update p accordingly. */
        /* Smallest breakpoint is identified. */
        for (i = 1; i <= n; ++i) {
          neggi = -g[i];
          if (iwhere[i] != 3 && iwhere[i] != -1) {
            /* if x(i) is not a constant and has bounds, */
            /* compute the difference between x(i) and its bounds. */
            if (nbd[i] <= 2) {
              tl = x[i] - l[i];
            }
            if (nbd[i] >= 2) {
              tu = u[i] - x[i];
            }
            /* If a variable is close enough to a bound */
            /* we treat it as at bound. */
            xlower = nbd[i] <= 2 && tl <= 0.;
            xupper = nbd[i] >= 2 && tu <= 0.;
            /*		reset iwhere(i). */
            iwhere[i] = 0;
            if (xlower) {
              if (neggi <= 0.) {
                iwhere[i] = 1;
              }
            } else if (xupper) {
              if (neggi >= 0.) {
                iwhere[i] = 2;
              }
            } else {
              if (std::abs(neggi) <= 0.) {
                iwhere[i] = -3;
              }
            }
          }
          pointr = *head;
          if (iwhere[i] != 0 && iwhere[i] != -1) {
            d[i] = 0.;
          } else {
            d[i] = neggi;
            f1 -= neggi * neggi;
            /* calculate p := p - W'e_i* (g_i). */
            i__2 = *col;
            for (j = 1; j <= i__2; ++j) {
              p[j] += wy[i + pointr * wy_dim1] * neggi;
              p[*col + j] += ws[i + pointr * ws_dim1] * neggi;
              pointr = pointr % m + 1;
            }
            if (nbd[i] <= 2 && nbd[i] != 0 && neggi < 0.) {
              /* x(i) + d(i) is bounded; compute t(i). */
              ++nbreak;
              iorder[nbreak] = i;
              t[nbreak] = tl / (-neggi);
              if (nbreak == 1 || t[nbreak] < bkmin) {
                bkmin = t[nbreak];
                ibkmin = nbreak;
              }
            } else if (nbd[i] >= 2 && neggi > 0.) {
              /* x(i) + d(i) is bounded; compute t(i). */
              ++nbreak;
              iorder[nbreak] = i;
              t[nbreak] = tu / neggi;
              if (nbreak == 1 || t[nbreak] < bkmin) {
                bkmin = t[nbreak];
                ibkmin = nbreak;
              }
            } else {/* x(i) + d(i) is not bounded. */
              --nfree;
              iorder[nfree] = i;
              if (std::abs(neggi) > 0.)
                bnded = false;
            }
          }
        }
        /* The indices of the nonzero components of d are now stored */
        /* in iorder(1),...,iorder(nbreak) and iorder(nfree),...,iorder(n). */
        /* The smallest of the nbreak breakpoints is in t(ibkmin)=bkmin. */
        if (*theta != 1.) {
          /* complete the initialization of p for theta not= one. */
          dscal(col, theta, &p[*col + 1], &c__1);
        }
        /* Initialize GCP xcp = x. */
        dcopy(&n, &x[1], &c__1, &xcp[1], &c__1);
        if (nbreak == 0 && nfree == n + 1) {
          /* is a zero vector, return with the initial xcp as GCP. */
          if (iprint > 100) {
            print_message("Cauchy X =  ");
            for(i = 1; i <= n; i++) print_message("%g ", xcp[i]);
            print_message("\n");
          }
          return;
        }
        /* Initialize c = W'(xcp - x) = 0. */
        for (j = 1; j <= col2; ++j)
          c[j] = 0.;

        /* Initialize derivative f2. */
        f2 = -(*theta) * f1;
        f2_org__ = f2;
        if (*col > 0) {
          bmv(m, &sy[sy_offset], &wt[wt_offset], col, &p[1], &v[1], info);
          if (*info != 0) return;

          f2 -= ddot(&col2, &v[1], &c__1, &p[1], &c__1);
        }
        dtm = -f1 / f2;
        tsum = 0.;
        *nint = 1;
        if (iprint >= 99) print_message("There are %d  breakpoints\n", nbreak);

        /* If there are no breakpoints, locate the GCP and return. */
        if (nbreak == 0) goto L888;

        nleft = nbreak;
        iter = 1;
        tj = 0.;
        /* ------------------- the beginning of the loop ------------------------- */
L777:
        /* Find the next smallest breakpoint; */
        /* compute dt = t(nleft) - t(nleft + 1). */
        tj0 = tj;
        if (iter == 1) {
          /* Since we already have the smallest breakpoint we need not do */
          /* heapsort yet. Often only one breakpoint is used and the */
          /* cost of heapsort is avoided. */
          tj = bkmin;
          ibp = iorder[ibkmin];
        } else {
          if (iter == 2) {
            /* Replace the already used smallest breakpoint with the */
            /* breakpoint numbered nbreak > nlast, before heapsort call. */
            if (ibkmin != nbreak) {
              t[ibkmin] = t[nbreak];
              iorder[ibkmin] = iorder[nbreak];
            }
          }
          /* Update heap structure of breakpoints */
          /* (if iter=2, initialize heap). */
          hpsolb(nleft, &t[1], &iorder[1], iter - 2);
          tj = t[nleft];
          ibp = iorder[nleft];
        }
        dt = tj - tj0;

        if (dt != 0 && iprint >=  100) {
          print_message("\nPiece    %3i f1, f2 at start point %11.4e %11.4e\n", *nint, f1, f2);
          print_message("Distance to the next break point =  %11.4e\n", dt);
          print_message("Distance to the stationary point =  %11.4e\n", dtm);
        }

        /* If a minimizer is within this interval, */
        /* locate the GCP and return. */
        if (dtm < dt) goto L888;

        /* Otherwise fix one variable and */
        /* reset the corresponding component of d to zero. */
        tsum += dt;
        --nleft;
        ++iter;
        dibp = d[ibp];
        d[ibp] = 0.;
        if (dibp > 0.) {
          zibp = u[ibp] - x[ibp];
          xcp[ibp] = u[ibp];
          iwhere[ibp] = 2;
        } else {
          zibp = l[ibp] - x[ibp];
          xcp[ibp] = l[ibp];
          iwhere[ibp] = 1;
        }
        if (iprint >= 100) print_message("Variable  %d  is fixed.\n", ibp);
        if (nleft == 0 && nbreak == n) {
          /* all n variables are fixed, */
          /* return with xcp as GCP. */
          dtm = dt;
          goto L999;
        }
        /* Update the derivative information. */
        ++(*nint);
        dibp2 = dibp * dibp;
        /* Update f1 and f2. */
        /* temporarily set f1 and f2 for col=0. */
        f1 += dt * f2 + dibp2 - *theta * dibp * zibp;
        f2 -= *theta * dibp2;
        if (*col > 0) {
          /*			    update c = c + dt*p. */
          daxpy(&col2, &dt, &p[1], &c__1, &c[1], &c__1);
          /* choose wbp, */
          /* the row of W corresponding to the breakpoint encountered. */
          pointr = *head;
          for (j = 1; j <= *col; ++j) {
            wbp[j] = wy[ibp + pointr * wy_dim1];
            wbp[*col + j] = *theta * ws[ibp + pointr * ws_dim1];
            pointr = pointr % m + 1;
          }
          /* compute (wbp)Mc, (wbp)Mp, and (wbp)M(wbp)'. */
          bmv(m, &sy[sy_offset], &wt[wt_offset], col, &wbp[1], &v[1], info);
          if (*info != 0) return;

          wmc = ddot(&col2,  &c[1], &c__1, &v[1], &c__1);
          wmp = ddot(&col2,  &p[1], &c__1, &v[1], &c__1);
          wmw = ddot(&col2,&wbp[1], &c__1, &v[1], &c__1);
          /* update p = p - dibp*wbp. */
          d__1 = -dibp;
          daxpy(&col2, &d__1, &wbp[1], &c__1, &p[1], &c__1);
          /* complete updating f1 and f2 while col > 0. */
          f1 += dibp * wmc;
          f2 += (2. * dibp * wmp - dibp2 * wmw);
        }
        if(f2 < (d__1 = *epsmch * f2_org__)) f2 = d__1;
        if (nleft > 0) {
          dtm = -f1 / f2;
          goto L777;
          /* to repeat the loop for unsearched intervals. */
        } else if (bnded) {
          f1 = 0.;
          f2 = 0.;
          dtm = 0.;
        } else {
          dtm = -f1 / f2;
        }
        /* ---the end of the loop --- */
L888:
        if (iprint >= 99) {
          print_message("\nGCP found in this segment\n");
          print_message("Piece    %3i f1, f2 at start point %11.4e %11.4e\n", *nint,f1,f2);
          print_message("Distance to the stationary point =  %11.4e\n", dtm);
        }

        if (dtm <= 0.) {
          dtm = 0.;
        }
        tsum += dtm;
        /* Move free variables (i.e., the ones w/o breakpoints) and */
        /* the variables whose breakpoints haven't been reached. */
        daxpy(&n, &tsum, &d[1], &c__1, &xcp[1], &c__1);
L999:
        /* Update c = c + dtm*p = W'(x^c - x) */
        /* which will be used in computing r = Z'(B(x^c - x) + g). */
        if (*col > 0) {
          daxpy(&col2, &dtm, &p[1], &c__1, &c[1], &c__1);
        }
        if (iprint >= 100) {
          print_message("Cauchy X =  ");
          for(i = 1; i <= n; i++) print_message("%g ", xcp[i]);
          print_message("\n");
        }

        if (iprint >= 99) print_message("\n--- exit CAUCHY---\n\n");
      }
      void cmprlb(int n, int m, double *x, double *g, double *ws, double *wy, double *sy,
          double *wt, double *z, double *r, double *wa,
          int *indx, double *theta, int *col, int *head, int *nfree, int *cnstnd, int *info) {
        int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset;
        int wt_dim1, wt_offset, Col, n_f;
        int i, j, k;
        double a1, a2;
        int pointr;
        /* Parameter adjustments */
        --indx;
        --r;
        --z;
        --g;
        --x;
        --wa;
        wt_dim1 = m;
        wt_offset = 1 + wt_dim1 * 1;
        wt -= wt_offset;
        sy_dim1 = m;
        sy_offset = 1 + sy_dim1 * 1;
        sy -= sy_offset;
        wy_dim1 = n;
        wy_offset = 1 + wy_dim1 * 1;
        wy -= wy_offset;
        ws_dim1 = n;
        ws_offset = 1 + ws_dim1 * 1;
        ws -= ws_offset;

        Col = *col;
        if (! (*cnstnd) && Col > 0) {
          for (i = 1; i <= n; ++i)
            r[i] = -g[i];
        }
        else {
          n_f = *nfree;
          for (i = 1; i <= n_f; ++i) {
            k = indx[i];
            r[i] = -(*theta) * (z[k] - x[k]) - g[k];
          }
          bmv(m, &sy[sy_offset], &wt[wt_offset], col,
              &wa[(m << 1) + 1], &wa[1], info);
          if (*info != 0) {
            *info = -8;
            return;
          }
          pointr = *head;
          for (j = 1; j <= Col; ++j) {
            a1 = wa[j];
            a2 = *theta * wa[Col + j];
            for (i = 1; i <= n_f; ++i) {
              k = indx[i];
              r[i] += wy[k + pointr * wy_dim1] * a1 +
                ws[k + pointr * ws_dim1] * a2;
            }
            pointr = pointr % m + 1;
          }
        }
      } 
      void errclb(int n, int m, double factr, double *l, double *u,
          int *nbd, char *task, int *info, int *k) {
        int i;
        /* Parameter adjustments */
        --nbd;
        --u;
        --l;

        /* Check the input arguments for errors. */
        if (n <= 0)
          strcpy(task, "ERROR: N .LE. 0");
        if (m <= 0)
          strcpy(task, "ERROR: M .LE. 0");
        if (factr < 0.)
          strcpy(task, "ERROR: FACTR .LT. 0");

        /* Check the validity of the arrays nbd(i), u(i), and l(i). */
        for (i = 1; i <= n; ++i) {
          if (nbd[i] < 0 || nbd[i] > 3) {
            strcpy(task, "ERROR: INVALID NBD");
            *info = -6;
            *k = i;
          }
          if (nbd[i] == 2) {
            if (l[i] > u[i]) {
              strcpy(task, "ERROR: NO FEASIBLE SOLUTION");
              *info = -7;
              *k = i;
            }
          }
        }
      }
      void formk(int n, int *nsub, int *ind, int * nenter, int *ileave,
          int *indx2, int *iupdat, int * updatd, double *wn,
          double *wn1, int m, double *ws, double *wy, double *sy,
          double *theta, int *col, int *head, int *info) {
        int wn_dim1, wn_offset, wn1_dim1, wn1_offset, ws_dim1, ws_offset;
        int wy_dim1, wy_offset, sy_dim1, sy_offset, i__1, i__2;
        int dend, pend;
        int upcl;
        double temp1, temp2, temp3, temp4;
        int i, k;
        int ipntr, jpntr, k1, m2, dbegin, is, js, iy, jy, pbegin, is1, js1, col2;
        /* Parameter adjustments */
        --indx2;
        --ind;
        sy_dim1 = m;
        sy_offset = 1 + sy_dim1 * 1;
        sy -= sy_offset;
        wy_dim1 = n;
        wy_offset = 1 + wy_dim1 * 1;
        wy -= wy_offset;
        ws_dim1 = n;
        ws_offset = 1 + ws_dim1 * 1;
        ws -= ws_offset;
        wn1_dim1 = 2 * m;
        wn1_offset = 1 + wn1_dim1 * 1;
        wn1 -= wn1_offset;
        wn_dim1 = 2 * m;
        wn_offset = 1 + wn_dim1 * 1;
        wn -= wn_offset;

        /* Form the lower triangular part of */
        /*	 WN1 = [Y' ZZ'Y	  L_a'+R_z'] */
        /*	       [L_a+R_z	  S'AA'S   ] */
        /* where L_a is the strictly lower triangular part of S'AA'Y */
        /* R_z is the upper triangular part of S'ZZ'Y. */
        if (*updatd) {
          if (*iupdat > m) {/* shift old part of WN1. */
            i__1 = m - 1;
            for (jy = 1; jy <= i__1; ++jy) {
              js = m + jy;
              i__2 = m - jy;
              dcopy(&i__2, &wn1[jy + 1 + (jy + 1)* wn1_dim1], &c__1,
                  &wn1[jy + jy * wn1_dim1], &c__1);
              dcopy(&i__2, &wn1[js + 1 + (js + 1)* wn1_dim1], &c__1,
                  &wn1[js + js * wn1_dim1], &c__1);
              i__2 = m - 1;
              dcopy(&i__2, &wn1[m + 2 + (jy + 1) * wn1_dim1], &c__1,
                  &wn1[m + 1 + jy * wn1_dim1], &c__1);
            }
          }
          /* put new rows in blocks (1,1), (2,1) and (2,2). */
          pbegin = 1;
          pend = *nsub;
          dbegin = *nsub + 1;
          dend = n;
          iy = *col;
          is = m + *col;
          ipntr = *head + *col - 1;
          if (ipntr > m) {
            ipntr -= m;
          }
          jpntr = *head;
          i__1 = *col;
          for (jy = 1; jy <= i__1; ++jy) {
            js = m + jy;
            temp1 = 0.;
            temp2 = 0.;
            temp3 = 0.;
            /* compute element jy of row 'col' of Y'ZZ'Y */
            for (k = pbegin; k <= pend; ++k) {
              k1 = ind[k];
              temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
            }
            /* compute elements jy of row 'col' of L_a and S'AA'S */
            for (k = dbegin; k <= dend; ++k) {
              k1 = ind[k];
              temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
              temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
            }
            wn1[iy + jy * wn1_dim1] = temp1;
            wn1[is + js * wn1_dim1] = temp2;
            wn1[is + jy * wn1_dim1] = temp3;
            jpntr = jpntr % m + 1;
          }
          /* put new column in block (2,1). */
          jy = *col;
          jpntr = *head + *col - 1;
          if (jpntr > m) {
            jpntr -= m;
          }
          ipntr = *head;
          i__1 = *col;
          for (i = 1; i <= i__1; ++i) {
            is = m + i;
            temp3 = 0.;
            /* compute element i of column 'col' of R_z */
            for (k = pbegin; k <= pend; ++k) {
              k1 = ind[k];
              temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
            }
            ipntr = ipntr % m + 1;
            wn1[is + jy * wn1_dim1] = temp3;
          }
          upcl = *col - 1;
        } else {
          upcl = *col;
        }
        /* modify the old parts in blocks (1,1) and (2,2) due to changes */
        /* in the set of free variables. */
        ipntr = *head;
        for (iy = 1; iy <= upcl; ++iy) {
          is = m + iy;
          jpntr = *head;
          for (jy = 1; jy <= iy; ++jy) {
            js = m + jy;
            temp1 = 0.;
            temp2 = 0.;
            temp3 = 0.;
            temp4 = 0.;
            for (k = 1; k <= *nenter; ++k) {
              k1 = indx2[k];
              temp1 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
              temp2 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
            }
            for (k = *ileave; k <= n; ++k) {
              k1 = indx2[k];
              temp3 += wy[k1 + ipntr * wy_dim1] * wy[k1 + jpntr * wy_dim1];
              temp4 += ws[k1 + ipntr * ws_dim1] * ws[k1 + jpntr * ws_dim1];
            }
            wn1[iy + jy * wn1_dim1] = wn1[iy + jy * wn1_dim1] + temp1 - temp3;
            wn1[is + js * wn1_dim1] = wn1[is + js * wn1_dim1] - temp2 + temp4;
            jpntr = jpntr % m + 1;
          }
          ipntr = ipntr % m + 1;
        }
        /* modify the old parts in block (2,1). */
        ipntr = *head;
        for (is = m + 1; is <= m + upcl; ++is) {
          jpntr = *head;
          for (jy = 1; jy <= upcl; ++jy) {
            temp1 = 0.;
            temp3 = 0.;
            for (k = 1; k <= *nenter; ++k) {
              k1 = indx2[k];
              temp1 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
            }
            for (k = *ileave; k <= n; ++k) {
              k1 = indx2[k];
              temp3 += ws[k1 + ipntr * ws_dim1] * wy[k1 + jpntr * wy_dim1];
            }
            if (is <= jy + m) {
              wn1[is + jy * wn1_dim1] +=  temp1 - temp3;
            } else {
              wn1[is + jy * wn1_dim1] += -temp1 + temp3;
            }
            jpntr = jpntr % m + 1;
          }
          ipntr = ipntr % m + 1;
        }
        /* Form the upper triangle of WN = [D+Y' ZZ'Y/theta	  -L_a'+R_z' ] */
        /*				       [-L_a +R_z	 S'AA'S*theta] */
        m2 = m << 1;
        i__1 = *col;
        for (iy = 1; iy <= i__1; ++iy) {
          is = *col + iy;
          is1 = m + iy;
          i__2 = iy;
          for (jy = 1; jy <= i__2; ++jy) {
            js = *col + jy;
            js1 = m + jy;
            wn[jy + iy * wn_dim1] = wn1[iy + jy * wn1_dim1] / *theta;
            wn[js + is * wn_dim1] = wn1[is1 + js1 * wn1_dim1] * *theta;
          }
          i__2 = iy - 1;
          for (jy = 1; jy <= i__2; ++jy) {
            wn[jy + is * wn_dim1] = -wn1[is1 + jy * wn1_dim1];
          }
          i__2 = *col;
          for (jy = iy; jy <= i__2; ++jy) {
            wn[jy + is * wn_dim1] = wn1[is1 + jy * wn1_dim1];
          }
          wn[iy + iy * wn_dim1] += sy[iy + iy * sy_dim1];
        }
        /*  Form the upper triangle of */
        /* WN= [  LL'		  L^-1(-L_a'+R_z')] */
        /*	[(-L_a +R_z)L'^-1   S'AA'S*theta  ] */
        /*  first Cholesky factor (1,1) block of wn to get LL' */
        /*		    with L' stored in the upper triangle of wn. */
        dpofa(&wn[wn_offset], &m2, col, info);
        if (*info != 0) {
          *info = -1;
          return;
        }
        /* then form L^-1(-L_a'+R_z') in the (1,2) block. */
        col2 = *col << 1;
        for (js = *col + 1; js <= col2; ++js) {
          dtrsl(&wn[wn_offset], &m2, col,
              &wn[js * wn_dim1 + 1], &c__11, info);
        }
        /* Form S'AA'S*theta + (L^-1(-L_a'+R_z'))'L^-1(-L_a'+R_z') in the */
        /* upper triangle of (2,2) block of wn. */
        for (is = *col + 1; is <= col2; ++is) {
          for (js = is; js <= col2; ++js) {
            wn[is + js * wn_dim1] +=
              ddot(col, &wn[is * wn_dim1 + 1], &c__1,
                  &wn[js * wn_dim1 + 1], &c__1);
          }
        }
        /* Cholesky factorization of (2,2) block of wn. */
        dpofa(&wn[*col + 1 + (*col + 1) * wn_dim1], &m2, col, info);
        if (*info != 0) *info = -2;
      }
      void formt(int m, double *wt, double *sy, double *ss,
          int *col, double *theta, int *info) {
        int wt_dim1, wt_offset, sy_dim1, sy_offset, ss_dim1, ss_offset, i__1;
        double ddum;
        int i, j, k;
        int k1;
        /* Parameter adjustments */
        ss_dim1 = m;
        ss_offset = 1 + ss_dim1 * 1;
        ss -= ss_offset;
        sy_dim1 = m;
        sy_offset = 1 + sy_dim1 * 1;
        sy -= sy_offset;
        wt_dim1 = m;
        wt_offset = 1 + wt_dim1 * 1;
        wt -= wt_offset;

        /* Form the upper half of  T = theta*SS + L*D^(-1)*L', */
        /* store T in the upper triangle of the array wt. */
        i__1 = *col;
        for (j = 1; j <= i__1; ++j) {
          wt[j * wt_dim1 + 1] = *theta * ss[j * ss_dim1 + 1];
        }
        for (i = 2; i <= i__1; ++i) {
          for (j = i; j <= i__1; ++j) {
            k1 = std::min(i,j) - 1;
            ddum = 0.;
            for (k = 1; k <= k1; ++k) {
              ddum += sy[i + k * sy_dim1] * sy[j + k * sy_dim1] / sy[k +
                k * sy_dim1];
            }
            wt[i + j * wt_dim1] = ddum + *theta * ss[i + j * ss_dim1];
          }
        }
        /* Cholesky factorize T to J*J' with */
        /* J' stored in the upper triangle of wt. */
        dpofa(&wt[wt_offset], &m, col, info);
        if (*info != 0) *info = -3;
      } 
      void freev(int n, int *nfree, int *indx,
          int *nenter, int *ileave, int *indx2, int *iwhere,
          int *wrk, int *updatd, int *cnstnd, int iprint, int *iter) {
        int i__1;
        int iact, i, k;
        /* Parameter adjustments */
        --iwhere;
        --indx2;
        --indx;

        *nenter = 0;
        *ileave = n + 1;
        if (*iter > 0 && *cnstnd) {/* count the entering and leaving variables. */
          i__1 = *nfree;
          for (i = 1; i <= i__1; ++i) {
            k = indx[i];
            if (iwhere[k] > 0) {
              --(*ileave);
              indx2[*ileave] = k;
              if (iprint >= 100)
                print_message("Variable %d leaves the set of free variables\n",
                    k);
            }
          }
          for (i = *nfree + 1; i <= n; ++i) {
            k = indx[i];
            if (iwhere[k] <= 0) {
              ++(*nenter);
              indx2[*nenter] = k;
              if (iprint >= 100)
                print_message("Variable %d enters the set of free variables\n",
                    k);
            }
            if (iprint >= 100)
              print_message("%d variables leave; %d variables enter\n",
                  n + 1 - *ileave, *nenter);
          }
        }
        *wrk = *ileave < n + 1 || *nenter > 0 || *updatd;
        /* Find the index set of free and active variables at the GCP. */
        *nfree = 0;
        iact = n + 1;
        for (i = 1; i <= n; ++i) {
          if (iwhere[i] <= 0) {
            ++(*nfree);
            indx[*nfree] = i;
          } else {
            --iact;
            indx[iact] = i;
          }
        }
        if (iprint >= 99)
          print_message("%d  variables are free at GCP on iteration %d\n",
              *nfree, *iter + 1);
      }
      void hpsolb(int n, double *t, int *iorder, int iheap) {
        double ddum;
        int i, j, k, indxin, indxou;
        double out;
        /* Parameter adjustments */
        --iorder;
        --t;

        if (iheap == 0) {
          /* Rearrange the elements t(1) to t(n) to form a heap. */
          for (k = 2; k <= n; ++k) {
            ddum = t[k];
            indxin = iorder[k];
            /* Add ddum to the heap. */
            i = k;
h_loop:
            if (i > 1) {
              j = i / 2;
              if (ddum < t[j]) {
                t[i] = t[j];
                iorder[i] = iorder[j];
                i = j;
                goto h_loop;
              }
            }
            t[i] = ddum;
            iorder[i] = indxin;
          }
        }
        /* Assign to 'out' the value of t(1), the least member of the heap, */
        /* and rearrange the remaining members to form a heap as */
        /* elements 1 to n-1 of t. */
        if (n > 1) {
          i = 1;
          out = t[1];
          indxou = iorder[1];
          ddum = t[n];
          indxin = iorder[n];
          /* Restore the heap */
Loop:
          j = i + i;
          if (j <= n - 1) {
            if (t[j + 1] < t[j]) {
              ++j;
            }
            if (t[j] < ddum) {
              t[i] = t[j];
              iorder[i] = iorder[j];
              i = j;
              goto Loop;
            }
          }
          t[i] = ddum;
          iorder[i] = indxin;
          /* Put the least member in t(n). */
          t[n] = out;
          iorder[n] = indxou;
        }
      }
      void lnsrlb(int n, double *l, double *u, int *nbd, double *x, double *f, double *fold,
          double *gd, double *gdold, double *g, double *d,
          double *r, double *t, double *z, double *stp,
          double *dnorm, double *dtd, double *xstep,
          double *stpmx, int *iter, int *ifun, int *iback, int *nfgv,
          int *info, char *task, int *boxed, int *cnstnd,
          char *csave, int *isave, double *dsave) {
        /* For dcsrch(): */
        const double stpmin = 0.;
        const double ftol = .001;
        const double gtol = .9;
        const double xtol = .1;

        double d1;
        int i;
        double a1, a2;
        /* Parameter adjustments */
        --z;
        --t;
        --r;
        --d;
        --g;
        --x;
        --nbd;
        --u;
        --l;

        if (strncmp(task, "FG_LN", 5) == 0) goto L556;

        *dtd = ddot(&n, &d[1], &c__1, &d[1], &c__1);
        *dnorm = sqrt(*dtd);
        /* Determine the maximum step length. */
        *stpmx = 1e10;
        if (*cnstnd) {
          if (*iter == 0) {
            *stpmx = 1.;
          } else {
            for (i = 1; i <= n; ++i) {
              a1 = d[i];
              if (nbd[i] != 0) {
                if (a1 < 0. && nbd[i] <= 2) {
                  a2 = l[i] - x[i];
                  if (a2 >= 0.) {
                    *stpmx = 0.;
                  } else if (a1 * *stpmx < a2) {
                    *stpmx = a2 / a1;
                  }
                } else if (a1 > 0. && nbd[i] >= 2) {
                  a2 = u[i] - x[i];
                  if (a2 <= 0.) {
                    *stpmx = 0.;
                  } else if (a1 * *stpmx > a2) {
                    *stpmx = a2 / a1;
                  }
                }
              }
            }
          }
        }
        if (*iter == 0 && ! (*boxed)) {
          d1 = 1. / *dnorm;
          *stp = std::min(d1,*stpmx);
        } else {
          *stp = 1.;
        }
        dcopy(&n, &x[1], &c__1, &t[1], &c__1);
        dcopy(&n, &g[1], &c__1, &r[1], &c__1);
        *fold = *f;
        *ifun = 0;
        *iback = 0;
        strcpy(csave, "START");
L556:
        *gd = ddot(&n, &g[1], &c__1, &d[1], &c__1);
        if (*ifun == 0) {
          *gdold = *gd;
          if (*gd >= 0.) {
            /* the directional derivative >=0. */
            /* Line search is impossible. */
            *info = -4;
            return;
          }
        }
        dcsrch(f, gd, stp,
            ftol, gtol, xtol,
            stpmin, *stpmx,
            csave, isave, dsave);
        *xstep = *stp * *dnorm;
        if (strncmp(csave, "CONV", 4) != 0 && strncmp(csave, "WARN", 4) != 0) {
          strcpy(task, "FG_LNSRCH");
          ++(*ifun);
          ++(*nfgv);
          *iback = *ifun - 1;
          if (*stp == 1.) {
            dcopy(&n, &z[1], &c__1, &x[1], &c__1);
          } else {
            for (i = 1; i <= n; ++i) x[i] = *stp * d[i] + t[i];
          }
        } else {
          strcpy(task, "NEW_X");
        }
      }
      void matupd(int n, int m, double *ws, double *wy, double *sy, double *ss, double *d,
          double *r, int *itail, int *iupdat, int *col,
          int *head, double *theta, double *rr, double *dr, double *stp, double *dtd) {
        int ws_dim1, ws_offset, wy_dim1, wy_offset, sy_dim1, sy_offset, ss_dim1, ss_offset, i__1, i__2;
        int j;
        int pointr;
        /* Parameter adjustments */
        --r;
        --d;
        ss_dim1 = m;
        ss_offset = 1 + ss_dim1 * 1;
        ss -= ss_offset;
        sy_dim1 = m;
        sy_offset = 1 + sy_dim1 * 1;
        sy -= sy_offset;
        wy_dim1 = n;
        wy_offset = 1 + wy_dim1 * 1;
        wy -= wy_offset;
        ws_dim1 = n;
        ws_offset = 1 + ws_dim1 * 1;
        ws -= ws_offset;

        /* Set pointers for matrices WS and WY. */
        if (*iupdat <= m) {
          *col = *iupdat;
          *itail = (*head + *iupdat - 2) % m + 1;
        } else {
          *itail = *itail % m + 1;
          *head = *head % m + 1;
        }
        /* Update matrices WS and WY. */
        dcopy(&n, &d[1], &c__1, &ws[*itail * ws_dim1 + 1], &c__1);
        dcopy(&n, &r[1], &c__1, &wy[*itail * wy_dim1 + 1], &c__1);
        /* Set theta=yy/ys. */
        *theta = *rr / *dr;
        /* Form the middle matrix in B. */
        /* update the upper triangle of SS, */
        /* and the lower triangle of SY: */
        if (*iupdat > m) {
          /* move old information */
          i__1 = *col - 1;
          for (j = 1; j <= i__1; ++j) {
            dcopy(&j, &ss[(j + 1) * ss_dim1 + 2], &c__1, &ss[j * ss_dim1 + 1], &c__1);
            i__2 = *col - j;
            dcopy(&i__2, &sy[j + 1 + (j + 1) * sy_dim1], &c__1, &sy[j + j * sy_dim1], &c__1);
          }
        }
        /* add new information: the last row of SY */
        /* and the last column of SS: */
        pointr = *head;
        i__1 = *col - 1;
        for (j = 1; j <= i__1; ++j) {
          sy[*col + j * sy_dim1] =
            ddot(&n, &d[1], &c__1, &wy[pointr * wy_dim1 + 1], &c__1);
          ss[j + *col * ss_dim1] =
            ddot(&n, &ws[pointr * ws_dim1 + 1], &c__1, &d[1], &c__1);
          pointr = pointr % m + 1;
        }
        if (*stp == 1.) {
          ss[*col + *col * ss_dim1] = *dtd;
        } else {
          ss[*col + *col * ss_dim1] = *stp * *stp * *dtd;
        }
        sy[*col + *col * sy_dim1] = *dr;
      }
      void projgr(int n, double *l, double *u,
          int *nbd, double *x, double *g, double *sbgnrm) {
        int i;
        double gi, d__1;

        *sbgnrm = 0.;
        for (i = 0; i < n; ++i) {
          gi = g[i];
          if (nbd[i] != 0) {
            if (gi < 0.) {
              if (nbd[i] >= 2) {
                if(gi < (d__1 = x[i] - u[i])) gi = d__1;
              }
            } else {
              if (nbd[i] <= 2) {
                if(gi > (d__1 = x[i] - l[i]))
                  gi = d__1;
              }
            }
          }
          if(*sbgnrm < (d__1 = std::abs(gi))) *sbgnrm = d__1;
        }
      }
      void subsm(int n, int m, int *nsub, int *ind,
          double *l, double *u, int *nbd, double *x,
          double *d, double *ws, double *wy, double *theta,
          int *col, int *head, int *iword, double *wv,
          double *wn, int /*iprint*/, int *info) {
        int ws_offset, wn_dim1, wn_offset;
        double alpha, dk, temp1, temp2;
        int i, j, k, m2, js, jy, pointr, ibd = 0, col2, ns;
        /* Parameter adjustments */
        --d;
        --u;
        --l;
        --x;
        --ind;
        --nbd;
        --wv;
        wn_dim1 = 2 * m;
        wn_offset = 1 + wn_dim1 * 1;
        wn -= wn_offset;
        /* ws[] and wy[] are both  [n x m ] :*/
        ws_offset = 1 + n * 1;
        ws -= ws_offset;
        wy -= ws_offset;

        ns = *nsub;
        if (ns <= 0) return;

        /* Compute wv = W'Zd. */
        pointr = *head;
        for (i = 1; i <= *col; ++i) {
          temp1 = 0.;
          temp2 = 0.;
          for (j = 1; j <= ns; ++j) {
            k = ind[j];
            temp1 += wy[k + pointr * n] * d[j];
            temp2 += ws[k + pointr * n] * d[j];
          }
          wv[i] = temp1;
          wv[*col + i] = *theta * temp2;
          pointr = pointr % m + 1;
        }
        /* Compute wv:=K^(-1)wv. */
        m2 = m << 1;
        col2 = *col << 1;
        dtrsl(&wn[wn_offset], &m2, &col2, &wv[1], &c__11, info);
        if (*info != 0) return;

        for (i = 1; i <= *col; ++i)
          wv[i] = -wv[i];

        dtrsl(&wn[wn_offset], &m2, &col2, &wv[1], &c__1, info);
        if (*info != 0) return;

        /* Compute d = (1/theta)d + (1/theta**2)Z'W wv. */
        pointr = *head;
        for (jy = 1; jy <= *col; ++jy) {
          js = *col + jy;
          for (i = 1; i <= ns; ++i) {
            k = ind[i];
            d[i] += (wy[k + pointr * n] * wv[jy] / *theta +
                ws[k + pointr * n] * wv[js]);
          }
          pointr = pointr % m + 1;
        }

        for (i = 1; i <= ns; ++i) d[i] /= *theta;

        /* Backtrack to the feasible region. */
        alpha = 1.;
        temp1 = alpha;
        for (i = 1; i <= ns; ++i) {
          k = ind[i];
          dk = d[i];
          if (nbd[k] != 0) {
            if (dk < 0. && nbd[k] <= 2) {
              temp2 = l[k] - x[k];
              if (temp2 >= 0.) {
                temp1 = 0.;
              } else if (dk * alpha < temp2) {
                temp1 = temp2 / dk;
              }
            } else if (dk > 0. && nbd[k] >= 2) {
              temp2 = u[k] - x[k];
              if (temp2 <= 0.) {
                temp1 = 0.;
              } else if (dk * alpha > temp2) {
                temp1 = temp2 / dk;
              }
            }
            if (temp1 < alpha) {
              alpha = temp1;
              ibd = i;
            }
          }
        }
        if (alpha < 1.) {
          dk = d[ibd];
          k = ind[ibd];
          if (dk > 0.) {
            x[k] = u[k];
            d[ibd] = 0.;
          } else if (dk < 0.) {
            x[k] = l[k];
            d[ibd] = 0.;
          }
        }
        for (i = 1; i <= ns; ++i) x[ind[i]] += alpha * d[i];

        *iword = (alpha < 1.) ? 1 : 0;
      }
      void dcsrch(double *f, double *g, double *stp, double ftol, double gtol, double xtol,
          double stpmin, double stpmax, char *task, int *isave, double *dsave) {
        int stage;
        double finit, ginit, width, ftest, gtest, stmin, stmax, width1, fm;
        double gm, fx, fy, gx, gy;
        int brackt;
        double fxm, fym, gxm, gym, stx, sty;
        /* Parameter adjustments */
        --dsave;
        --isave;

        /* Initialization block. */
        if (strncmp(task, "START", 5) == 0) {
          /* Check the input arguments for errors. */
          if (*stp < stpmin)	strcpy(task, "ERROR: STP .LT. STPMIN");
          if (*stp > stpmax)	strcpy(task, "ERROR: STP .GT. STPMAX");
          if (*g >= 0.)		strcpy(task, "ERROR: INITIAL G .GE. ZERO");
          if (ftol < 0.)		strcpy(task, "ERROR: FTOL .LT. ZERO");
          if (gtol < 0.)		strcpy(task, "ERROR: GTOL .LT. ZERO");
          if (xtol < 0.)		strcpy(task, "ERROR: XTOL .LT. ZERO");
          if (stpmin < 0.)	        strcpy(task, "ERROR: STPMIN .LT. ZERO");
          if (stpmax < stpmin)	strcpy(task, "ERROR: STPMAX .LT. STPMIN");

          /* Exit if there are errors on input. */
          if (strncmp(task, "ERROR", 5) == 0) return;

          /* Initialize local variables. */
          brackt = false;
          stage = 1;
          finit = *f;
          ginit = *g;
          gtest = ftol * ginit;
          width = stpmax - stpmin;
          width1 = width / .5;
          /* The variables stx, fx, gx contain the values of the step, */
          /* function, and derivative at the best step. */
          /* The variables sty, fy, gy contain the value of the step, */
          /* function, and derivative at sty. */
          /* The variables stp, f, g contain the values of the step, */
          /* function, and derivative at stp. */
          stx = 0.;	fx = finit;	gx = ginit;
          sty = 0.;	fy = finit;	gy = ginit;
          stmin = 0.;
          stmax = *stp + *stp * 4.;
          strcpy(task, "FG");
          goto L1000;
        } else {
          /* Restore local variables. */
          if (isave[1] == 1) {
            brackt = true;
          } else {
            brackt = false;
          }
          stage = isave[2];
          ginit = dsave[1];
          gtest = dsave[2];
          gx = dsave[3];
          gy = dsave[4];
          finit = dsave[5];
          fx = dsave[6];
          fy = dsave[7];
          stx = dsave[8];
          sty = dsave[9];
          stmin = dsave[10];
          stmax = dsave[11];
          width = dsave[12];
          width1 = dsave[13];
        }
        /* If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
        /* algorithm enters the second stage. */
        ftest = finit + *stp * gtest;
        if (stage == 1 && *f <= ftest && *g >= 0.) {
          stage = 2;
        }
        /*	Test for warnings. */
        if (brackt && (*stp <= stmin || *stp >= stmax))
          strcpy(task, "WARNING: ROUNDING ERRORS PREVENT PROGRESS");
        if (brackt && stmax - stmin <= xtol * stmax) strcpy(task, "WARNING: XTOL TEST SATISFIED");
        if (*stp == stpmax && *f <= ftest && *g <= gtest) strcpy(task, "WARNING: STP = STPMAX");
        if (*stp == stpmin && (*f > ftest || *g >= gtest)) strcpy(task, "WARNING: STP = STPMIN");
        /*	Test for convergence. */
        if (*f <= ftest && std::abs(*g) <= gtol * (-ginit)) strcpy(task, "CONVERGENCE");
        /*	Test for termination. */
        if (strncmp(task, "WARN", 4) == 0 || strncmp(task, "CONV", 4) == 0) goto L1000;

        /* A modified function is used to predict the step during the */
        /* first stage if a lower function value has been obtained but */
        /* The decrease is not sufficient. */
        if (stage == 1 && *f <= fx && *f > ftest) {
          /* Define the modified function and derivative values. */
          fm = *f - *stp * gtest;
          fxm = fx - stx * gtest;
          fym = fy - sty * gtest;
          gm = *g - gtest;
          gxm = gx - gtest;
          gym = gy - gtest;
          /* Call dcstep to update stx, sty, and to compute the new step. */
          dcstep(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, &fm, &gm, &brackt, &
              stmin, &stmax);
          /* Reset the function and derivative values for f. */
          fx = fxm + stx * gtest;
          fy = fym + sty * gtest;
          gx = gxm + gtest;
          gy = gym + gtest;
        } else {
          /* Call dcstep to update stx, sty, and to compute the new step. */
          dcstep(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, &stmin, &
              stmax);
        }
        /* Decide if a bisection step is needed. */
        if (brackt) {
          if (std::abs(sty - stx) >= width1 * .66) {
            *stp = stx + (sty - stx) * .5;
          }
          width1 = width;
          width = std::abs(sty - stx);
        }
        /* Set the minimum and maximum steps allowed for stp. */
        if (brackt) {
          stmin = std::min(stx,sty);
          stmax = std::max(stx,sty);
        } else {
          stmin = *stp + (*stp - stx) * 1.1;
          stmax = *stp + (*stp - stx) * 4.;
        }
        /* Force the step to be within the bounds stpmax and stpmin. */
        if(*stp < stpmin) *stp = stpmin;
        if(*stp > stpmax) *stp = stpmax;

        /* If further progress is not possible, let stp be the best */
        /* point obtained during the search. */
        if ((brackt && (*stp <= stmin || *stp >= stmax)) ||
            (brackt && (stmax - stmin <= xtol * stmax))) {
          *stp = stx;
        }
        /* Obtain another function and derivative. */
        strcpy(task, "FG");
L1000:
        /* Save local variables. */
        if (brackt) {
          isave[1] = 1;
        } else {
          isave[1] = 0;
        }
        isave[2] = stage;
        dsave[1] = ginit;
        dsave[2] = gtest;
        dsave[3] = gx;
        dsave[4] = gy;
        dsave[5] = finit;
        dsave[6] = fx;
        dsave[7] = fy;
        dsave[8] = stx;
        dsave[9] = sty;
        dsave[10] = stmin;
        dsave[11] = stmax;
        dsave[12] = width;
        dsave[13] = width1;
      }
      void dcstep(double *stx, double *fx, double *dx,
          double *sty, double *fy, double *dy, double *stp,
          double *fp, double *dp, int *brackt, double *stpmin, double *stpmax) {
        double d__1, d__2;
        double sgnd, stpc, stpf, stpq, p, q, gamm, r__, s, theta;

        sgnd = *dp * (*dx / std::abs(*dx));
        /* First case: A higher function value. The minimum is bracketed. */
        /* If the cubic step is closer to stx than the quadratic step, the */
        /* cubic step is taken, otherwise the average of the cubic and */
        /* quadratic steps is taken. */
        if (*fp > *fx) {
          theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
          /* Computing MAX */
          d__1 = std::abs(theta), d__2 = std::abs(*dx),
               d__1 = std::max(d__1,d__2), d__2 = std::abs(*dp);
          s = std::max(d__1,d__2);
          /* Computing 2nd power */
          d__1 = theta / s;
          gamm = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
          if (*stp < *stx) {
            gamm = -gamm;
          }
          p = gamm - *dx + theta;
          q = gamm - *dx + gamm + *dp;
          r__ = p / q;
          stpc = *stx + r__ * (*stp - *stx);
          stpq = *stx + *dx / ((*fx - *fp) / (*stp - *stx) + *dx) / 2. * (*stp
              - *stx);
          if (std::abs(stpc - *stx) < std::abs(stpq - *stx)) {
            stpf = stpc;
          } else {
            stpf = stpc + (stpq - stpc) / 2.;
          }
          *brackt = true;
          /* Second case: A lower function value and derivatives of opposite */
          /* sign. The minimum is bracketed. If the cubic step is farther from */
          /* stp than the secant step, the cubic step is taken, otherwise the */
          /* secant step is taken. */
        } else if (sgnd < 0.) {
          theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
          /* Computing MAX */
          d__1 = std::abs(theta), d__2 = std::abs(*dx),
               d__1 = std::max(d__1,d__2), d__2 = std::abs(*dp);
          s = std::max(d__1,d__2);
          /* Computing 2nd power */
          d__1 = theta / s;
          gamm = s * sqrt(d__1 * d__1 - *dx / s * (*dp / s));
          if (*stp > *stx) {
            gamm = -gamm;
          }
          p = gamm - *dp + theta;
          q = gamm - *dp + gamm + *dx;
          r__ = p / q;
          stpc = *stp + r__ * (*stx - *stp);
          stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
          if (std::abs(stpc - *stp) > std::abs(stpq - *stp)) {
            stpf = stpc;
          } else {
            stpf = stpq;
          }
          *brackt = true;
          /* Third case: A lower function value, derivatives of the same sign, */
          /* and the magnitude of the derivative decreases. */
        } else if (std::abs(*dp) < std::abs(*dx)) {
          /* The cubic step is computed only if the cubic tends to infinity */
          /* in the direction of the step or if the minimum of the cubic */
          /* is beyond stp. Otherwise the cubic step is defined to be the */
          /* secant step. */
          theta = (*fx - *fp) * 3. / (*stp - *stx) + *dx + *dp;
          /* Computing MAX */
          d__1 = std::abs(theta), d__2 = std::abs(*dx),
               d__1 = std::max(d__1,d__2), d__2 = std::abs(*dp);
          s = std::max(d__1,d__2);
          /* The case gamm = 0 only arises if the cubic does not tend */
          /* to infinity in the direction of the step. */
          /* Computing MAX */
          /* Computing 2nd power */
          d__1 = theta / s;
          d__1 = d__1 * d__1 - *dx / s * (*dp / s);
          gamm = d__1 < 0 ? 0. : s * sqrt(d__1);
          if (*stp > *stx) {
            gamm = -gamm;
          }
          p = gamm - *dp + theta;
          q = gamm + (*dx - *dp) + gamm;
          r__ = p / q;
          if (r__ < 0. && gamm != 0.) {
            stpc = *stp + r__ * (*stx - *stp);
          } else if (*stp > *stx) {
            stpc = *stpmax;
          } else {
            stpc = *stpmin;
          }
          stpq = *stp + *dp / (*dp - *dx) * (*stx - *stp);
          if (*brackt) {
            /* A minimizer has been bracketed. If the cubic step is */
            /* closer to stp than the secant step, the cubic step is */
            /* taken, otherwise the secant step is taken. */
            if (std::abs(stpc - *stp) < std::abs(stpq - *stp)) {
              stpf = stpc;
            } else {
              stpf = stpq;
            }
            d__1 = *stp + (*sty - *stp) * .66;
            if (*stp > *stx) {
              stpf = std::min(d__1,stpf);
            } else {
              stpf = std::max(d__1,stpf);
            }
          } else {
            /* A minimizer has not been bracketed. If the cubic step is */
            /* farther from stp than the secant step, the cubic step is */
            /* taken, otherwise the secant step is taken. */
            if (std::abs(stpc - *stp) > std::abs(stpq - *stp)) {
              stpf = stpc;
            } else {
              stpf = stpq;
            }
            stpf = std::min(*stpmax,stpf);
            stpf = std::max(*stpmin,stpf);
          }
          /* Fourth case: A lower function value, derivatives of the */
          /* same sign, and the magnitude of the derivative does not */
          /* decrease. If the minimum is not bracketed, the step is either */
          /* stpmin or stpmax, otherwise the cubic step is taken. */
        } else {
          if (*brackt) {
            theta = (*fp - *fy) * 3. / (*sty - *stp) + *dy + *dp;
            /* Computing MAX */
            d__1 = std::abs(theta), d__2 = std::abs(*dy), d__1 = std::max(d__1,d__2), d__2 =
              std::abs(*dp);
            s = std::max(d__1,d__2);
            /* Computing 2nd power */
            d__1 = theta / s;
            gamm = s * sqrt(d__1 * d__1 - *dy / s * (*dp / s));
            if (*stp > *sty) {
              gamm = -gamm;
            }
            p = gamm - *dp + theta;
            q = gamm - *dp + gamm + *dy;
            r__ = p / q;
            stpc = *stp + r__ * (*sty - *stp);
            stpf = stpc;
          } else if (*stp > *stx) {
            stpf = *stpmax;
          } else {
            stpf = *stpmin;
          }
        }
        /* Update the interval which contains a minimizer. */
        if (*fp > *fx) {
          *sty = *stp;
          *fy = *fp;
          *dy = *dp;
        } else {
          if (sgnd < 0.) {
            *sty = *stx;
            *fy = *fx;
            *dy = *dx;
          }
          *stx = *stp;
          *fx = *fp;
          *dx = *dp;
        }
        /* Compute the new step. */
        *stp = stpf;
      } 
      void pvector(const char *title, double *x, int n) {
        int i;
        print_message("%s ", title);
        for (i = 0; i < n; i++) print_message("%g ", x[i]);
        print_message("\n");
      }
      void prn1lb(int n, int m, double *l, double *u, double *x, int iprint, double epsmch) {
        if (iprint >=  0) {
          print_message("N = %d, M = %d machine precision = %g\n", n, m, epsmch);
          if (iprint >= 100){
            pvector("L =", l, n);
            pvector("X0 =",x, n);
            pvector("U =", u, n);
          }
        }
      }
      void prn2lb(int n, double *x, double *f, double *g, int iprint,
          int iter, int /*nfgv*/, int /*nact*/, double sbgnrm,
          int /*nint*/, char* /*word*/, int /*iword*/, int iback,
          double /*stp*/, double xstep) {
        if (iprint >=  99) {
          print_message("LINE SEARCH %d times; norm of step = %g\n", iback, xstep);
          if (iprint > 100) {
            pvector("X =", x, n);
            pvector("G =", g, n);
          }
        } else if (iprint > 0 && iter%iprint == 0) {
          print_message("iter: %d , f: %.5g , |gr|: %.5g , p|x|: %.5g\n",
              iter, *f, sbgnrm,_rgl_term);
        }
      }
      void prn3lb(int n, double *x, double *f, char *task, int iprint,
          int info, int iter, int nfgv, int nintol, int nskip,
          int nact, double sbgnrm, int /*nint*/,
          char* /*word*/, int /*iback*/, double /*stp*/, double /*xstep*/,
          int k) {
        if(strncmp(task, "CONV", 4) == 0) {
          if (iprint >= 0) print_message("iterations: %d\n"
              "function evaluations: %d\n"
              "segments explored during Cauchy searches: %d\n"
              "BFGS updates skipped: %d\n"
              "active bounds at final generalized Cauchy point: %d\n"
              "norm of the final projected gradient: %g\n"
              "final function value: %g\n",
              iter, nfgv, nintol, nskip, nact, sbgnrm, *f);

          if (iprint >= 100) pvector("X =", x, n);
          if (iprint >= 1) print_message("F = %g\n", *f);
        }
        if (iprint >= 0) {
          switch(info) {
            case -1: print_message("Matrix in 1st Cholesky factorization in formk is not Pos. Def.\n"); break;
            case -2: print_message("Matrix in 2st Cholesky factorization in formk is not Pos. Def.\n"); break;
            case -3: print_message("Matrix in the Cholesky factorization in formt is not Pos. Def.\n"); break;
            case -4: print_message("Derivative >= 0, backtracking line search impossible\n"); break;
            case -5: print_message("l(%d) > u(%d). No feasible solution\n", k, k); break;
            case -6: print_message("Input nbd(%d) is invalid\n", k); break;
            case -7: print_message("Warning:  more than 10 function and gradient evaluations in the last line search\n"); break;
            case -8: print_message("The triangular system is singular\n"); break;
            case -9: print_message("Line search cannot locate an adequate point after 20 function and gradient evaluations\n"); break;
            default: break;
          }
        }
      }
      int dpofa(double *a, int *lda, int *n, int *info) {
        const double eps = 1e-14;
        int a_dim1, a_offset, i__1, i__2, i__3;
        double d__1;
        int j, k;
        double s, t;
        int jm1;
        /* Parameter adjustments */
        a_dim1 = *lda;
        a_offset = 1 + a_dim1;
        a -= a_offset;

        /* begin block with ...exits to 40 */
        i__1 = *n;
        for (j = 1; j <= i__1; ++j) {
          *info = j;
          s = 0.;
          jm1 = j - 1;
          if (jm1 < 1) goto L20;

          i__2 = jm1;
          for (k = 1; k <= i__2; ++k) {
            i__3 = k - 1;
            t = a[k + j * a_dim1] - ddot(&i__3, &a[k * a_dim1 + 1], &c__1, &
                a[j * a_dim1 + 1], &c__1);
            t /= a[k + k * a_dim1];
            a[k + j * a_dim1] = t;
            s += t * t;
          }
L20:
          s = a[j + j * a_dim1] - s;
          /* ......exit */
          /*  if (s .le. 0.0d0) go to 40 */
          if (s <= eps * (d__1 = a[j + j * a_dim1], std::abs(d__1))) goto L40;

          a[j + j * a_dim1] = sqrt(s);
        }
        *info = 0;
L40:
        return 0;
      }
      int dtrsl(double *t, int *ldt, int *n, double *b, int *job, int *info) {
        int t_dim1, t_offset, i__1, i__2;
        int j, jj, case__;
        double temp;
        /* begin block permitting ...exits to 150 */
        /* check for zero diagonal elements. */
        /* Parameter adjustments */
        t_dim1 = *ldt;
        t_offset = 1 + t_dim1;
        t -= t_offset;
        --b;

        i__1 = *n;
        for (*info = 1; *info <= i__1; ++(*info)) {
          if (t[*info + *info * t_dim1] == 0.) goto L150;
          /* ......exit */
        }
        *info = 0;
        /* determine the task and go to it. */
        case__ = 1;
        if (*job % 10 != 0) case__ = 2;
        if (*job % 100 / 10 != 0) case__ += 2;

        switch (case__) {
          case 1:  goto L20;
          case 2:  goto L50;
          case 3:  goto L80;
          case 4:  goto L110;
        }

        /* Case 1 (job = 00): */
        /*  solve t*x=b for t lower triangular */
L20:
        b[1] /= t[t_dim1 + 1];
        if (*n >= 2) {
          i__1 = *n;
          for (j = 2; j <= i__1; ++j) {
            temp = -b[j - 1];
            i__2 = *n - j + 1;
            daxpy(&i__2, &temp, &t[j + (j - 1) * t_dim1], &c__1, &b[j], &
                c__1);
            b[j] /= t[j + j * t_dim1];
          }
        }
        goto L140;

        /* Case 2 (job = 01): */
        /*  solve t*x=b for t upper triangular. */
L50:
        b[*n] /= t[*n + *n * t_dim1];
        if (*n >= 2) {
          i__1 = *n;
          for (jj = 2; jj <= i__1; ++jj) {
            j = *n - jj + 1;
            temp = -b[j + 1];
            daxpy(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c__1, &b[1], &c__1);
            b[j] /= t[j + j * t_dim1];
          }
        }
        goto L140;

        /* Case 3 (job = 10): */
        /*  solve trans(t)*x=b for t lower triangular. */
L80:
        b[*n] /= t[*n + *n * t_dim1];
        if (*n >= 2) {
          i__1 = *n;
          for (jj = 2; jj <= i__1; ++jj) {
            j = *n - jj + 1;
            i__2 = jj - 1;
            b[j] -= ddot(&i__2, &t[j + 1 + j * t_dim1], &c__1, &b[j + 1], &c__1);
            b[j] /= t[j + j * t_dim1];
          }
        }
        goto L140;

        /* Case 4 (job = 11): */
        /*  solve trans(t)*x=b for t upper triangular. */
L110:
        b[1] /= t[t_dim1 + 1];
        if (*n >= 2) {
          i__1 = *n;
          for (j = 2; j <= i__1; ++j) {
            i__2 = j - 1;
            b[j] -= ddot(&i__2, &t[j * t_dim1 + 1], &c__1, &b[1], &c__1);
            b[j] /= t[j + j * t_dim1];
          }
        }

L140:
L150:
        return 0;
      } 
      int daxpy(int *n, double *da, double *dx, int *incx, double *dy, int *incy) {
        int i__1;
        int i__, m, ix, iy, mp1;
        /* constant times a vector plus a vector. */
        /* uses unrolled loops for increments equal to one. */
        /* Parameter adjustments */
        --dy;
        --dx;

        if (*n <= 0) return 0;
        if (*da == 0.) return 0;
        if (*incx == 1 && *incy == 1) goto L20;

        /* code for unequal increments or equal increments */
        /*   not equal to 1 */
        ix = 1;
        iy = 1;
        if (*incx < 0) {
          ix = (-(*n) + 1) * *incx + 1;
        }
        if (*incy < 0) {
          iy = (-(*n) + 1) * *incy + 1;
        }
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          dy[iy] += *da * dx[ix];
          ix += *incx;
          iy += *incy;
        }
        return 0;
        /* code for both increments equal to 1 */
        /* clean-up loop */
L20:
        m = *n % 4;
        if (m == 0) goto L40;

        i__1 = m;
        for (i__ = 1; i__ <= i__1; ++i__) {
          dy[i__] += *da * dx[i__];
        }
        if (*n < 4) return 0;

L40:
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 4) {
          dy[i__] += *da * dx[i__];
          dy[i__ + 1] += *da * dx[i__ + 1];
          dy[i__ + 2] += *da * dx[i__ + 2];
          dy[i__ + 3] += *da * dx[i__ + 3];
        }
        return 0;
      }
      int dcopy(int *n, double *dx, int *incx, double *dy, int *incy) {
        int i__1;
        int i__, m, ix, iy, mp1;
        /* copies a vector, x, to a vector, y. */
        /* uses unrolled loops for increments equal to one. */
        /* Parameter adjustments */
        --dy;
        --dx;

        if (*n <= 0) return 0;
        if (*incx == 1 && *incy == 1) goto L20;

        /*  code for unequal increments or equal increments */
        /*    not equal to 1 */
        ix = 1;
        iy = 1;
        if (*incx < 0) {
          ix = (-(*n) + 1) * *incx + 1;
        }
        if (*incy < 0) {
          iy = (-(*n) + 1) * *incy + 1;
        }
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          dy[iy] = dx[ix];
          ix += *incx;
          iy += *incy;
        }
        return 0;
        /* code for both increments equal to 1 */
        /* clean-up loop */
L20:
        m = *n % 7;
        if (m == 0) goto L40;

        i__1 = m;
        for (i__ = 1; i__ <= i__1; ++i__) {
          dy[i__] = dx[i__];
        }
        if (*n < 7) return 0;
L40:
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 7) {
          dy[i__] = dx[i__];
          dy[i__ + 1] = dx[i__ + 1];
          dy[i__ + 2] = dx[i__ + 2];
          dy[i__ + 3] = dx[i__ + 3];
          dy[i__ + 4] = dx[i__ + 4];
          dy[i__ + 5] = dx[i__ + 5];
          dy[i__ + 6] = dx[i__ + 6];
        }
        return 0;
      }
      double ddot(int *n, double *dx, int *incx, double *dy, int *incy) {
        int i__1;
        double ret_val;
        int i__, m, ix, iy, mp1;
        double dtemp;
        /* forms the dot product of two vectors. */
        /* uses unrolled loops for increments equal to one. */
        /* Parameter adjustments */
        --dy;
        --dx;

        ret_val = 0.;
        dtemp = 0.;
        if (*n <= 0) return ret_val;
        if (*incx == 1 && *incy == 1) goto L20;

        /* code for unequal increments or equal increments */
        /* not equal to 1 */
        ix = 1;
        iy = 1;
        if (*incx < 0) {
          ix = (-(*n) + 1) * *incx + 1;
        }
        if (*incy < 0) {
          iy = (-(*n) + 1) * *incy + 1;
        }
        i__1 = *n;
        for (i__ = 1; i__ <= i__1; ++i__) {
          dtemp += dx[ix] * dy[iy];
          ix += *incx;
          iy += *incy;
        }
        ret_val = dtemp;
        return ret_val;
        /* code for both increments equal to 1 */
        /* clean-up loop */
L20:
        m = *n % 5;
        if (m == 0) goto L40;

        i__1 = m;
        for (i__ = 1; i__ <= i__1; ++i__) {
          dtemp += dx[i__] * dy[i__];
        }
        if (*n < 5) goto L60;
L40:
        mp1 = m + 1;
        i__1 = *n;
        for (i__ = mp1; i__ <= i__1; i__ += 5) {
          dtemp = dtemp + dx[i__] * dy[i__] + dx[i__ + 1] * dy[i__ + 1] 
            + dx[i__ + 2] * dy[i__ + 2] + dx[i__ + 3] * dy[i__ + 3] + dx[i__ + 4] * dy[i__ + 4];
        }
L60:
        ret_val = dtemp;
        return ret_val;
      }
      int dscal(int *n, double *da, double *dx, int *incx) {
        int i__1, i__2;
        int i__, m, mp1, nincx;
        /* scales a vector by a constant. */
        /* uses unrolled loops for increment equal to one. */
        /* Parameter adjustments */
        --dx;

        if (*n <= 0 || *incx <= 0) return 0;
        if (*incx == 1) goto L20;

        /* code for increment not equal to 1 */
        nincx = *n * *incx;
        i__1 = nincx;
        i__2 = *incx;
        for (i__ = 1; i__2 < 0 ? i__ >= i__1 : i__ <= i__1; i__ += i__2) {
          dx[i__] = *da * dx[i__];
        }
        return 0;
        /* code for increment equal to 1 */
        /* clean-up loop */
L20:
        m = *n % 5;
        if (m == 0) goto L40;

        i__2 = m;
        for (i__ = 1; i__ <= i__2; ++i__) {
          dx[i__] = *da * dx[i__];
        }
        if (*n < 5) return 0;
L40:
        mp1 = m + 1;
        i__2 = *n;
        for (i__ = mp1; i__ <= i__2; i__ += 5) {
          dx[i__] = *da * dx[i__];
          dx[i__ + 1] = *da * dx[i__ + 1];
          dx[i__ + 2] = *da * dx[i__ + 2];
          dx[i__ + 3] = *da * dx[i__ + 3];
          dx[i__ + 4] = *da * dx[i__ + 4];
        }
        return 0;
      }
  };

}
#endif /* optimizer_h */
