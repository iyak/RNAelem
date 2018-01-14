//
//  struct.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2018/01/05.
//  Copyright © 2018 Kiryu Lab. All rights reserved.
//
#ifndef struct_h
#define struct_h

#ifndef M_PI
#define M_PI (3.1415926535897932384626433832795028841971693993751058)
#endif

#include "util.hpp"
namespace iyak{
  class RNAplot {
    V _angle{};
    VI _ct{};
    int _g=3;/* gap length */
    ostream* _ofs;
    double const _r=15.,_m=10.;
    string static const ps_head,ps_main;
    void build_ct(string const& rss){
      /* ct notation of RNA secondary structure
       * eg: (.(...*)). => 8,1,7,3,4,5,-1,2,0,9 */
      VI stack{};
      for(int i=0;i<size(rss);++i){
        int j=i;
        if('('==rss[i])stack.push_back(i);
        else if(')'==rss[i]){
          j=stack.back();stack.pop_back();
          _ct[j]=i;
        }
        if('*'==rss[i])_ct.push_back(-1);
        else _ct.push_back(j);
      }
    }
    V build_angle(int b,int e){
      V angle{};
      int s=0;while(_ct[b+s]==e-s and _ct[e-s]==b+s)++s;/* calc stem length */
      int i=b+s,npoly=2;while(i<e-s+1){/* calc loop polygon */
        int j=_ct[i];
        npoly+=i<j?2:1;
        i=i<j?j+1:i+1;
      }
      double const extn=2.*M_PI/npoly;
      angle.insert(angle.end(),max(0,s-1),M_PI/2.);
      i=b+s;int k=1;while(i<e-s+1){
        auto j=_ct[i];
        if(i<j){
          angle.push_back(M_PI-extn*k++);/* 1st base dont have angle */
          V sub=build_angle(i,j);
          for(auto& v:sub)v+=M_PI-extn*k;
          angle.insert(angle.end(),sub.begin(),sub.end());
          i=j+1;++k;
        }else{
          angle.push_back(M_PI-extn*k++);
          ++i;
        }
      }
      angle.push_back(M_PI-extn*k++);
      angle.insert(angle.end(),max(0,s-1),-M_PI/2.);
      return angle;
    }
    VV calc_coor(V const& angle){
      double x=0,y=0;
      VV xy{{x,y}};
      for(int i=(size(_ct)-1==_ct[0]?0:1);
          i<size(angle)-(size(_ct)-1==_ct[0]?0:1);
          ++i){
        double xd=_r*cos(angle[i]),yd=_r*sin(angle[i]);
        x+=xd;y+=yd;
        xy.push_back({x,y});
      }
      return xy;
    }
    void make_straight(){
      double dmin=0,wmin=inf;
      for(int d=0;d<90;d+=1){
        for(auto sign:{1,-1}){
          V angle(_angle);
          for(auto& a:angle)a+=M_PI*double(d*sign)/180;
          auto xy=calc_coor(angle);
          double xmin=inf,xmax=-inf;
          for(auto& xyi:xy){
            if(xyi[0]<xmin)xmin=xyi[0];
            if(xmax<xyi[0])xmax=xyi[0];
          }
          if(xmax-xmin<wmin){wmin=xmax-xmin;dmin=d*sign;};
        }
      }
      for(auto& a:_angle)a+=M_PI*double(dmin)/180;
    }
    V calc_bbox(VV const& xy){
      double xmin=inf,xmax=-inf,ymin=inf,ymax=-inf;
      for(auto const& xyi:xy){
        if(xyi[0]<xmin)xmin=xyi[0];
        if(xmax<xyi[0])xmax=xyi[0];
        if(xyi[1]<ymin)ymin=xyi[1];
        if(ymax<xyi[1])ymax=xyi[1];
      }
      return {xmin,ymin,xmax,ymax};
    }
    void set_seq_rss(string& seq,string& rss){
      seq.resize(size(rss),' ');
      /* remove *s at both ends */
      auto b=rss.find_first_not_of("*"),e=rss.find_last_not_of("*");
      rss=rss.substr(b,e-b+1);seq=seq.substr(b,e-b+1);
      /* remove adjascent *s */
      for(int i=1;i<size(rss);++i)
        if('*'==rss[i-1] and '*'==rss[i]){rss.erase(i,1);seq.erase(i,1);--i;}
      /* triple *s */
      for(int i=0;i<size(rss);++i){
        if('*'==rss[i]){
          rss.insert(i,string(_g-1,'*'));
          seq.insert(i,string(_g-1,' '));
          i+=_g-1;
        }
      }
    }
  public:
    void set_ostream(ostream& ofs){_ofs=&ofs;}
    void no_plot(){
      (*_ofs)<<
      "%!PS-Adobe-3.0 EPSF-3.0\n"
      "%%BoundingBox:0 0 62 62\n"
      "/Courier findfont 8 scalefont setfont\n"
      "2 28 moveto (no structure) show showpage\n";
    }
    void plot(string seq,string rss){
      set_seq_rss(seq,rss);
      build_ct(rss);
      _angle=build_angle(0,size(rss)-1);
      make_straight();
      auto xy=calc_coor(_angle);
      auto bbox=calc_bbox(xy);
      (*_ofs)<<"%!PS-Adobe-3.0 EPSF-3.0\n%%BoundingBox:"
      <<paste1(bbox[0]-_m,bbox[1]-_m,bbox[2]+_m,bbox[3]+_m)<<"\n"
      <<ps_head;
      (*_ofs)<<"/xy [\n";
      for(auto& xyi:xy)(*_ofs)<<"["<<paste(xyi," ")<<"]\n";
      (*_ofs)<<"] def\n";
      (*_ofs)<<"/pairs [\n";
      for(int i=0;i<size(_ct);++i)(*_ofs)<<(i<_ct[i]?"["+paste1(i,_ct[i])+"]\n":"");
      (*_ofs)<<"] def\n";
      (*_ofs)<<"/unpairs [";
      for(int i=0;i<size(_ct);++i)(*_ofs)<<(i==_ct[i]?to_str(i)+" ":"");
      (*_ofs)<<"] def\n";
      (*_ofs)<<"/gaps [";
      for(int i=size(_ct)-1;0<=i;--i)if(-1==_ct[i]){(*_ofs)<<i-2<<" ";i-=2;}
      (*_ofs)<<"] def\n";
      (*_ofs)<<"/seq ("<<seq<<") def\n";
      (*_ofs)<<ps_main;
    }
  };
  string const RNAplot::ps_head=
  "/mapa {\n"
  "  6 dict begin\n"
  "  /f exch def /a1 exch def /n a1 length def /a2 n array def\n"
  "  0 1 n 1 sub {/i exch def /e a1 i get def a2 i e f put} for\n"
  "  a2 end\n"
  "} bind def\n"
  "/drawpairs {\n"
  "  gsave newpath 0.7 setgray 3 setlinewidth [7 4.01] 7 setdash\n"
  "  pairs {aload pop xy exch get aload pop moveto xy exch get aload pop lineto} forall\n"
  "  stroke grestore\n"
  "} bind def\n"
  "/fxy {exch xy exch get aload pop 3 -1 roll exec} bind def % int proc fxy -\n"
  "/drawbackbone {\n"
  "  gsave newpath 0.2 setgray 0.8 setlinewidth\n"
  "  0 {moveto} fxy -1 gaps aload pop 1 {2 copy eq {3 add {moveto} fxy 4 add}\n"
  "    {dup {lineto} fxy 1 add} ifelse dup xy length ge {exit} if} loop\n"
  "  stroke grestore\n"
  "} bind def\n"
  "/drawgaps {\n"
  "  gsave newpath 0.2 setgray 0.8 setlinewidth [1 3] 0 setdash\n"
  "  gaps {1 sub dup {moveto} fxy 1 add dup {lineto} fxy\n"
  "    2 add dup {moveto} fxy 1 add {lineto} fxy} forall\n"
  "  stroke grestore\n"
  "} bind def\n"
  "/drawloop {\n"
  "  gsave newpath 0 setgray 0.8 setlinewidth\n"
  "  unpairs {xy exch get aload pop 2 copy moveto 0.5 0 360 arc} forall\n"
  "  stroke grestore\n"
  "} bind def\n"
  "/drawbases {\n"
  "  gsave newpath 0 setgray 0.8 setlinewidth\n"
  "  0 {moveto} fxy 0 {4 0 360 arc gsave stroke grestore fill} fxy\n"
  "  0 xy {{2.3 sub} mapa aload pop moveto dup seq exch 1 getinterval false charpath\n"
  "    gsave 2.5 setlinewidth 1 setgray stroke grestore fill 1 add} forall\n"
  "  stroke grestore\n"
  "} bind def\n";
  string const RNAplot::ps_main=
  "/Courier findfont 8 scalefont setfont 1 setlinejoin 1 setlinecap\n"
  "drawbackbone drawgaps drawpairs drawloop drawbases showpage\n";
}

#endif
