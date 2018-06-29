//
//  logo.hpp
//  RNAelem
//
//  Created by Hiroshi Miyake on 2017/01/31.
//  Copyright © 2017 Kiryu Lab. All rights reserved.
//

#ifndef logo_h
#define logo_h

#include<iostream>
#include<sstream>
#include<string>
#include<vector>
#include<cmath>
#include<unordered_map>
#include<ft2build.h>
#include"util.hpp"
#include FT_FREETYPE_H
#include FT_OUTLINE_H
#include FT_BBOX_H

#ifndef DEFAULT_FONT
#define DEFAULT_FONT "You must manually pass font file (.ttf) path"
#endif

namespace iyak {
  umap<uint32_t, string> _color_map {
    /* ascii can replace utf8 */
    {'A',"#339541"}, {'a',"#339541"},
    {'G',"#F5C000"}, {'g',"#F5C000"},
    {'C',"#545FFF"}, {'c',"#545FFF"},
    {'U',"#D21010"}, {'u',"#D21010"},
    {'T',"#D21010"}, {'t',"#D21010"},
    {'B',"#BEAED4"},
    {'O',"#FFFF99"},
    {'H',"#F0027F"},
    {'I',"#FDC086"},
    {'M',"#7FC97F"},
    {'L',"#386CB0"}, {'(',"#386CB0"},
    {'R',"#386CB0"}, {')',"#386CB0"},
  };

  class utf8 {
    vector<uint32_t> _codes;
    string _s;
    VS _vs{};
    int nunits(unsigned char byte) {
      return
      byte<0x80? 1:
      byte<0xe0? 2:
      byte<0xf0? 3:
      byte<0xf8? 4:
      -1;
    }
    uint32_t unit(string s) {
      check(size(s)==nunits(s[0]), "wrong unit:", s);
      if (1==size(s)) return s[0];
      uint32_t c = (s[0]&(0xff>>(size(s)+1))) << (6*(size(s)-1));
      for (int i=1; i<size(s); ++i) {
        c |= (s[i]&0x3f) << (6*(size(s)-1-i));
      }
      return c;
    }
  public:
    int len() {return size(_codes);}
    string str() {return _s;}
    string str(int i){return _vs[i];}
    uint32_t code(int i) {return _codes.at(i);}
    vector<uint32_t> codes() {return _codes;}
    utf8(string s) {
      _s = s;
      for (int i=0; i<size(s);) {
        int n = nunits(s[i]);
        _vs.push_back(s.substr(i,n));
        _codes.push_back(unit(_vs.back()));
        i += n;
      }
    }
  };

  class ft_face {
    FT_Face _face = nullptr;
    FT_Library _lib = nullptr;
  public:
    ft_face(string font) {
      auto e=FT_Init_FreeType(&_lib);
      check(!e, "fail init lib");
      e = FT_New_Face(_lib, font.c_str(), 0, &_face);
      check(!e, "fail init lib:", font);
    }
    ~ft_face() {FT_Done_Face(_face);FT_Done_FreeType(_lib);}
    FT_Face get() {return _face;}
  };

  class RNAlogoAlph {
    string _alph;
    string _font;
    string _color;
    vector<uptr<ft_face>> _faces;

    int _h = 0;
    int _w = 0;

    void load_glyph(ft_face& f, FT_ULong u, int bold=25) {
      auto i = FT_Get_Char_Index(f.get(), u);
      auto err = FT_Load_Glyph(f.get(), i, FT_LOAD_NO_SCALE|FT_LOAD_NO_BITMAP);
      check(!err, "fail load glyph", u);
      err = FT_Outline_Embolden(&o(f), bold);
      check(!err, "fail embolden:", _alph);
      mapply(f, {{1,0},{0,-1}}); /* flip */
      double rx = double(m(f).horiBearingX) / m(f).width;
      double ry = double(m(f).horiBearingY) / m(f).height;
      mapply(f, {{rx,0},{0,ry}}); /* shrink */
    }
    void mapply (ft_face& f, VV m) {
      FT_Fixed scaler = 1<<16;
      FT_Matrix mm {
        .xx=static_cast<FT_Fixed>(round(m[0][0]*scaler)),
        .xy=static_cast<FT_Fixed>(round(m[0][1]*scaler)),
        .yx=static_cast<FT_Fixed>(round(m[1][0]*scaler)),
        .yy=static_cast<FT_Fixed>(round(m[1][1]*scaler)),
      };
      FT_Outline_Transform(&o(f), &mm);
    }
    struct bbox {
      int x,y,w,h;
      bbox(FT_Pos a,FT_Pos b,FT_Pos c,FT_Pos d):
      x((int)a),y((int)b),w((int)c),h((int)d){}
    };
    bbox calc_bb (ft_face& f) {
      FT_BBox bb;
      auto err = FT_Outline_Get_BBox(&o(f), &bb);
      check(!err, "fail get bb:", _alph);
      return bbox(bb.xMin,bb.yMin,bb.xMax-bb.xMin,bb.yMax-bb.yMin);
    }
    void shift_bb(ft_face& f, int x, int y) {
      FT_Outline_Translate(&o(f), x, y);
    }
    void set_pos(ft_face& f, int x, int y) {
      auto bb = calc_bb(f);
      shift_bb(f, x-bb.x, y-bb.y);
      bb = calc_bb(f);
    }
    FT_Outline& o(ft_face& f) {return f.get()->glyph->outline;}
    FT_Glyph_Metrics& m(ft_face& f) {return f.get()->glyph->metrics;}
  public:
    RNAlogoAlph(string const& a,double x,double y,double w,double h,
                string const& color="black",string const& font=DEFAULT_FONT,
                int bold=25)
    :_alph(a),_font(font),_color(color){
      auto utf = utf8(_alph).codes();
      _h=0;
      _w=999*size(utf);
      for (int i=0; i<size(utf); ++i) {
        _faces.push_back(uptrize(new ft_face(_font)));
        auto& f = _faces.back();
        load_glyph(*f, utf[i], bold);
        auto bb = calc_bb(*f);
        mapply(*f, {{999./bb.w,0.},{0.,1.}});
        _h = max(_h, bb.h);
        set_pos(*f, 999*i, 0);
      }
      for (auto& f: _faces) {
        auto bb = calc_bb(*f);
        shift_bb(*f, 0, _h-bb.h);
        mapply(*f, {{w/_w,0.},{0.,h/_h}});
        shift_bb(*f, x, y);
      }
    }
    friend ostream& operator<<(ostream&,RNAlogoAlph const&);
  };
  ostream& operator<<(ostream& os,RNAlogoAlph const& alph){
    for (int i=0; i<size(alph._faces); ++i){
      auto& f = alph._faces[i];
      os<<"<path d=\"\n";
      FT_Outline_Funcs callback{
        .move_to=[](FT_Vector const* to,
                    void* ptr){
          *(ostream*)ptr<<paste1("M",to->x,to->y,"\n");
          return 0;
        },
        .line_to=[](FT_Vector const* to,
                    void* ptr){
          *(ostream*)ptr<<paste1("L",to->x,to->y,"\n");
          return 0;
        },
        .conic_to=[](FT_Vector const* c,
                     FT_Vector const* to,
                     void* ptr){
          *(ostream*)ptr<<paste1("Q",c->x,c->y,",",to->x,to->y,"\n");
          return 0;
        },
        .cubic_to=[](FT_Vector const* c0,
                     FT_Vector const* c1,
                     FT_Vector const* to,
                     void* ptr){
          *(ostream*)ptr<<paste1("C",c0->x,c0->y,",",c1->x,c1->y,
                                 ",",to->x,to->y,"\n");
          return 0;
        }
      };
      FT_Outline_Decompose(&(f->get()->glyph->outline), &callback,&os);
      os<<"\" fill=\""+alph._color+"\" />\n";
    }
    return os;
  }
  struct block{
    string chrs;
    vector<RNAlogoAlph> chr_paths;
    double val;
    VS colors/* len(chrs)==len(colors) */;
    string font=DEFAULT_FONT;
    block(string s):chrs(s){
      auto utf=utf8(s).codes();
      for(auto u:utf)colors.push_back(_color_map[u]);
    }
  };
  struct column{
    vector<block> blocks;
    string meta="";
    string meta_font=DEFAULT_FONT;
  };
  struct logo_data{
    vector<column> columns;
    string title="";
    double vmax=2.;
  };

  class RNAlogo {
    logo_data _data;
    ostream* _out;
    ostream& out(){return *_out;}

    int _colw = 100;
    int _space = 5;
    int _rowh = 500;
    int _titleh = 50;
    int _yaxisw = 50;
    int _yrulerl = 10;
    int _xaxish = 50;
    int _metah = 50;

    double _scale;
    int spaced(int a) {return a;}
    template<class...T>
    int spaced(int a,T...b) {return a + spaced(b...) + (0<a? _space:0);}

    void put_svg_viewbox() {
      out()<<paste1
      (0,0,
       spaced(_yaxisw, _yrulerl, (_colw+_space)*size(_data.columns)),
       spaced(_titleh, _rowh, _metah, _xaxish));
    }

    void put_svg_path() {
      for (auto& column:_data.columns)
        for (auto& block:column.blocks)
          for(auto& path:block.chr_paths)
            out()<<path;
    }

    void put_svg_axes() {
      if (0<_yaxisw) {
        /* y axis */
        out()<<paste1
        (
         "<path d=\"M",
         spaced(_yaxisw, _yrulerl),
         spaced(_titleh, _rowh - _scale*(int(_data.vmax))),
         "L",
         spaced(_yaxisw, _yrulerl),
         spaced(_titleh, _rowh),
         "\" stroke=\"black\" stroke-width=\"1.5\"/>\n"
         );

        for (int i=0; i<=_data.vmax+1e-10; ++i) {
          /* y ruler */
          out()<<paste1
          (
           "<path d=\"M",
           _yaxisw + _space,
           spaced(_titleh, _rowh-round(_scale*i)),
           "L",
           spaced(_yaxisw, _yrulerl),
           spaced(_titleh, _rowh-round(_scale*i)),
           "\" stroke=\"black\" stroke-width=\"4\"/>\n"
           );
          /* y ruler text */
          out()<<paste0
          (
           "<text text-anchor=\"end\" x=\"", _yaxisw,
           "\" y=\"", spaced(_titleh, _rowh-round(_scale*i)),
           "\" font-size=\"", _yaxisw,
           "\" alignment-baseline=\"middle",
           "\">", to_str(i), "</text>\n"
           );
        }
      }
      if (0<_xaxish) {
        /* x axis text */
        for (int i=0; i<size(_data.columns); ++i) {
          out()<<paste0
          (
           "<text text-anchor=\"middle\" x=\"",
           spaced(_yaxisw, _yrulerl, (_colw+_space)*i + _colw/2),
           "\" y=\"", spaced(_titleh, _rowh, _metah, _xaxish/2),
           "\" font-size=\"", _xaxish,
           "\" alignment-baseline=\"baseline"
           "\">", to_str(i+1), "</text>\n"
           );
        }
      }
    }

    void put_svg_title() {
      if (0<_titleh) {
        out()<<paste0
        (
         "<text text-anchor=\"middle\" x=\"",
         spaced(_yaxisw, _yrulerl, ((_colw+_space)*size(_data.columns)-_space)/2),
         "\" y=\"", _titleh,
         "\" font-size=\"", _titleh,
         "\" alignment-baseline=\"baseline",
         "\">", _data.title, "</text>\n"
         );
      }
    }

    void put_svg_meta(){
      if(0<_metah){
        for(int i=0; i<size(_data.columns); ++i){
          auto utf=utf8(_data.columns[i].meta);
          auto w=double(_colw)/utf.len();
          for(int j=0;j<utf.len();++j){
            out()<<
            RNAlogoAlph(utf.str(j),
                        spaced(_yaxisw,_yrulerl,(_colw+_space)*i)+(w*j)+w*0.1,
                        spaced(_titleh,_rowh)+_space+_metah*0.1,
                        w*0.8,_metah*0.8,"black",_data.columns[i].meta_font);
          }
        }
      }
    }

    void put_svg() {
      string s=
      "<svg xmlns=\"http://www.w3.org/2000/svg\"\n"
      "viewBox=\"$VIEWBOX\">\n"
      "$PATHTAG"
      "$AXESTAG"
      "$TITLETAG"
      "$METATAG"
      "</svg>\n";
      size_t i=0;
      while(i<s.size()){
        if('$'==s[i]){
            if("$VIEWBOX"==s.substr(i,8)){put_svg_viewbox();i+=8;}
            else if("$PATHTAG"==s.substr(i,8)){put_svg_path();i+=8;}
            else if("$AXESTAG"==s.substr(i,8)){put_svg_axes();i+=8;}
            else if("$METATAG"==s.substr(i,8)){put_svg_meta();i+=8;}
            else if("$TITLETAG"==s.substr(i,9)){put_svg_title();i+=8;}
            else{out()<<s[i++];}
        }
        else{out()<<s[i++];}
      }
    }

    V v2bit (V const& w) {
      V v(w);
      double sum = 0;
      for (auto vv: v) sum += vv;
      if (0!=sum) for (auto& vv: v) vv /= sum;
      double signal = log2(size(v));
      for (auto vv: v) if (0!=vv) signal += vv * log2(vv);
      for (auto& vv: v) vv = vv * signal;
      return v;
    }
  public:
    /* setters */
    void set_column_width(int column_width) {_colw = column_width;}
    void set_row_height(int row_height) {_rowh = row_height;}
    void set_title_height(int title_height) {_titleh = title_height;}
    void set_x_axis_height(int x_axis_height) {_xaxish = x_axis_height;}
    void set_y_axis_width(int y_axis_width) {_yaxisw = y_axis_width;}
    void set_meta_height(int meta_height) {_metah = meta_height;}
    void set_space(int space) {_space = space;}

    logo_data& data(){return _data;}
    void map_color(std::string s, string c) {
      auto utf = utf8(s).codes();
      for (auto u: utf) _color_map[u] = c;
    }
    void plot(){
      _data.vmax=0;
      for(auto& c:_data.columns)
        if(0<size(c.blocks) and _data.vmax<log2(size(c.blocks)))
          _data.vmax=log2(size(c.blocks));
      _scale=_rowh/_data.vmax;
      _titleh=max(_titleh,_yaxisw);
      auto x=spaced(_yaxisw,_yrulerl)+_space;
      for(int i=0;i<size(_data.columns);++i){
        auto& c=_data.columns[i];
        std::sort(c.blocks.begin(),c.blocks.end(),
                  [](block const& x,block const& y){return x.val<y.val;});
        V vals{};for(auto& b:c.blocks)vals.push_back(b.val);
        auto vv=v2bit(vals);
        auto y=spaced(_titleh,_rowh);
        for(int j=0;j<size(c.blocks);++j){
          auto& b=c.blocks[j];
          b.val=vv[j];
          y-=b.val*_scale;
          auto utf=utf8(b.chrs);
          for(int k=0;k<utf.len();++k){
            double w=_colw/utf.len();
            b.chr_paths.emplace_back(utf.str(k),x+k*w,y,w,b.val*_scale,
                                     b.colors[k],b.font,100);
          }
        }
        x+=_colw+_space;
      }
      put_svg();
    }
    void set_ostream(ostream& out){_out=&out;}
  };
}

#endif /* logo_h */
