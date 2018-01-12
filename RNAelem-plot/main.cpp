//
//  main.cpp
//  RNAelem-plot
//
//  Created by Hiroshi Miyake on 2018/01/12.
//  Copyright © 2018 Hiroshi Miyake. All rights reserved.
//
#include"util.hpp"
#include"struct.hpp"

int main(int argc,char *argv[]){
  iyak::RNAplot plot;
  plot.set_ostream(std::cout);
  try{plot.plot(argv[1]/*seq*/,argv[2]/*rss*/);}
  catch(...){plot.no_plot();}
  return 0;
}
