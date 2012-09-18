#ifndef dpdxGraph_H
#define dpdxGraph_H

#include <TGraph.h>
#include <TLegend.h>

class dpdxGraph {
 public:

  dpdxGraph(const char*, int);
  ~dpdxGraph();

  void addPoint(double&, double&);
  void draw(TLegend*);

 private:

  TGraph * theGraph;
  std::vector<double> xVec;  
  std::vector<double> yVec;  

  TString TSname;
  int iCol;

};

#endif 
