#include <iostream>
#include "dpdxGraph.h"

dpdxGraph::dpdxGraph(const char* name, int icol){

  std::cout << " >>>>>>>>>>><<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<Constructing dpdxGraph:" << name << std::endl; 

  xVec.clear();
  yVec.clear();
  TSname.Append(name);
  iCol = icol;

}
  
dpdxGraph::~dpdxGraph()
{

  std::cout << " Cleaning dpdxGraph object..." << std::endl;

}

void dpdxGraph::addPoint(double& x, double& y){

  //  std::cout << " adding " << x << " " << y << " " << xVec.size() << std::endl;

  xVec.push_back(x);
  yVec.push_back(y);

}

void dpdxGraph::draw(TLegend *leg){


  double* xArray = &xVec[0];
  double* yArray = &yVec[0];
  theGraph = new TGraph(xVec.size(), xArray, yArray);

  theGraph->SetMarkerStyle(20);
  theGraph->SetMarkerSize(0.7);
  theGraph->SetMarkerColor(iCol);
  theGraph->Draw("psame");

  leg->AddEntry(theGraph,TSname.Data(),"p");

}
