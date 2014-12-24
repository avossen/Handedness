#ifndef STYLE_SETTER_H
#define STYLE_SETTER_H

#include "TStyle.h"
#include <TROOT.h>


void setStyleOpts()
{
  gStyle->SetOptFit(1);
  gROOT->SetStyle("Plain");
  //gStyle->SetOptTitle(kFALSE);
  //  gStyle->SetOptStat(kFALSE);
  gStyle->SetPalette(1,0);

  gStyle->SetCanvasColor(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptTitle(0);

  gROOT->ForceStyle();
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadBottomMargin(0.14);

  gStyle->SetLabelOffset(0.00,"x");
  gStyle->SetTitleOffset(1.,"x");

  gStyle->SetMarkerSize(1.3);
  gStyle->SetLineColor(1);

  gStyle->SetLabelSize(0.06);
  gStyle->SetTitleXSize(0.07);

  gStyle->SetTitleYSize(0.07);

  gStyle->SetTextSize(0.07);

 
}


#endif
