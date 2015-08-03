#ifndef configure_histos_h
#define onfigure_histos_h

#include <iostream>
using namespace std;

#include <TRandom3.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMath.h>
#include <TPad.h>
#include <TPaveLabel.h>
#include <TObjArray.h>
#include <TLorentzVector.h>
#include <TChain.h>
#include <TPaveText.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TGraphAsymmErrors.h>
#include <TVectorD.h>
#include <TColor.h>
#include <TBox.h>
#include <TLine.h>
#include <TMarker.h>
#include <TStopwatch.h>
#include <TMinuit.h>
#include <TArrow.h>





void drawcutborders(TH1 *gr,Double_t x1, Int_t x1_dir)
{
 
  
  Double_t x[2], y[2];
  x[0] = x1;
  x[1] = x1;
  y[0] = gr->GetMinimum()*0.01;
  y[1] = gr->GetMaximum()*10.0;
  TGraph *bla = new TGraph(2, x, y);
  bla->SetLineStyle(1);
  bla->SetLineWidth(1.75);
  bla->SetLineColor(46);
  
  bla->Draw("L same");

  Double_t arrow_y_pos = (gr->GetMaximum() - gr->GetMinimum())*0.5 + gr->GetMinimum();
  Double_t delta_x = 0.0;
  if(gPad->GetLogx() == 0)
    {
      delta_x = 0.04*TMath::Abs(gr->GetBinLowEdge(gr->GetXaxis()->GetLast()) +  gr->GetBinWidth(gr->GetXaxis()->GetLast()) - gr->GetBinLowEdge(gr->GetXaxis()->GetFirst()));
      if(x1_dir<0)
	delta_x = delta_x * -1.0;
 
    }
  else
    {
      if(x1_dir<0)
	delta_x =x1 * 0.8 - x1;
      else
	delta_x =x1 * 1.2 - x1;

    }
  
 
  TArrow *arrow = new TArrow(x1, arrow_y_pos, x1+delta_x, arrow_y_pos, 0.01, ">");
  arrow->SetLineColor(46);
  arrow->SetLineWidth(1.75);
  arrow->Draw("same >");
 
}




Double_t getmax(TH1D *h1, TH1D *h2, TH1D *h3)
{
  
  if(h1->GetMaximum() > h2->GetMaximum() && h1->GetMaximum() > h3->GetMaximum())
    return h1->GetMaximum();
  if(h2->GetMaximum() > h3->GetMaximum() && h2->GetMaximum() > h1->GetMaximum())
    return h2->GetMaximum();
  if(h3->GetMaximum() > h1->GetMaximum() && h3->GetMaximum() > h2->GetMaximum())
    return h3->GetMaximum();
 
return 0.; 
}
Double_t getmin(TH1D *h1, TH1D *h2, TH1D *h3)
{
  
  if(h1->GetMaximum() <= h2->GetMaximum() && h1->GetMaximum() <= h3->GetMaximum())
    return h1->GetMinimum();
  if(h2->GetMaximum() < h3->GetMaximum() && h2->GetMaximum() < h1->GetMaximum())
    return h2->GetMinimum();
  if(h3->GetMaximum() < h1->GetMaximum() && h3->GetMaximum() < h2->GetMaximum())
    return h3->GetMinimum();
  
return 0.;
}





TString outputdirectory = "/data/zenhamburg1b/schorner/jbehr/results";
void pad_cfg(Int_t maxd = 2, //2
	     Short_t borderMode = 0, //0
	     Double_t leftMargin = 0.16, //0.12
	     Double_t rightMargin = 0.05,//0.05
	     Double_t bottomMargin = 0.12)//0.12
{
  gStyle->SetPalette(1,0);
  gPad->SetTickx(1);
  gPad->SetTicky(1);
  gPad->SetBorderMode(borderMode);
  gPad->SetLeftMargin(leftMargin);
  gPad->SetRightMargin(rightMargin);
  gPad->SetBottomMargin(0.15);
  gPad->SetTopMargin(0.07); //0.001 //0.1
  TGaxis::SetMaxDigits(maxd);

  maxd = 2; //2
  borderMode = 0; //0
  leftMargin = 0.16; //0.12
  rightMargin = 0.05;//0.05
  bottomMargin = 0.12;//0.12
}

void drawunityline(TH1 *gr)
{
  Double_t x[2], y[2];
  x[0] = gr->GetXaxis()->GetXmin();
  x[1] = gr->GetXaxis()->GetXmax();
  y[0] = 1.0;
  y[1] = 1.0;
  TGraph *bla = new TGraph(2, x, y);
  bla->SetLineStyle(1);
  bla->SetLineWidth(1);
  bla->SetLineColor(kBlack);
  bla->Draw("L same");
}


void drawvalueline(TH1 *gr, Double_t v)
{
  Double_t x[2], y[2];
  x[0] = gr->GetXaxis()->GetXmin();
  x[1] = gr->GetXaxis()->GetXmax();
  y[0] = v;
  y[1] = v;
  TGraph *bla = new TGraph(2, x, y);
  bla->SetLineStyle(1);
  bla->SetLineWidth(3);
  bla->SetLineColor(kBlack);
  bla->Draw("L same");
}




void drawverticalunityline(TH1 *gr)
{
  Double_t x[2], y[2];
  x[0] = 1.0;
  x[1] = 1.0;
  y[0] = 0.0;
  y[1] = gr->GetMaximum()*2.0;
  TGraph *bla = new TGraph(2, x, y);
  bla->SetLineStyle(1);
  bla->SetLineWidth(1);
  bla->Draw("L same");
}

void pad_max()
{
  //
  // set correct maximum for y-axis in a pad with superimposed histograms
  //
/*   TList *fPrimitives = gPad->GetListOfPrimitives(); */
/*   TIter iter(fPrimitives); */
/*   TObject *obj; */
/*   TH1 *h_1; */
/*   Bool_t foundfirstHisto = kFALSE; */
/*   while ((obj = (TObject *)iter.Next())) { */
/*     if(obj->InheritsFrom("TH1")) { */
/*       if(foundfirstHisto == kFALSE) { */
/* 	h_1 = (TH1*)obj; */
/* 	foundfirstHisto = kTRUE; */
/*       } */
/*       else { */
/* 	TH1 *h_i; */
/* 	h_i = (TH1*)obj; */
/* 	Double_t max_i = h_i->GetMaximum(); */
/* 	if(max_i > h_1->GetMaximum()) h_1->SetMaximum(1.1 * max_i); */
/*       } */
/*     } */
/*   } */
}

void histo_cfg(TH1 *thisHisto, const char* xTitle, const char* yTitle,
	       const char* title = "")
{
  //thisHisto->SetLineWidth(1.0);
  pad_cfg();
  thisHisto->SetTitle(title);
  gStyle->SetTitleX(.2); //.15
  gStyle->SetTitleW(.95);
  gStyle->SetTitleFontSize(.035);//0.04
  gStyle->SetTitleStyle(0);
  //gStyle->SetTitleAlign(22);
  gStyle->SetTitleBorderSize(0);
  thisHisto->GetXaxis()->SetNoExponent();
  
  //thisHisto->GetYaxis()->SetNoExponent();
  thisHisto->GetXaxis()->SetLabelSize(0.07);
  thisHisto->GetXaxis()->SetTitle(xTitle);

  thisHisto->GetXaxis()->SetTitleSize(0.07); //was 0.1
  thisHisto->GetXaxis()->SetTitleOffset(0.9); //.83

  

  //thisHisto->GetXaxis()->CenterTitle();
  thisHisto->GetYaxis()->SetLabelSize(0.07);
  thisHisto->GetYaxis()->SetTitle(yTitle);
 /*  thisHisto->GetYaxis()->SetTitleOffset(1.0); //1.2 //0.7 */
  thisHisto->GetYaxis()->SetTitleOffset(0.9); //0.74
  thisHisto->GetYaxis()->SetTitleSize(0.07);
  //thisHisto->GetYaxis()->CenterTitle();
  thisHisto->SetStats(kFALSE);
  thisHisto->GetYaxis()->SetNdivisions(505);
  thisHisto->GetXaxis()->SetNdivisions(505);

}


void histo_cfg_subpads(TH1 *thisHisto, const char* xTitle, const char* yTitle,
	       const char* title = "")
{
  //thisHisto->SetLineWidth(1.0);
  
  thisHisto->SetTitle(title);
  gStyle->SetTitleX(.15);
  gStyle->SetTitleW(.95);
  gStyle->SetTitleFontSize(.04);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  thisHisto->GetXaxis()->SetNoExponent();
  //thisHisto->GetYaxis()->SetNoExponent();
  thisHisto->GetXaxis()->SetLabelSize(0.03);
  thisHisto->GetXaxis()->SetTitle(xTitle);
  thisHisto->GetXaxis()->SetTitleOffset(0.7); //1.0

  thisHisto->GetXaxis()->SetTitleSize(0.08);

  //thisHisto->GetXaxis()->CenterTitle();
  thisHisto->GetYaxis()->SetLabelSize(0.03);
  thisHisto->GetYaxis()->SetTitle(yTitle);
  thisHisto->GetYaxis()->SetTitleOffset(0.7); //1.2 //0.7
  thisHisto->GetYaxis()->SetTitleSize(0.045);
  //thisHisto->GetYaxis()->CenterTitle();
  thisHisto->SetStats(kFALSE);
  thisHisto->GetYaxis()->SetNdivisions(505);
  thisHisto->GetXaxis()->SetNdivisions(505);

}

void tprofile_cfg(TProfile *thisHisto, const char* xTitle, const char* yTitle,
	       const char* title = "")
{
  pad_cfg();
  thisHisto->SetTitle(title);
  gStyle->SetTitleX(.15);
  gStyle->SetTitleW(.95);
  gStyle->SetTitleFontSize(.04);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleBorderSize(0);
  thisHisto->GetXaxis()->SetLabelSize(0.03);
  thisHisto->GetXaxis()->SetTitle(xTitle);
  thisHisto->GetXaxis()->SetTitleOffset(0.7); //1.0
  thisHisto->GetXaxis()->SetTitleSize(0.06);
  //thisHisto->GetXaxis()->CenterTitle();
  thisHisto->GetYaxis()->SetLabelSize(0.03);
  thisHisto->GetYaxis()->SetTitle(yTitle);
  thisHisto->GetYaxis()->SetTitleOffset(0.7); //1.2
  thisHisto->GetYaxis()->SetTitleSize(0.06);
  //thisHisto->GetYaxis()->CenterTitle();
  thisHisto->SetStats(kFALSE);
}

void leg_cfg(TLegend *thisLegend)
{
  thisLegend->SetFillStyle(0);
  thisLegend->SetBorderSize(0);
}

void histo_norm_1(TH1F *thisHisto)
{
  Double_t scale = 1./thisHisto->Integral();
  thisHisto->Scale(scale);
}

void histo_norm_bincontent(TH1 *thisHisto)
{
  Double_t width = thisHisto->GetBinWidth(1);
  thisHisto->Scale(1./width);
}

void histo_norm_entries(TH1 *thisHisto)
{
  Double_t entries = thisHisto->GetEntries();
  thisHisto->Scale(1/entries);
}

void histo_norm_entries_and_bin(TH1 *thisHisto)
{
  histo_norm_entries(thisHisto);
  histo_norm_bincontent(thisHisto);
}

void bisect_draw(TH1 *thisHisto, Style_t lstyle = 2, Width_t lwidth = 1)
{
  //
  // draw a bisecting line on a given histogram
  //
  Double_t x[2], y[2];
  x[0] = thisHisto->GetXaxis()->GetXmin();
  x[1] = thisHisto->GetXaxis()->GetXmax();
  y[0] = thisHisto->GetMinimum();
  y[1] = thisHisto->GetMaximum();
  TGraph *gr = new TGraph(2, x, y);
  gr->SetLineStyle(lstyle);
  gr->SetLineWidth(lwidth);
  gr->Draw("L same");
}

void drawzeroline(TH1 *gr)
{
  Double_t x[2], y[2];
  x[0] = gr->GetXaxis()->GetXmin();
  x[1] = gr->GetXaxis()->GetXmax();
  y[0] = 0.0;
  y[1] = 0.0;
  TGraph *bla = new TGraph(2, x, y);
  bla->SetLineStyle(1);
  bla->SetLineWidth(1);
  bla->Draw("L same");
}

void shifthisto(TH1* histo, TString drawoptions, Double_t mul = 1.0, Double_t add = 0.0 )
{
  //histo->SetMarkerSize(0.6);
  
  //Int_t xfirst = histo->GetXaxis()->GetFirst();
  Int_t xlast = histo->GetXaxis()->GetLast();

  TGraphErrors *gr;
  gr = new TGraphErrors(histo);
  
 

  for(Int_t i = 0; i < xlast; i++)
    {
      Double_t x = histo->GetBinCenter(i+1) * mul + add;
      Double_t xerror = 0.0;
      Double_t y = histo->GetBinContent(i+1);
      Double_t yerror = histo->GetBinError(i+1);
      gr->RemovePoint(i);
      if(y==0.0)
	{
	  gr->SetPoint(i, x, -99999.9);
	  gr->SetPointError(i, xerror, yerror);
	}
      else
	{
	  gr->SetPoint(i, x, y);
	  gr->SetPointError(i, xerror, yerror);
    	}
      
	
    }
  gr->SetMarkerSize(1.3);
  gr->Draw(drawoptions);
 
}


//void ratio_plot(TH1 *hist1, TH1 *hist2, TString x_Axis, Float_t lowli, Float_t upli)
void ratio_plot( TH1D *hist1, TH1D *hist2, TString x_Axis, TString y_Axis)
{
  gPad->SetFrameBorderMode(0);

  Int_t logx = gPad->GetLogx();
   gPad->SetTickx(1);
  gPad->SetTicky(1);

  gPad->SetLogx(0);
  // gPad->SetLogy(0);

  TH1D* ratio = (TH1D*)hist1->Clone(); 
  ratio->Divide(hist1,hist2, 1.0, 1.0);
  ratio->SetLineColor(kBlack);



  gPad->SetFillStyle(0);
  gPad->SetFillColor(0);
 
  hist1->SetTitle("");
  hist2->SetTitle("");
  
  hist1->SetStats(kFALSE);
  hist1->GetXaxis()->SetLabelSize(0);
  hist1->GetXaxis()->SetTitleSize(0);
  hist1->GetYaxis()->SetLabelSize(0.04);
  hist1->GetYaxis()->SetLabelOffset(0.01);
  hist1->GetYaxis()->SetTitleSize(0.03);
  //hist1->GetYaxis()->SetTitleOffset(0.001);
  
  hist1->GetXaxis()->SetNdivisions(505);
  hist1->GetYaxis()->SetTitle(y_Axis);
  hist1->GetXaxis()->SetNoExponent();
  //hist1->GetYaxis()->SetNoExponent();
  hist1->GetYaxis()->SetTitleOffset(1.5);
  
  hist1->DrawCopy("ehist");


  hist2->DrawCopy("p same");
  
  
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.4); //36
  gPad->SetRightMargin(0.05);
  gPad->SetLeftMargin(0.15);

  TPad *rPad;
  rPad = new TPad("rPad","",0,0,1,1);
  rPad->SetFrameBorderMode(0);
  rPad->SetFillStyle(0);
  rPad->SetFillColor(0);
  rPad->SetBorderMode(0);
  rPad->SetBorderSize(0);
  rPad->Draw();
  rPad->cd();
  rPad->SetLogy(0);
  rPad->SetLogx(logx);
  gPad->SetLogx(0);
      gPad->SetLogy(0);
  ratio->SetTitle("");
  ratio->SetMaximum(1.13);
  ratio->SetMinimum(0.87);
  ratio->GetXaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetNdivisions(505);

  ratio->SetLineWidth(2);

  ratio->SetStats(kFALSE);
  ratio->GetXaxis()->SetTitle(x_Axis);
  ratio->GetXaxis()->SetTitleSize(0.09);
  ratio->GetXaxis()->SetTitleOffset(0.7);
  ratio->GetXaxis()->SetLabelSize(0.06);
  //ratio->GetYaxis()->CenterTitle();
  ratio->GetYaxis()->SetTitle("ratio");
  ratio->GetYaxis()->SetTitleSize(0.04);//06
  ratio->GetYaxis()->SetTitleOffset(1.12);
  ratio->GetYaxis()->SetLabelSize(0.04); // 0.06
  ratio->GetYaxis()->SetLabelOffset(0.01);
  ratio->GetYaxis()->SetTickLength(0.1);
  ratio->GetXaxis()->SetNoExponent();
  //ratio->GetYaxis()->SetNoExponent();

  ratio->DrawCopy("hist");
  drawunityline(ratio);
  rPad->SetTopMargin(0.6);
  rPad->SetBottomMargin(0.18);
  rPad->SetRightMargin(0.05);
  rPad->SetLeftMargin(0.15); //????
  rPad->SetTickx(1);
  rPad->SetTicky(1);


}

void format_tprofile(TProfile *histo,Int_t markerstyle = 20, Int_t markercolor = 2)
{
  histo->SetMarkerStyle(markerstyle);
  histo->SetMarkerSize(0.6);
  histo->SetMarkerColor(markercolor);
}

void format_data(TH1 *histo,Int_t markerstyle = 20, Int_t markercolor = 4)
{
  histo->SetMarkerStyle(markerstyle);
  histo->SetMarkerSize(0.6);//0.9
  histo->SetMarkerColor(markercolor);
}

void format_mc(TH1 *histo,Int_t linestyle = 1, Int_t linecolor = 3)
{
  histo->SetLineWidth((Short_t)2.0); //3
  histo->SetLineStyle(linestyle);
  histo->SetLineColor(linecolor);
  histo->SetFillStyle(0);
  histo->SetFillColor(17);
}

void norm_format_data(TH1 *histo, Int_t markerstyle = 20, Int_t markercolor = 2)
{
  // histo->Sumw2();
  histo->Scale(1.0/histo->Integral());
  histo->SetMarkerStyle(markerstyle);
  histo->SetMarkerSize(0.6);
  histo->SetMarkerColor(markercolor);
}

void norm_format_mc(TH1 *histo, Int_t linestyle = 1, Int_t linecolor = 3)
{
  //histo->Sumw2();
  histo->Scale(1.0/histo->Integral());
  //histo->SetLineWidth(3);//2.5
  histo->SetLineWidth(2.0);//2.5

  histo->SetLineStyle(linestyle);
  histo->SetLineColor(linecolor);
  histo->SetFillStyle(0);
  histo->SetFillColor(17);
}

Double_t correlationtest(TH2* histo)
{
  Int_t xfirst = histo->GetXaxis()->GetFirst();
  Int_t yfirst = histo->GetYaxis()->GetFirst();
  Int_t xlast = histo->GetXaxis()->GetLast();
  Int_t ylast = histo->GetYaxis()->GetLast();
  
  if(xlast >= 200 || ylast >= 200)
    return 0.0;
  
  Double_t u[200]; // u_i = summe ueber alle y werte
  Double_t v[200]; // u_i = summe ueber alle x werte
  
  for(Int_t i = 0; i < 200; i++)
    {
      u[i] = 0.0;
      v[i] = 0.0;
    }
  
  Double_t sumy = 0.0;
  Double_t sumx = 0.0;
  
  for(Int_t i = xfirst; i <= xlast; i++)
    {
      sumy = 0.0;
      for(Int_t t = yfirst; t <= ylast; t++)
	{
	  sumy += histo->GetBinContent(i,t);;
	}
      u[i-1] = sumy;
     }

  for(Int_t z = yfirst; z <= ylast; z++)
    {
      sumx = 0.0;
      for(Int_t k = xfirst; k <= xlast; k++)
	{
	  sumx += histo->GetBinContent(k,z);;
	}
      v[z-1] = sumx;
    }
  
  Double_t sumn = 0.0;
  for(Int_t l = 0; l < xlast; l++)
    sumn += u[l];
  
  Double_t chi = 0.0;
  Int_t empty_cells = 0;
   for(Int_t i = xfirst; i <= xlast; i++)
    {
      for(Int_t t = yfirst; t <= ylast; t++)
	{
	  if(u[i-1] != 0.0 && v[t-1] != 0.0)
	    {
	      Double_t n = u[i-1] * v[t-1] / sumn;
	      Double_t f = histo->GetBinContent(i,t);
	      chi += (f - n) * (f - n) / n;
	    }
	  else
	    empty_cells++;
	}
    }
   cout << "empty cells: " << empty_cells << endl;
   Int_t ndf = ( xlast - 1 ) * ( ylast - 1 );
   Double_t prob = TMath::Prob(chi, ndf);
   return prob;
}

/* TH1D fitslicesy_gaus(TH2F *h) */
/* { */
/*   h->FitSlicesY(); */
  
/*   TString title = h->GetName(); */

/*   TString ti1 = title + "_1"; */
/*   TString ti2 = title + "_2"; */

/*   TH1D *mean = (TH1D*)gDirectory->Get(ti1); */
/*   TH1D *sigma = (TH1D*)gDirectory->Get(ti2); */

/*   Int_t xlast = mean->GetXaxis()->GetLast(); */
/*   Int_t xfirst = mean->GetXaxis()->GetFirst(); */

/*   for(Int_t i = xfirst; i <= xlast; i++) */
/*     { */
/*       Double_t error = sigma->GetBinContent(i); */
/*       mean->SetBinError(i,error); */
/*     } */

/*   return *mean; */
/* } */


Double_t rounderror(Double_t value, Double_t error, Int_t k = 2)
{
  Double_t start = 1e20;
  Bool_t ok = kTRUE;
  Double_t ex = 0.0;
  
  for(Int_t i = 0; i < 40 && ok; i++)
    {
      ex = start*pow(10.0,-1*i);
      Int_t tmp = static_cast<Int_t>(error/ex);
      if(tmp != 0)
 	{
	  ok = kFALSE;
	}
    }
  Double_t blaexp = pow(10.0,-(k-1));
  Double_t tempo = static_cast<Int_t>(value/(ex*blaexp)+0.5)*ex*blaexp;
  return tempo;
}



/* class newtlegend : public TLegend */
/* { */
/*  public: */
/*   newtlegend(Double_t x1, Double_t y1, Double_t x2, Double_t y2) : TLegend(x1,y1,x2,y2) */
/*     { */
/*       SetBorderSize(1); */
/*       SetFillColor(19); */
/*       SetFillStyle(0); */
/*       cout << "constructor " << endl; */
/*     } */
/*   void AddData(TH1 *h, TString text, TString option) */
/*     { */
/*       data = h; */
/*       AddEntry(h,text,option); */
/*     } */
/*   void AddMC(TH1 *h, TString text, TString option) */
/*     { */
/*       cout << "1" << endl; */
/*       Double_t chi2 =0.0; */
/*       Int_t ndf = 0, igood = 0; */
/*       Double_t p = data->Chi2TestX(h,chi2, ndf, igood,"UU NORM P CHI2/NDF"); */
/*       TString tmp; */
/*       char tmpstring[20]; */
/*       sprintf(tmpstring,": #chi^{2} / ndf = %1.2f / %1.0i",chi2,ndf); */
/*       tmp = tmpstring; */
/*       cout << text << " " << tmp << endl; */
/*       TString t = ""; */
/*       t = text + tmp; */
/*       AddEntry(h,t,option); */
/*       cout << GetNRows() << endl; */
/*     } */
/*   virtual void Draw() */
/*     { */
/*       cout << "5" << endl; */
/*       Draw(); */
/*       // AppendPad(""); */
/*     } */
/*  private: */
/*   TH1 *data; */
/* }; */


#endif
