/******************************************************************************************
 * Original code by Joerg Behr & Igor Rubinskiy                                           *
 * 2012                                                                                   *
 *           - functions to plot 6 Mimosa26 planes residuals                              *
 *             and extract plane intrinsic resolution                                     *
 *             along with telescope pointing resolution                                   *
 *                                                                                        *
 *                                                                                        *
 ******************************************************************************************/

// rewrite as of 03.07.2013: Thomas Eichhorn


#include "configure_histos.h"
#include "AnaTel.h"
#include "TF1.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <list>
#include <ctime>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <utility>
#include <map>
#include <sstream>

//Root headers
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TDirectory.h"
#include "TMinuit.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTree.h"
#include "TObject.h"
#include "TH1.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TGaxis.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TStopwatch.h"
#include "TLatex.h"

TFile* _outputFile;


TH1D *delta0 = new TH1D("temp0","dXY 0",1000, -1, 5);
TH1D *delta1 = new TH1D("temp1","dXY 1",1000, -1, 5);
TH1D *delta2 = new TH1D("temp2","dXY 2",1000, -1, 5);
TH1D *delta3 = new TH1D("temp3","dXY 3",1000, -1, 5);
TH1D *delta4 = new TH1D("temp4","dXY 4",1000, -1, 5);
TH1D *delta5 = new TH1D("temp5","dXY 5",1000, -1, 5);
TH1D *delta6 = new TH1D("temp6","dX 0 and 5",1000, -1, 5);
TH1D *delta7 = new TH1D("temp7","dX 1 and 4",1000, -1, 5);
TH1D *delta8 = new TH1D("temp8","dX 2 and 3",1000, -1, 5);
TH1D *delta9 = new TH1D("temp9","dY 0 and 5",1000, -1, 5);
TH1D *delta10 = new TH1D("temp10","dY 1 and 4",1000, -1, 5);
TH1D *delta11 = new TH1D("temp11","dY 2 and 3",1000, -1, 5);



// Intrinsic resolution
Double_t m26_resolution =0.;
Double_t m26_res_error =0.;
Double_t global_plot_error = 0;


// Beam energy
Double_t global_beam = 0.;
Double_t global_spread = 0.;
Double_t global_thickness = 0.;

// Average noise and efficiency
Double_t avgnoise = 0.;
Double_t avgeffi = 0.;
Double_t avgnoise_error = 0.;
Double_t avgeffi_error = 0.;

// Average clustersize
Double_t avgclustersize =0;
Double_t avgclustersize_error =0;

Double_t avgmeas =0;
Double_t avgmeas_error =0;

// Plot options
Bool_t plot_residuals = false ;

// Verbosity
Bool_t verbose0 = false;
Bool_t verbose1 = false;
Bool_t verbose = true;



// Prediction graphs for the smilie plot
#define ngraphs 6

// Total plane count
#define nplanes 6

// Plane position for the smilie plot
Double_t posx[nplanes] = { 0.0, 150.0, 300.0, 450.0, 600.0, 750.0 };
Double_t posx_error[nplanes] = { 2.5,2.5,2.5,2.5,2.5,2.5};


//char telescopebuild[50];
string telescopebuild;
int planedistance;

// The observed resolution & error
Double_t obsresol_x[nplanes];
Double_t obsresol_y[nplanes];
Double_t obsresol_error_x[nplanes] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };
Double_t obsresol_error_y[nplanes] = { 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 };

// Put a cut on the threshold, 0 for all
Int_t global_thresh = 0;

// Cut in clustersize, set to 0 for all
#define clusterlimit 0

// What subsensor information do we want? Format: A-D for submatrix, 1-7 for clustersize, empty for all.
// e.g.: submask = "A1" or submask = "3" or submask = "C"
TString submask="";

// 2 Dimensional chi2 (X, Y) with uncertainty in the fit for beam, spread and sensor thickness
class  MyFunctionObject2D_beam {
  public:
    MyFunctionObject2D_beam(AnaTel *t, Double_t *arr1, Double_t *err1, Double_t *arr2, Double_t *err2)
    {
      tl = t;
      measured1 = arr1;
      error1    = err1;
      measured2 = arr2;
      error2    = err2;
    }
    double getchi2_2D(Double_t *par, Double_t *p)
    {
      // cout << "chi2call" << endl;
      tl->SetResolution(par[0]);
      //    tl->SetBeam(par[1], par[2]);
      tl->SetBeam(par[1], 0.0);
      tl->SetThickness(par[3]);
      Double_t chi2 =0.0;
      for(Int_t j = 0; j < nplanes; j++)
	//for(Int_t j = 1; j < nplanes-1; j++) // hj
      {
	float w = tl->GetWidth(j,0)*1000.0;
	chi2 += pow( (w - measured1[j])/(error1[j]) ,2) + pow( (w - measured2[j])/(error2[j]) ,2) ; 
	//if(j == 0) cout << "sigma_hat 0 = " <<  w ;
	//if(j == 3) cout << "  sigma_hat 3 = " <<  w ;
	if(j == nplanes -1) cout << " -> chi2/11 = " << chi2/11. << endl;

      }
      return chi2;

    }
    double getchi2_1D(Double_t *par, Double_t *p)
    {
      tl->SetResolution(par[0]);
      //tl->SetBeam(par[1], par[2]);
      tl->SetBeam(par[1], 0.0);
      tl->SetThickness(par[3]);
      Double_t chi2 =0.0;
      for(Int_t j = 0; j < nplanes; j++ )
      {
	chi2 +=  pow( (tl->GetWidth(j,0)*1000.0 - measured1[j])/(error1[j]) ,2) ; 
      }
      return chi2;
    }

    AnaTel *tl;
    Double_t *measured1;
    Double_t *error1   ;
    Double_t *measured2;
    Double_t *error2   ;
};

/////////////////////////////////////////////////////////////////////////////////////////////
MyFunctionObject2D_beam* fcn_beam=0;

// use MINUIT: 
void fcn_wrapper(int &npar, double *gin, double &f, double *par, int iflag)
{

  double p[4];
  p[0] = 0.0;
  p[1] = 0.0;
  p[2] = 0.0;
  p[3] = 0.0;
  f = fcn_beam->getchi2_2D(par,p);
}
//////////////////////////////////////////////////////////////////////////////////////////////

// 1 Dimensional chi2 (X or Y)
class  MyFunctionObject1D {
  public:
    MyFunctionObject1D(AnaTel *t, Double_t *arr1, Double_t *err1)
    {
      tl = t;
      measured1 = arr1;
      error1    = err1;
    }
    double getchi2_1D(double resolution)
    {
      tl->SetResolution(resolution);
      Double_t chi2 = 0.0;
      for(Int_t j = 0; j < nplanes; j++)
      {
	chi2 += pow( (tl->GetWidth(j,0)*1000.0 - measured1[j])/(error1[j]) ,2);
      }
      return chi2;
    }
    double operator() (double *x, double *p)
    {
      Double_t chi2 = 0.0;
      for(Int_t j = 0; j < nplanes; j++)
      {
	chi2 = getchi2_1D(x[0]) + p[0];
      }
      return chi2;
    }
    AnaTel *tl;
    Double_t *measured1;
    Double_t *error1;
};


// 2 Dimensional chi2 (X,Y)
class  MyFunctionObject2D {
  public:
    MyFunctionObject2D(AnaTel *t, Double_t *arr1, Double_t *arr2, Double_t *err1, Double_t *err2)
    {
      tl = t;
      measured1 = arr1;
      measured2 = arr2;
      error1 = err1;
      error2 = err2;
    }
    double getchi2_2D(double resolution)
    {



      tl->SetResolution(resolution);
      Double_t chi2 =0.0;
      for(Int_t j = 0; j < nplanes; j++)
      {
	chi2 += pow((tl->GetWidth(j,0)*1000.0 - measured1[j])/(error1[j]) ,2) + pow( (tl->GetWidth(j,0)*1000.0 - measured2[j])/(error2[j]) ,2);


      }
      return chi2;
    }
    double operator() (double *x, double *p)
    {
      Double_t chi2 = 0.0;
      for(Int_t j = 0; j < nplanes; j++)
      {
	chi2 = getchi2_2D(x[0]) + p[0];
      }
      return chi2;
    }
    AnaTel *tl;
    Double_t *measured1;
    Double_t *measured2;
    Double_t *error1;
    Double_t *error2;
};


// Get the telescope resolution from the measured input via k-factor
Double_t resol_estimate(const Double_t measured, const Int_t pos)
{
  std::vector<Double_t> vec_z;
  for(int i=0; i < nplanes; i++)
  {
    vec_z.push_back(posx[i]);
  }

  const Double_t offset = vec_z[pos];
  const Double_t n = (Double_t)vec_z.size();
  Double_t s2sum = 0.0;
  Double_t sum = 0.0;
  for(Int_t i =0; i < nplanes;i++)
  {
    if(i != pos)
    {
      s2sum += pow((vec_z[i] - offset), 2);
      sum += (vec_z[i] - offset);
    }
  }
  const Double_t k = s2sum / (n*s2sum-pow(sum,2));

  Double_t tel_res = k /(1+k) * measured * measured;
  Double_t single_point = TMath::Sqrt(measured * measured - k * tel_res);
  return TMath::Sqrt(tel_res);
}

// Get the telescope resolution from the measured input via k-factor
Double_t get_k(const Int_t pos)
{
  std::vector<Double_t> vec_z;
  for(int i=0; i < nplanes; i++)
  {
    vec_z.push_back(posx[i]);
  }

  const Double_t offset = vec_z[pos];
  const Double_t n = (Double_t)vec_z.size();
  Double_t s2sum = 0.0;
  Double_t sum = 0.0;
  for(Int_t i =0; i < nplanes;i++)
  {
    if(i != pos)
    {
      s2sum += pow((vec_z[i] - offset), 2);
      sum += (vec_z[i] - offset);
    }
  }
  const Double_t k = s2sum / ((n-1)*s2sum-pow(sum,2));

  return k;
}

// Run minuit
void run_global( Double_t ebeam, Double_t *obsresol_x, Double_t* obsresol_error_x, Double_t* obsresol_y, Double_t* obsresol_error_y )
{
  cout << "Run global" << endl;
  // Create telescope
  AnaTel *tlb = new AnaTel(telescopebuild.c_str());
  cout << "Tscope created" << endl;
  // Give original beam energy with spread // assume spread to be negligible
  tlb->SetBeam( ebeam, 0.0 );
  cout << " E = " << ebeam << endl;
  cout << " Pointing reso estimation at plane 0 using given initial reso (" << *(tlb->GetResolution()) << ") = " << tlb->GetError(0,0) << endl;
  cout << " Pointing reso estimation at plane 3 using given initial reso (" << *(tlb->GetResolution()) << ") = " << tlb->GetError(3,0) << endl;

  fcn_beam = new MyFunctionObject2D_beam(tlb, obsresol_x, obsresol_error_x, obsresol_y, obsresol_error_y );
  /*
     printf( "Plane 0  %8.5f %8.5f \n", tlb->GetWidth(0,0), tlb->GetWidth(0,1) );
     printf( "Plane 1  %8.5f %8.5f \n", tlb->GetWidth(1,0), tlb->GetWidth(1,1) );
     printf( "Plane 2  %8.5f %8.5f \n", tlb->GetWidth(2,0), tlb->GetWidth(2,1) );
     printf( "Plane 3  %8.5f %8.5f \n", tlb->GetWidth(3,0), tlb->GetWidth(3,1) );
     printf( "Plane 4  %8.5f %8.5f \n", tlb->GetWidth(4,0), tlb->GetWidth(4,1) );
     printf( "Plane 5  %8.5f %8.5f \n", tlb->GetWidth(5,0), tlb->GetWidth(5,1) );
   */

  /*
     bool firstminuitcall = true;
     if(firstminuitcall)
     {
     gSystem->Load("libMinuit");//is this really needed?
     firstminuitcall = false;
     }*/

  //initialize TMinuit with a maximum of 4 params
  TMinuit *gMinuit = new TMinuit(4);

  // set print level (-1 = quiet, 0 = normal, 1 = verbose)
  gMinuit->SetPrintLevel(1);

  // give the function
  gMinuit->SetFCN(fcn_wrapper);
  double arglist[10];
  int ierflg = 0;

  // minimization strategy (1 = standard, 2 = slower)
  arglist[0] = 2;
  gMinuit->mnexcm("SET STR",arglist,2,ierflg);

  // set error definition (1 = for chi square)
  arglist[0] = 1;
  gMinuit->mnexcm("SET ERR",arglist,1,ierflg);

  // Set starting values and step sizes for parameters
  Double_t vstart[4] = {0.0050, ebeam,   0.0, 0.05 };
  Double_t step[4]   = {0.0001,   0.1, 0.001, 0.001};

  /* orig:
     gMinuit->mnparm(0, "m26resol" , vstart[0], step[0],    0.0001,      0.020, ierflg);
     gMinuit->mnparm(1, "pbeam"    , vstart[1], step[1], 0.1*ebeam,    2*ebeam, ierflg);
     gMinuit->mnparm(2, "spread"   , vstart[2], step[2],      0.00,        0.1, ierflg);
     gMinuit->mnparm(3, "thickness", vstart[3], step[3],      0.05,      0.100, ierflg);

   */
  gMinuit->mnparm(0, "m26resol" , vstart[0], step[0],    0.001,      0.020, ierflg);
  gMinuit->mnparm(1, "pbeam"    , vstart[1], step[1], 0.9*ebeam,  1.1*ebeam, ierflg);
  gMinuit->mnparm(2, "spread"   , vstart[2], step[2],     0.00,        0.1, ierflg);
  gMinuit->mnparm(3, "thickness", vstart[3], step[3],     0.04,        0.06, ierflg);

  //gMinuit->FixParameter(0);
  gMinuit->FixParameter(1);
  gMinuit->FixParameter(2);
  gMinuit->FixParameter(3);

  // Now ready for minimization step
  arglist[0] = 100 ;
  arglist[1] = 0.001 ;
  gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
  //gMinuit->mnexcm("IMPROVE", arglist ,1,ierflg);
  //gMinuit->mnexcm("IMPROVE", arglist ,1,ierflg);

  //gMinuit->Release(2);
  //gMinuit->Release(3);

  ////  Now ready for minimization step
  //arglist[0] = 8000;
  //arglist[1] = 1.0;
  //gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);


  //           calculate errors using MINOS. do we need this?
  //           arglist[0] = 2000;
  //           arglist[1] = 0.1;
  //           gMinuit->mnexcm("MINOS",arglist,1,ierflg);

  //   get results from migrad
  double par[4];
  double perr[4];
  gMinuit->GetParameter(0, par[0], perr[0]);
  gMinuit->GetParameter(1, par[1], perr[1]);
  gMinuit->GetParameter(2, par[2], perr[2]);
  gMinuit->GetParameter(3, par[3], perr[3]);



  /*          printf("%s 0 %8.4f %8.4f\n", par[0], perr[0] );
	      printf("%s 1 %8.4f %8.4f\n", par[1], perr[1] );
	      printf("%s 2 %8.4f %8.4f\n", par[2], perr[2] );
	      printf("%s 3 %8.4f %8.4f\n", par[3], perr[3] );

	      printf("%s chi2 %8.4f \n", fcn_beam->getchi2_2D(par, perr) );
	      printf("%s tel_resolution at 0 %8.4f \n", resol_estimate(par[0],0));
	      printf("%s tel_resolution at 1 %8.4f \n", resol_estimate(par[0],1));
	      printf("%s tel_resolution at 2 %8.4f \n", resol_estimate(par[0],2));
	      printf("%s tel_resolution at 3 %8.4f \n", resol_estimate(par[0],3));
	      printf("%s tel_resolution at 4 %8.4f \n", resol_estimate(par[0],4));
	      printf("%s tel_resolution at 5 %8.4f \n", resol_estimate(par[0],5));
   */

  cout << "Resolution was " << m26_resolution*1000.0 << " mu m" << endl;
  m26_resolution = par[0];
  cout << "Resolution is " << m26_resolution*1000.0 << " mu m" << endl;
  cout << "Resolution error was " << m26_res_error*1000.0 << " mu m" << endl;
  m26_res_error  = perr[0];
  cout << "Resolution error is " << m26_res_error*1000.0 << " mu m" << endl;
  cout << "Beam was " << global_beam << endl;
  global_beam      = par[1];
  cout << "global beam now " << global_beam << endl;
  cout << "Spread was " << global_spread << endl;
  global_spread    = par[2];
  cout << "Spread is " << global_spread << endl;
  cout << "Thickness was " << global_thickness << endl;
  global_thickness = par[3];
  cout << "Thickness is " << global_thickness << endl;

  Double_t tempdist = posx[1] - posx[0];
  Double_t scatterer= 0.0136/global_beam * sqrt(global_thickness/93.66 + tempdist/304200.0 + 0.1/286.0) * (1.+0.038*std::log(global_thickness/93.66 + tempdist/304200.0 + 0.05/286.0)) ;

  cout << "Scattering is " << scatterer << endl;
}


// Fitting of each file -> resolution
void fitter(Int_t runnumber, Double_t ebeam)
{

  int nominalbeam = ebeam;

  //  TString filename("histograms/run00");
  TString filename("../filteredhistos/run00");
  if (runnumber <= 999)
    filename+="0";
  if (runnumber <= 99)
    filename+="0";
  filename+=runnumber;
  filename+="-fitter.root";
  TFile *f6 = new TFile(filename);
  TH1D *h_m26[12];
  f6->cd();

  // If the clustersize is limited, put clustersize 1 into the original histogram
  if(clusterlimit >0)
    submask="1";

  // Load file
  h_m26[0] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftX"+submask))->Clone();
  h_m26[1] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftY"+submask))->Clone();
  h_m26[2] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftX"+submask))->Clone();
  h_m26[3] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftY"+submask))->Clone();
  h_m26[4] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftX"+submask))->Clone();
  h_m26[5] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftY"+submask))->Clone();
  h_m26[6] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftX"+submask))->Clone();
  h_m26[7] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftY"+submask))->Clone();
  h_m26[8] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftX"+submask))->Clone();
  h_m26[9] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftY"+submask))->Clone();
  h_m26[10] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftX"+submask))->Clone();
  h_m26[11] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftY"+submask))->Clone();

  // Add the clustersizes 2-4 into the histos
  if(clusterlimit >0)
  {
    cout << "Clustersize is limited. This is only implemented up to 4!" << endl;
    TH1D *h_m26_2[12];
    TH1D *h_m26_3[12];
    TH1D *h_m26_4[12];
    h_m26_2[0] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftX2"))->Clone();
    h_m26_2[1] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftY2"))->Clone();
    h_m26_2[2] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftX2"))->Clone();
    h_m26_2[3] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftY2"))->Clone();
    h_m26_2[4] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftX2"))->Clone();
    h_m26_2[5] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftY2"))->Clone();
    h_m26_2[6] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftX2"))->Clone();
    h_m26_2[7] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftY2"))->Clone();
    h_m26_2[8] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftX2"))->Clone();
    h_m26_2[9] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftY2"))->Clone();
    h_m26_2[10] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftX2"))->Clone();
    h_m26_2[11] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftY2"))->Clone();

    h_m26_3[0] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftX3"))->Clone();
    h_m26_3[1] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftY3"))->Clone();
    h_m26_3[2] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftX3"))->Clone();
    h_m26_3[3] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftY3"))->Clone();
    h_m26_3[4] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftX3"))->Clone();
    h_m26_3[5] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftY3"))->Clone();
    h_m26_3[6] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftX3"))->Clone();
    h_m26_3[7] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftY3"))->Clone();
    h_m26_3[8] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftX3"))->Clone();
    h_m26_3[9] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftY3"))->Clone();
    h_m26_3[10] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftX3"))->Clone();
    h_m26_3[11] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftY3"))->Clone();

    h_m26_4[0] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftX4"))->Clone();
    h_m26_4[1] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftY4"))->Clone();
    h_m26_4[2] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftX4"))->Clone();
    h_m26_4[3] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftY4"))->Clone();
    h_m26_4[4] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftX4"))->Clone();
    h_m26_4[5] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftY4"))->Clone();
    h_m26_4[6] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftX4"))->Clone();
    h_m26_4[7] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftY4"))->Clone();
    h_m26_4[8] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftX4"))->Clone();
    h_m26_4[9] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftY4"))->Clone();
    h_m26_4[10] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftX4"))->Clone();
    h_m26_4[11] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftY4"))->Clone();

    // Add the clustersize 2-4 histos into the "original" one and continue
    for(int i=0;i<12;i++)
    {
      h_m26_3[i]->Add(h_m26_4[i]);
      h_m26_2[i]->Add(h_m26_3[i]);
      h_m26[i]->Add(h_m26_2[i]);
    }
  }

  // Canvas and range
  for(int i=0;i<12;i++)
    h_m26[i]->GetXaxis()->SetRangeUser(-0.05,0.05);
  TCanvas *canv;
  canv = new TCanvas("m26fitter","m26fitter",900,10,600,800); 
  gStyle->SetPadBorderMode(0);
  gStyle->SetOptStat(0);
  canv->SetFillColor(0);
  canv->Divide(2,6);

  // TCanvas *canv2;
  // canv2 = new TCanvas("m26fitter2","m26fitter2",900,10,600,800); 
  // canv2->SetFillColor(0);
  // canv2->Divide(2,3);

  if (verbose1)
  {
    cout << " " << endl;
    cout << "Reading residuals from file..." << endl;
    cout << " " << endl;
  }

  // Sort into X or Y
  for(Int_t i = 0; i < nplanes*2; i++ )
  {
    canv->cd(i+1);
    gStyle->SetErrorX(0);
    gPad->SetLogx(0);
    gPad->SetLogy(0);

    int j = 6;
    if(i == 0 || i == 1)
      j = 0;
    if(i == 2 || i == 3)
      j = 1;
    if(i == 4 || i == 5)
      j = 2;
    if(i == 6 || i == 7)
      j = 3;
    if(i == 8 || i == 9)
      j = 4;
    if(i == 10 || i == 11)
      j = 5;

    TString st = "";
    char tmpstring[500];
    sprintf(tmpstring, "Sensor %1.1i", j);
    st = tmpstring;

    if(plot_residuals)
    {
      TString xtitle;
      if(i == 0 || i == 2 || i == 4 ||  i == 6 || i == 8 || i == 10)
	xtitle = "(x_{pred.} - x_{meas.}) / mm";
      else
	xtitle = "(y_{pred.} - y_{meas.}) / mm";

      format_mc(h_m26[i], 1, kBlue);
      //histo_cfg(h_m26[i],xtitle , "tracks", st);
      h_m26[i]->GetYaxis()->SetNoExponent();
      h_m26[i]->GetYaxis()->SetTitleOffset(1.68);
      h_m26[i]->GetXaxis()->SetTitleOffset(0.97);
      h_m26[i]->GetYaxis()->SetTitleSize(0.035);
      h_m26[i]->GetXaxis()->SetTitleSize(0.035);
      h_m26[i]->GetYaxis()->SetLabelSize(0.035);
      h_m26[i]->GetXaxis()->SetLabelSize(0.035);
      h_m26[i]->GetXaxis()->SetLabelOffset(0.004);
      gPad->SetLeftMargin(0.2238535);
      h_m26[i]->SetMaximum(h_m26[i]->GetMaximum()*1.3);
      h_m26[i]->DrawCopy("hist");
    }

    // Fit the residuals
    TF1 *f1 = new TF1("f1","gaus");
    TF1 *f2 = new TF1("f2","gaus");
    f2->SetLineWidth(2);
    f2->SetLineStyle(1);
    f2->SetLineColor(kBlack);
    h_m26[i]->Fit(f1,"EMQI0","");
    Double_t mean_1 = f1->GetParameter(1);
    Double_t sigma_1 = f1->GetParameter(2);
    Double_t sigma_error_1 = f1->GetParError(2);

    // Repeat within 2 sigma

    float nsig = 2.5;
    h_m26[i]->Fit(f2,"EQMI","", (mean_1 - nsig*sigma_1), (mean_1 + nsig*sigma_1));
    Double_t mean_2 = f2->GetParameter(1);

    if (verbose0)
      cout << " Measured residual mean plane " << j << " is " << (mean_2 *1000.0) << " mu m" << endl;

    Double_t sigma_2 = f2->GetParameter(2) * 1000.0;
    cout << " measured width in um = " << sigma_2 << endl;
    //Double_t sigma_error_2 = f2->GetParError(2) * 1000.0;
    Double_t sigma_error_2 = 0.05 * sigma_2;

    // Now get resolution estimate
    Double_t tel_resol = resol_estimate(sigma_2,j);

    if (verbose0)
      cout << " Measured residual width plane " << j << " is " <<  tel_resol << " mu m" << endl;

    // Fill fitted resolution into array
    if(i == 0 || i == 2 || i == 4 ||  i == 6 || i == 8 || i == 10)
    {
      obsresol_x[j] = sigma_2;
      obsresol_error_x[j] = sigma_error_2;
      //   obsresol_x[j] = sigma_1*1000.0;
      //  obsresol_error_x[j] = sigma_error_1*1000.0;
    }
    else
    {
      obsresol_y[j] = sigma_2; 
      obsresol_error_y[j] = sigma_error_2;
      //  obsresol_y[j] = sigma_1*1000.0; 
      //  obsresol_error_y[j] = sigma_error_1*1000.0;
    }

    // Plot this
    if(plot_residuals)
    {
      TPaveText *label;
      label = new TPaveText(0.2588087,0.7992021,0.6363255,0.8823138, "brNDC");
      label->SetTextAlign(11);
      label->SetTextFont(22);
      label->SetTextSize(0.08);
      label->SetFillStyle(0);
      label->SetBorderSize(0);
      TString sigmatext = "";
      char tmpstring2[500];
      sprintf(tmpstring2, "#sigma = (%1.2f #pm %1.2f) #mum", sigma_2, sigma_error_2);
      sigmatext = tmpstring2;
      label->AddText(sigmatext);
      label->Draw();
    }
  }

  TString outputname="run_";
  outputname+=runnumber;
  outputname+="_";
  outputname+=submask;
  //canv->Print("pics/"+outputname+"_measured_residuals.eps");
  _outputFile->cd();
  canv->Write();
  canv->Close();

  // Run minimization with energy and measured resolutions
  run_global( ebeam, obsresol_x, obsresol_error_x, obsresol_y, obsresol_error_y );
  ebeam = global_beam;
  if(verbose0)
    cout << "Beam energy after first minimization now at " << global_beam << " GeV" << endl;



  /*
   * 
   * 
   * 

  // Create 4 new "telescopes"
  AnaTel *tl = new AnaTel(telescopebuild);
  tl->SetBeam( ebeam );
  AnaTel *tl2 = new AnaTel(telescopebuild);
  tl2->SetBeam( ebeam );
  AnaTel *tl2x = new AnaTel(telescopebuild);
  tl2x->SetBeam( ebeam );
  AnaTel *tl2y = new AnaTel(telescopebuild);
  tl2y->SetBeam( ebeam );
  tl->SetThickness( global_thickness );
  tl2->SetThickness( global_thickness );
  tl2x->SetThickness( global_thickness );
  tl2y->SetThickness( global_thickness );
  tl->SetBeam( global_beam, global_spread);   
  tl2->SetBeam( global_beam, global_spread);  
  tl2x->SetBeam( global_beam, global_spread);
  tl2y->SetBeam( global_beam, global_spread);


  AnaTel *tsmile = new AnaTel(telescopebuild);
  tsmile->SetBeam( ebeam );
  tsmile->SetThickness( global_thickness );
  tsmile->SetBeam( global_beam, global_spread);
  tsmile->SetResolution(m26_resolution);

  // predicted resolution, fill pre with posx and prediction_y
  TGraph*  pre[ngraphs];
  Double_t prediction_y[nplanes];
  // = {
  //tl->GetWidth(0,0)*1000.0,
  //tl->GetWidth(1,0)*1000.0,
  //tl->GetWidth(2,0)*1000.0,
  //tl->GetWidth(3,0)*1000.0,
  //tl->GetWidth(4,0)*1000.0
  //};

  // if(verbose)
  //  cout << "Resolution estimate: " << tl->GetWidth(0,0) << endl;

  // The starting Mimosa26 singlepoint resolution
  Double_t singlepoint_resolution = 0.001;

  // Start predicting
  if (verbose1)
  {
  cout << "Calculating Resolution predictions" << endl;
  cout << " " << endl;
  }

  for(Int_t icurve = 0; icurve <  ngraphs; icurve++)
  {
  tl->SetResolution(singlepoint_resolution);
  if(verbose0)
  cout << " Resolution prediction on plane : " ;
  for(Int_t j = 0; j < nplanes; j++)
  {
  prediction_y[j] = tl->GetWidth(j,0)*1000.0;
  if(verbose0)
  printf(" %d: %8.3f", j, prediction_y[j] );
  }
  if(verbose0)
  cout << endl; 
  pre[icurve] = new TGraph(nplanes, posx, prediction_y);
  pre[icurve]->SetLineStyle(icurve+1);
  pre[icurve]->SetLineWidth(1.3);
  singlepoint_resolution += 0.0010; 
}
/*
// hacks
AnaTel *teltest = new AnaTel(telescopebuild);
MyFunctionObject2D *fobjtest = new MyFunctionObject2D(teltest, obsresol_x, obsresol_y, obsresol_error_x, obsresol_error_y );

teltest->SetThickness( global_thickness );
teltest->SetBeam( global_beam, global_spread);  
teltest->SetResolution( m26_resolution);

Double_t tempchi2  = fobjtest->getchi2_2D(m26_resolution );

cout << "the chi2 is " << tempchi2 << endl;




cout << " " << endl;

// cout << " observed " << obsresol_x[0] << endl;
// cout << " observed error " << obsresol_error_x[0] << endl;
cout << " 26 fit " << m26_resolution << endl;
cout << " 26 er " << m26_res_error << endl;


double the_resi = m26_resolution*1000.0;
double the_error = 0.1;

//m26_res_error*1000.0;


double temptemp = 0;


for (int i = 0; i< 6;i++)
{
//    temptemp+= (sqrt(obsresol_x[i]*obsresol_x[i]/(1+get_k(i))) - m26_resolution*1000.0)*(sqrt(obsresol_x[i]*obsresol_x[i]/(1+get_k(i))) - m26_resolution*1000.0)/ m26_res_error/1000.0/m26_res_error/1000.0 + (sqrt(obsresol_y[i]*obsresol_y[i]/(1+get_k(i))) - m26_resolution*1000.0)*(sqrt(obsresol_y[i]*obsresol_y[i]/(1+get_k(i))) - m26_resolution*1000.0)/ m26_res_error/1000.0/m26_res_error/1000.0;

cout << " observed " << obsresol_x[i] << endl;
cout << " observed error " << obsresol_error_x[i] << endl;

cout << " k is " << get_k(i) << endl;
cout << " sqrt ob^2/1+k is " << sqrt(obsresol_x[i]*obsresol_x[i]/(1+get_k(i))) << endl;

// temptemp += pow(((sqrt((obsresol_x[i]*obsresol_x[i]-obsresol_error_x[i]*obsresol_error_x[i])/(1+get_k(i))) - the_resi)/the_error) , 2);
// temptemp += pow(((sqrt((obsresol_y[i]*obsresol_y[i]-obsresol_error_y[i]*obsresol_error_y[i])/(1+get_k(i))) - the_resi)/the_error) , 2);

temptemp += pow(((sqrt((obsresol_x[i]*obsresol_x[i]+obsresol_error_x[i]*obsresol_error_x[i])/(1+get_k(i))) - the_resi)/the_error) , 2);
temptemp += pow(((sqrt((obsresol_y[i]*obsresol_y[i]+obsresol_error_y[i]*obsresol_error_y[i])/(1+get_k(i))) - the_resi)/the_error) , 2);

}

cout << endl;
cout << " CHI2/ndf " << temptemp/11.0 << endl;
cout << endl;





// Plot the chi2
TCanvas *canv_chi2 = new TCanvas("chi2","chi2",1400,10,500,500);
canv_chi2->Divide(1,1,0.02,0.02);
//  cout << "New resolution estimate: " << tl->GetWidth(0,0) << endl;
Double_t chi2_axis_Xmin = 0.0001;
Double_t chi2_axis_Xmax = 0.007;
Double_t epsilon = 0.1;
Int_t    niter = 10;

MyFunctionObject2D *fobj = new MyFunctionObject2D(tl2, obsresol_x, obsresol_y, obsresol_error_x, obsresol_error_y );
MyFunctionObject1D *fobjX = new MyFunctionObject1D(tl2x, obsresol_x, obsresol_error_x );
MyFunctionObject1D *fobjY = new MyFunctionObject1D(tl2y, obsresol_y, obsresol_error_y );

TF1* fchi2      = new TF1("fchi2",      fobj,  chi2_axis_Xmin, chi2_axis_Xmax, 1, "MyFunctionObject2D");
TF1* fchi2X     = new TF1("fchi2X",     fobjX, chi2_axis_Xmin, chi2_axis_Xmax, 1, "MyFunctionObject1D");
TF1* fchi2Y     = new TF1("fchi2Y",     fobjY, chi2_axis_Xmin, chi2_axis_Xmax, 1, "MyFunctionObject1D");


fchi2->SetLineColor(45);
fchi2->SetNpx(100);
fchi2->Draw("L");

fchi2X->SetLineColor(25);
fchi2X->SetNpx(100);
fchi2X->Draw("LSAME");

fchi2Y->SetLineColor(65);
fchi2Y->SetNpx(100);
fchi2Y->Draw("LSAME");

Double_t resolmin2D      = fchi2->GetMinimumX(chi2_axis_Xmin, chi2_axis_Xmax, epsilon, niter);
Double_t resolmin1DX     = fchi2X->GetMinimumX(chi2_axis_Xmin, chi2_axis_Xmax, epsilon, niter);
Double_t resolmin1DY     = fchi2Y->GetMinimumX(chi2_axis_Xmin, chi2_axis_Xmax, epsilon, niter);

Double_t chimin2D  = fobj->getchi2_2D(resolmin2D );
Double_t chimin1D_X = fobjX->getchi2_1D(resolmin1DX);
Double_t chimin1D_Y = fobjY->getchi2_1D(resolmin1DY);

Double_t  ymin = 0.;
Double_t  ymax = chimin2D*2.;

fchi2->SetMinimum(ymin);
fchi2X->SetMinimum(ymin);
fchi2Y->SetMinimum(ymin);

fchi2->SetMaximum(ymax);
fchi2X->SetMaximum(ymax);
fchi2Y->SetMaximum(ymax);

if(verbose)
  cout << "Minimal 2D resolution: " << resolmin2D << " . Minimal 2D Chi2: " << chimin2D << endl;

  cout << "hackchi2 " << endl;


  // starting values only:
  Double_t resol_p =  resolmin2D;
  Double_t resol_n = -resolmin2D/3.;
  Double_t resolX_p =  resolmin1DX;
  Double_t resolX_n = -resolmin1DX/3.;
  Double_t resolY_p =  resolmin1DY;
  Double_t resolY_n = -resolmin1DY/3.;

  Double_t value_x_p = fchi2->GetX( chimin2D+1., resolmin2D, resolmin2D+5., epsilon, niter  );
  //   Double_t value_x_n = fobj->GetX( chimin2D+1.);
  //cout << value_x_p <<  endl;

  Double_t fraction_step = TMath::Abs( value_x_p-resolmin2D ) /10.;

  Int_t nstep     = 100; // step to find the plus and minus range of the resolmin error bars
  Double_t dstep  = 0.01  ; // delta = 1/nstep and dstep/nstep sets the "resolution" of the error bars scan.
if(fraction_step < dstep/nstep )
  fraction_step = dstep/nstep;

  for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  double valuep= fobj->getchi2_2D(resolmin2D + delta );
  if( valuep>= chimin2D+1.0)
  {
    resol_p = + delta ;
    if(verbose)
      printf("p %d %8.3f %8.3f \n", i, resol_p,  valuep);
    break;
  }
}

for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  if(resolmin2D-delta<1e-7)
    break;
  double valuen= fobj->getchi2_2D(resolmin2D - delta );
  if( valuen>= chimin2D+1.0)
  {
    resol_n = -delta;
    if(verbose)
      printf("n %d %8.3f %8.3f \n", i, resol_n,  valuen);
    break;
  }
}

for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  double valuep= fobjX->getchi2_1D(resolmin1DX + delta );
  if( valuep>= chimin1D_X+1.0)
  {
    resolX_p = + delta ;
    if(verbose)
      printf("p %d %8.3f %8.3f \n", i, resolX_p,  valuep);
    break;
  }
}

for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  if(resolmin1DX-delta<1e-7)
    break;
  double valuen= fobjX->getchi2_1D(resolmin1DX-delta );
  if( valuen>= chimin1D_X+1.0)
  {
    resolX_n = -delta;
    if(verbose)
      printf("n %d %8.3f %8.3f \n", i, resolX_n,  valuen);
    break;
  }
}

for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  double valuep= fobjY->getchi2_1D(resolmin1DY + delta );
  if( valuep>= chimin1D_Y+1.0)
  {
    resolY_p = + delta ;
    if(verbose)
      printf("p %d %8.3f %8.3f \n", i, resolY_p,  valuep);
    break;
  }
}

for(int i=0;i<nstep;i++)
{
  double delta = i*fraction_step;
  if(resolmin1DY-delta<1e-7) break;
  double valuen= fobjY->getchi2_1D(resolmin1DY-delta );
  if( valuen>= chimin1D_Y+1.0)
  {
    resolY_n = -delta;
    if(verbose)
      printf("n %d %8.3f %8.3f \n", i, resolY_n,  valuen);
    break;
  }
}


TGraphErrors *x_direction = new TGraphErrors(nplanes, posx, obsresol_x, posx_error, obsresol_error_x);
TGraphErrors *y_direction = new TGraphErrors(nplanes, posx, obsresol_y, posx_error, obsresol_error_y);

resolmin2D=m26_resolution;

tl->SetResolution(resolmin2D);

// Get errors for t1
Double_t prediction_error[6];

// for the plot
float average_error = 0.0;


for(Int_t j = 0; j < nplanes; j++)
{
  tl->SetResolution(resolmin2D);
  prediction_y[j] = tl->GetWidth(j,0)*1000.0;
  tl->SetResolution(resolmin2D-0.0001);
  const Double_t up = TMath::Abs(tl->GetWidth(j,0)*1000.0 - prediction_y[j]);
  tl->SetResolution(resolmin2D+0.0001);
  const Double_t down = TMath::Abs(tl->GetWidth(j,0)*1000.0 - prediction_y[j]);
  prediction_error[j] = 0.5*(up+down);
  average_error += prediction_error[j];

}

for (int j = 0; j<6;j++)
cout << "asdfasdfasdfasdf " << get_k(j) << endl;


// fail!
// Errorband for the smilie plot
cout << "asdfasdf " << tsmile->GetWidth(0,0)*1000.0 << endl;
cout << "asdfasdf " << tsmile->GetWidth(1,0)*1000.0 << endl;
cout << "asdfasdf " << tsmile->GetWidth(2,0)*1000.0 << endl;
cout << "asdfasdf " << tsmile->GetWidth(3,0)*1000.0 << endl;


// fix error
for (int j=0;j<6;j++)
{
  //prediction_error[j] = prediction_error[j] / (sqrt(12.0));



  cout << endl;
  cout << " !!! " << endl;
  cout << endl;
  cout << "Prediction error on plane " << j << " is " << prediction_error[j] << endl;
  cout << endl;

  float tempup = 0.0;
  float tempdn = 0.0;

  tsmile->SetResolution(m26_resolution);

  float deltam26 = 0.1;


  double localscatter = sqrt(tsmile->GetWidth(j,0)*1000.0*tsmile->GetWidth(j,0)*1000.0 - m26_resolution*1000.0*m26_resolution*1000.0*(1 + get_k(j)));
  float deltams = 0.1*localscatter;
  cout << "scatter here " << localscatter << endl;

  float thedelta = sqrt(pow(((1+get_k(j))*m26_resolution*1000.0*deltam26/(tsmile->GetWidth(j,0)*1000.0)),2) + pow((localscatter*deltams/(tsmile->GetWidth(j,0)*1000.0)),2));

  cout << " THE DELTA " << thedelta << endl;
  cout << "error times k upper " << sqrt((m26_resolution*1000.0+prediction_error[j])*(m26_resolution*1000.0+prediction_error[j])*(1 + get_k(j)) + localscatter*localscatter) << endl;
  cout << "error times k lower " << sqrt((m26_resolution*1000.0-prediction_error[j])*(m26_resolution*1000.0-prediction_error[j])*(1 + get_k(j)) + localscatter*localscatter) << endl;



  tsmile->SetResolution(m26_resolution+prediction_error[j]/1000.0);
  tempup = tsmile->GetWidth(j,0)*1000.0;
  cout << " max " << m26_resolution*1000.0+prediction_error[j] << "    " << tsmile->GetWidth(j,0)*1000.0 << endl;
  tsmile->SetResolution(m26_resolution-prediction_error[j]/1000.0);
  cout << " max " << m26_resolution*1000.0-prediction_error[j] << "    " << tsmile->GetWidth(j,0)*1000.0 << endl;
  tempdn = tsmile->GetWidth(j,0)*1000.0;


  prediction_error[j] = (tempup - tempdn)/2.0;
  cout << "Prediction error on plane " << j << " is " << prediction_error[j] <<  " in sig meas! " << endl;
  cout << endl;


  prediction_error[j] = thedelta;

}




average_error = average_error / nplanes;

global_plot_error = average_error;


cout << "global_plot_error " << global_plot_error << endl;

global_plot_error = global_plot_error/ (sqrt(12.0));

cout << "global_plot_error / sqrt 12 " << global_plot_error << endl;


// Get fiterrors on t1 for plane 2
tl->SetResolution(resolmin2D);
Double_t fiterror_central =  tl->GetError(2,0);
tl->SetResolution(resolmin2D+0.0001);
Double_t fiterror_up = TMath::Abs(tl->GetError(2,0) - fiterror_central);
if(verbose)
  cout << "fiterror_up " << fiterror_up  << endl;
  tl->SetResolution(resolmin2D-0.0001);
  Double_t fiterror_down = TMath::Abs(tl->GetError(2,0) - fiterror_central);
if(verbose)
  cout << "fiterror_down " << fiterror_down  << endl;
  Double_t fiterror_error = 0.5*(fiterror_up + fiterror_down);

  canv_chi2->Print("pics/"+outputname+"_chi2.eps");
  _outputFile->cd();
  canv_chi2->Write();
  canv_chi2->Close();


  //TGraphErrors *sys_prediction = new TGraphErrors(nplanes, posx, prediction_y, posx_error, prediction_error);
  TGraphErrors *sys_prediction = new TGraphErrors(nplanes, posx, prediction_y);

  sys_prediction->SetLineWidth(3);

  TGraphErrors *sys_prediction2 = new TGraphErrors(nplanes, posx, prediction_y, posx_error, prediction_error);
  sys_prediction2->SetFillColor(kOrange);

  //


  TH1D *h_axis;

if (planedistance == 150)
{
  h_axis = new TH1D("h_axis","Intrinsic Resolution Calculation - #Delta_{z} = 150 mm",1, -75, 825.0);
  //h_axis = new TH1D("h_axis","Intrinsic Resolution Calculation - #Delta_{z} = 150 mm",100, -75, 825.0*2);
  cout << "planedistance 150" << endl;
} else {
  h_axis = new TH1D("h_axis","Intrinsic Resolution Calculation - #Delta_{z} = 20 mm",1, -10, 110.0);
  //h_axis = new TH1D("h_axis","Intrinsic Resolution Calculation - #Delta_{z} = 20 mm",100, -10, 110.0*2);
  cout << "planedistance 20" << endl;
}
TCanvas *smilie = new TCanvas("smilie","smilie",10,10,800,600);
smilie->SetRightMargin(0.3);
// Do smilie plot, x/y_direction has the measured residuals, pre the predictions and sys_prediction the errorband around the fit
//  gStyle->SetPadBorderMode(0);
//  gStyle->SetOptStat(0);
smilie->SetFillColor(0);
//  smilie->Divide(1,1);
//  smilie->cd(1);
gStyle->SetErrorX(0);

gPad->SetLogx(0);
gPad->SetLogy(0);
x_direction->SetMarkerStyle(21);
x_direction->SetMarkerColor(kRed);
x_direction->SetMarkerSize(2.5);
y_direction->SetMarkerStyle(22);
y_direction->SetMarkerColor(kBlue);
y_direction->SetMarkerSize(2.5);
h_axis->GetXaxis()->SetLabelFont(42);
h_axis->GetXaxis()->SetLabelSize(0.035);
h_axis->GetXaxis()->SetTitleSize(0.035);
h_axis->GetXaxis()->SetTitleFont(42);
h_axis->GetXaxis()->SetTitle("z in mm");
h_axis->GetYaxis()->SetLabelFont(42);
h_axis->GetYaxis()->SetLabelSize(0.035);
h_axis->GetYaxis()->SetTitleSize(0.035);
h_axis->GetYaxis()->SetTitleFont(42);
h_axis->GetYaxis()->SetTitle("#sigma_{meas} in #mum");
//histo_cfg(h_axis, "z in mm","#sigma_{meas.} in #mum","");
h_axis->SetMinimum(0.0);

// The y-axis can be adjusted to geometry too
if(telescopebuild == "AnaTel.geom")
h_axis->SetMaximum(30.);
if(telescopebuild == "AnaTel_thin.geom")

if (planedistance == 20)
  h_axis->SetMaximum(15.);
if (planedistance == 150)
  h_axis->SetMaximum(30.);

  h_axis->Draw("hist");

  sys_prediction2->Draw("LE3 same");
  sys_prediction->Draw("L same");
  x_direction->Draw("p same");
  y_direction->Draw("p same");

  // Draw predictions
  //for(Int_t i = 0; i < ngraphs; i++)
  for(Int_t i = 1; i < 5; i++)
  pre[i]->Draw("L same");
  /*
  // Legend for x and y
  TLegend *leg = new TLegend(0.19,0.70,0.50,0.80);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->AddEntry(x_direction,"x direction","p");
  leg->AddEntry(y_direction,"y direction","p");
  leg->SetTextSize(0.03);
  leg->Draw();

  // Create second legend with the predictions
  TLegend *leg2 = new TLegend(0.7,0.6,0.95,0.9);
  leg2->SetBorderSize(0);
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  //leg2->SetHeader("Prediction for intrinsic resolution:");
  char clabel[100]="";
  sprintf(clabel,  "#sigma_{M26} = (%4.2f #pm %4.2f) #mum", resolmin2D*1000., global_plot_error);
  //sprintf(clabel,  "#sigma_{DUT} = (%4.2f #pm %4.2f) #mum", resolmin2D*1000., m26_res_error*1000.);
  leg2->AddEntry(sys_prediction, clabel ,"l");
  //leg2->AddEntry(pre[0],"#sigma_{M26} = 1.0 #mum","l");
  leg2->AddEntry("","Prediction for intrinsic resolution:","h");
  leg2->AddEntry(pre[1],"#sigma_{M26} = 2.0 #mum","l");
  leg2->AddEntry(pre[2],"#sigma_{M26} = 3.0 #mum","l");
  leg2->AddEntry(pre[3],"#sigma_{M26} = 4.0 #mum","l");
  leg2->AddEntry(pre[4],"#sigma_{M26} = 5.0 #mum","l");
  //leg2->AddEntry(pre[5],"#sigma_{M26} = 6.0 #mum","l");
  leg2->SetTextSize(0.03);
  leg2->Draw();

  // Create third legend with some additional information for the plot
  TLegend *leg3 = new TLegend(0.7,0.10,0.95,0.30);
  leg3->SetBorderSize(0);
  leg3->SetFillColor(0);
  leg3->SetFillStyle(0);
  char beambuffer[100]="";
  sprintf(beambuffer, "DESY %d GeV e^{-}", nominalbeam);
  leg3->AddEntry((TObject*)0, beambuffer, "");
  TString threshlabel = "M26 threshold: ";
  threshlabel += global_thresh;
  leg3->AddEntry((TObject*)0, threshlabel, "");
  char chi2buffer[100]="";
  sprintf(chi2buffer, "The #Chi^{2} of the fit is %4.2f", (chimin2D/12.0));
  //leg3->AddEntry((TObject*)0, chi2buffer, "");
  leg3->AddEntry(x_direction,"x direction","p");
  leg3->AddEntry(y_direction,"y direction","p");
  leg3->SetTextSize(0.030);
  leg3->Draw();

  // Output
  smilie->Print("pics/"+outputname+"_smilie.eps");
  _outputFile->cd();
  smilie->Write();
  smilie->Close();
   */

  }

// Fitting of each file -> noise and efficiency
void noise(Int_t runnumber)
{
  TString filename("../filteredhistos/run00");
  if (runnumber <= 999)
    filename+="0";
  if (runnumber <= 99)
    filename+="0";
  filename+=runnumber;
  filename+="-fitter.root";
  TFile *f6 = new TFile(filename);
  TH1D *h_m26[12];
  TH1D *h_m26_eff[12];
  f6->cd();

  // Load file
  h_m26[0] = (TH1D*) (f6->Get("DUTHisto00/DUTnoiseX"))->Clone();
  h_m26[1] = (TH1D*) (f6->Get("DUTHisto00/DUTnoiseY"))->Clone();
  h_m26[2] = (TH1D*) (f6->Get("DUTHisto01/DUTnoiseX"))->Clone();
  h_m26[3] = (TH1D*) (f6->Get("DUTHisto01/DUTnoiseY"))->Clone();
  h_m26[4] = (TH1D*) (f6->Get("DUTHisto02/DUTnoiseX"))->Clone();
  h_m26[5] = (TH1D*) (f6->Get("DUTHisto02/DUTnoiseY"))->Clone();
  h_m26[6] = (TH1D*) (f6->Get("DUTHisto03/DUTnoiseX"))->Clone();
  h_m26[7] = (TH1D*) (f6->Get("DUTHisto03/DUTnoiseY"))->Clone();
  h_m26[8] = (TH1D*) (f6->Get("DUTHisto04/DUTnoiseX"))->Clone();
  h_m26[9] = (TH1D*) (f6->Get("DUTHisto04/DUTnoiseY"))->Clone();
  h_m26[10] = (TH1D*) (f6->Get("DUTHisto05/DUTnoiseX"))->Clone();
  h_m26[11] = (TH1D*) (f6->Get("DUTHisto05/DUTnoiseY"))->Clone();

  h_m26_eff[0] = (TH1D*) (f6->Get("DUTHisto00/DUTeffiX"))->Clone();
  h_m26_eff[1] = (TH1D*) (f6->Get("DUTHisto00/DUTeffiY"))->Clone();
  h_m26_eff[2] = (TH1D*) (f6->Get("DUTHisto01/DUTeffiX"))->Clone();
  h_m26_eff[3] = (TH1D*) (f6->Get("DUTHisto01/DUTeffiY"))->Clone();
  h_m26_eff[4] = (TH1D*) (f6->Get("DUTHisto02/DUTeffiX"))->Clone();
  h_m26_eff[5] = (TH1D*) (f6->Get("DUTHisto02/DUTeffiY"))->Clone();
  h_m26_eff[6] = (TH1D*) (f6->Get("DUTHisto03/DUTeffiX"))->Clone();
  h_m26_eff[7] = (TH1D*) (f6->Get("DUTHisto03/DUTeffiY"))->Clone();
  h_m26_eff[8] = (TH1D*) (f6->Get("DUTHisto04/DUTeffiX"))->Clone();
  h_m26_eff[9] = (TH1D*) (f6->Get("DUTHisto04/DUTeffiY"))->Clone();
  h_m26_eff[10] = (TH1D*) (f6->Get("DUTHisto05/DUTeffiX"))->Clone();
  h_m26_eff[11] = (TH1D*) (f6->Get("DUTHisto05/DUTeffiY"))->Clone();



  // Add the values for all 6 sensor planes, in x and y and divide by 12

  Double_t noisevalue = 0;
  Double_t noisevalue_error = 0;
  Double_t efficiency = 0;
  Double_t efficiency_error = 0;
  for(int i=4;i<8;i++)
  {
    noisevalue += h_m26[i]->GetMean(2);
    noisevalue_error += h_m26[i]->GetMeanError(2);
    efficiency += h_m26_eff[i]->GetMean(2);
    efficiency_error += h_m26_eff[i]->GetMeanError(2);
  }
  avgnoise = noisevalue / 4.0;
  avgnoise_error = noisevalue_error / 4.0;
  avgeffi = efficiency / 4.0;
  avgeffi_error = efficiency_error / 4.0;




  /*
     Double_t noisevalue = 0;
     Double_t noisevalue_error = 0;
     Double_t efficiency = 0;
     Double_t efficiency_error = 0;

     noisevalue = h_m26[4]->GetMean(2);
     noisevalue_error = h_m26[4]->GetMeanError(2);
     efficiency = h_m26_eff[4]->GetMean(2);
     efficiency_error = h_m26_eff[4]->GetMeanError(2);

     avgnoise = noisevalue;
     avgnoise_error = noisevalue_error;
     avgeffi = efficiency;
     avgeffi_error = efficiency_error;
   */

}

void getclusize(Int_t runnumber)
{
  TString filename("../filteredhistos/run00");
  if (runnumber <= 999)
    filename+="0";
  if (runnumber <= 99)
    filename+="0";
  filename+=runnumber;
  filename+="-fitter.root";
  TFile *f6 = new TFile(filename);
  TH1D *h_m26[12];
  TH1D *h_m26_eff[12];
  f6->cd();

  // Load file
  h_m26[0] = (TH1D*) (f6->Get("DUTHisto00/clusterSizeX"))->Clone();
  h_m26[1] = (TH1D*) (f6->Get("DUTHisto00/clusterSizeY"))->Clone();
  h_m26[2] = (TH1D*) (f6->Get("DUTHisto01/clusterSizeX"))->Clone();
  h_m26[3] = (TH1D*) (f6->Get("DUTHisto01/clusterSizeY"))->Clone();
  h_m26[4] = (TH1D*) (f6->Get("DUTHisto02/clusterSizeX"))->Clone();
  h_m26[5] = (TH1D*) (f6->Get("DUTHisto02/clusterSizeY"))->Clone();
  h_m26[6] = (TH1D*) (f6->Get("DUTHisto03/clusterSizeX"))->Clone();
  h_m26[7] = (TH1D*) (f6->Get("DUTHisto03/clusterSizeY"))->Clone();
  h_m26[8] = (TH1D*) (f6->Get("DUTHisto04/clusterSizeX"))->Clone();
  h_m26[9] = (TH1D*) (f6->Get("DUTHisto04/clusterSizeY"))->Clone();
  h_m26[10] = (TH1D*) (f6->Get("DUTHisto05/clusterSizeX"))->Clone();
  h_m26[11] = (TH1D*) (f6->Get("DUTHisto05/clusterSizeY"))->Clone();

  // Add the values for all 6 sensor planes, in x and y and divide by 12
  Double_t clustervalue = 0;
  Double_t clustervalue_error = 0;
  for(int i=0;i<12;i++)
  {
    clustervalue += h_m26[i]->GetMean(1);
    clustervalue_error += h_m26[i]->GetMeanError(1);
  }
  avgclustersize = clustervalue / 12.0;
  avgclustersize_error = clustervalue_error / 12.0;
}


void getpointing(Int_t runnumber, float sigm26)
{
  TString filename("../filteredhistos/run00");
  if (runnumber <= 999)
    filename+="0";
  if (runnumber <= 99)
    filename+="0";
  filename+=runnumber;
  filename+="-fitter.root";
  TFile *f6 = new TFile(filename);
  TH1D *h_m26[12];
  TH1D *h_m26_eff[12];
  f6->cd();

  cout << "File is " << runnumber << endl;

  // Load file
  h_m26[0] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftX"))->Clone();
  h_m26[1] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftY"))->Clone();
  h_m26[2] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftX"))->Clone();
  h_m26[3] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftY"))->Clone();
  h_m26[4] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftX"))->Clone();
  h_m26[5] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftY"))->Clone();
  h_m26[6] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftX"))->Clone();
  h_m26[7] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftY"))->Clone();
  h_m26[8] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftX"))->Clone();
  h_m26[9] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftY"))->Clone();
  h_m26[10] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftX"))->Clone();
  h_m26[11] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftY"))->Clone();

  Double_t sigmeas = 0;
  Double_t sigmeas_error = 0;
  for (int i = 4; i<6; i++)
  {
    TF1 *f1 = new TF1("f1","gaus");
    TF1 *f2 = new TF1("f2","gaus");
    f2->SetLineWidth(2);
    f2->SetLineStyle(1);
    f2->SetLineColor(kBlack);
    h_m26[i]->Fit(f1,"EMQI0","");
    Double_t mean_1 = f1->GetParameter(1);
    Double_t sigma_1 = f1->GetParameter(2);



    // Repeat within 2 sigma
    h_m26[i]->Fit(f2,"EQMI","", (mean_1 - 2.0*sigma_1), (mean_1 + 2.0*sigma_1));
    Double_t mean_2 = f2->GetParameter(1);
    Double_t sigma_2 = f2->GetParameter(2);
    Double_t sigma_2e = f2->GetParError(2);

    cout << "measured width: " << sigma_2*1000.0 << endl;

    sigmeas += sigma_2*1000.0;
    sigmeas_error += sigma_2e*1000.0;
  }


  //sigm26 = 3.42;
  //Double_t sigm26_e = 0.12;
  Double_t sigm26_e = 0.035;
  /*
     if (runnumber == 291)
     {
     sigm26 = 3.55;
     cout << "asdfasdf" << endl; 
     }

     if (runnumber == 755)
     {
     sigm26 = 3.18;
     cout << "asdf" << endl;  

     }*/

  // These lines for pointing at plane 3:

  /*
     avgmeas = sqrt( (sigmeas / 2.0)*(sigmeas / 2.0) - sigm26*sigm26);
  //avgmeas_error = sqrt( (sigmeas_error / 2.0)*(sigmeas_error / 2.0) + sigm26_e*sigm26_e);
  avgmeas_error = sqrt( (sigmeas_error / 2.0 * sigmeas / 2.0 / avgmeas)*(sigmeas_error / 2.0 * sigmeas / 2.0 / avgmeas) + (sigm26*sigm26_e/avgmeas)*(sigm26*sigm26_e/avgmeas)     );

  cout << "sigmeas is " << sigmeas/2.0 << " pm " << sigmeas_error/2.0 << endl;
  cout << "pointing res at plane 3 is " << avgmeas << " pm  " << avgmeas_error << endl;
  cout << endl;
   */

  // these lines for extrapolation to telescope center:


  float kfive = 0.2209302326;

  avgmeas = sqrt( (sigmeas / 2.0)*(sigmeas / 2.0) + sigm26*sigm26*((1.0/6.0) - kfive - 1.0) );
  avgmeas_error = sqrt((sigmeas/2.0*sigmeas_error/2.0/avgmeas)*(sigmeas/2.0*sigmeas_error/2.0/avgmeas) + (sigm26*sigm26_e/avgmeas*((1.0/6.0) - kfive - 1.0))*(sigm26*sigm26_e/avgmeas*((1.0/6.0) - kfive - 1.0)));

  cout << "sigmeas is " << sigmeas/2.0 << " pm " << sigmeas_error/2.0 << endl;
  cout << "pointing res at center telescope is " << avgmeas << " pm  " << avgmeas_error << endl;
  cout << endl;


}


void histoplot(Int_t runnumber)
{
  TString filename("../filteredhistos/run00");
  if (runnumber <= 999)
    filename+="0";
  if (runnumber <= 99)
    filename+="0";
  filename+=runnumber;
  filename+="-fitter.root";
  TFile *f6 = new TFile(filename);
  TH1D *h_m26[12];
  TH1D *h_m26_eff[12];
  f6->cd();

  cout << "File is " << runnumber << endl;

  // Load file
  h_m26[0] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftX"))->Clone();
  h_m26[1] = (TH1D*) (f6->Get("DUTHisto00/DUTshiftY"))->Clone();
  h_m26[2] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftX"))->Clone();
  h_m26[3] = (TH1D*) (f6->Get("DUTHisto01/DUTshiftY"))->Clone();
  h_m26[4] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftX"))->Clone();
  h_m26[5] = (TH1D*) (f6->Get("DUTHisto02/DUTshiftY"))->Clone();
  h_m26[6] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftX"))->Clone();
  h_m26[7] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftY"))->Clone();
  h_m26[8] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftX"))->Clone();
  h_m26[9] = (TH1D*) (f6->Get("DUTHisto04/DUTshiftY"))->Clone();
  h_m26[10] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftX"))->Clone();
  h_m26[11] = (TH1D*) (f6->Get("DUTHisto05/DUTshiftY"))->Clone();

  double val[12];


  for (int i = 0; i<12; i++)
  {
    TF1 *f1 = new TF1("f1","gaus");
    TF1 *f2 = new TF1("f2","gaus");
    f2->SetLineWidth(2);
    f2->SetLineStyle(1);
    f2->SetLineColor(kBlack);
    h_m26[i]->Fit(f1,"EMQI0","");
    Double_t mean_1 = f1->GetParameter(1);
    Double_t sigma_1 = f1->GetParameter(2);



    // Repeat within 2 sigma
    h_m26[i]->Fit(f2,"EQMI","", (mean_1 - 2.0*sigma_1), (mean_1 + 2.0*sigma_1));
    Double_t mean_2 = f2->GetParameter(1);
    Double_t sigma_2 = f2->GetParameter(2);
    Double_t sigma_2e = f2->GetParError(2);

    cout << "measured width: " << sigma_2*1000.0 << endl;

    val[i] = sigma_2*1000.0;

    /*
       if (i==0 || i == 1 || i == 10 || i == 11)
       {
       delta0->Fill(sigma_2*1000.0);
       }
       if (i==2 || i == 3 || i == 8 || i == 9)
       {
       delta1->Fill(sigma_2*1000.0);
       }
       if (i==4 || i == 5 || i == 6 || i == 7)
       {
       delta2->Fill(sigma_2*1000.0);
       }*/

  }

  delta0->Fill(fabs(val[0]-val[1]));
  delta1->Fill(fabs(val[2]-val[3]));
  delta2->Fill(fabs(val[4]-val[5]));
  delta3->Fill(fabs(val[6]-val[7]));
  delta4->Fill(fabs(val[8]-val[9]));
  delta5->Fill(fabs(val[10]-val[11]));





  delta6->Fill(fabs(val[0]-val[10]));
  delta7->Fill(fabs(val[2]-val[8]));
  delta8->Fill(fabs(val[4]-val[6]));
  delta9->Fill(fabs(val[1]-val[11]));
  delta10->Fill(fabs(val[3]-val[9]));
  delta11->Fill(fabs(val[5]-val[7]));





}



// Let's go!
int main()
{

  gROOT->SetBatch();
  //gSystem->Load("libMinuit");
  //gSystem->Load("lib/libGBL.so");
  //gSystem->AddIncludePath("include/");


  _outputFile = new TFile("output.root", "RECREATE");
  _outputFile->cd();


  // Initial telescope constructor name, AnaTel.geom for 150mm, AnaTel_thin.geom for 20mm data
  //char telescopebuild[50];
  //int tempint = 0;

  //telescopebuild = "AnaTel.geom";


  // Initialises the run vectors
  std::vector<int> run2;
  std::vector<int> run3;
  std::vector<int> run4;
  std::vector<int> run5;
  std::vector<int> run2th;
  std::vector<int> run3th;
  std::vector<int> run5th;
  std::vector<int> run120;

  std::vector<int> testing;

  testing.push_back(291);
  testing.push_back(294);
  testing.push_back(295);
  testing.push_back(296);
  testing.push_back(297);
  testing.push_back(298);
  testing.push_back(299);
  testing.push_back(300);

  // Fills the run vectors with runnumbers
  if(true)
  {
    // energy 2 GeV:
    //run2.push_back(5037);	// thr 3
    run2.push_back(5041);	// thr 3
    run2.push_back(5040);	// thr 4
    run2.push_back(5041);	// thr 5
    run2.push_back(5043);	// thr 6
    run2.push_back(5045);	// thr 7
    run2.push_back(5049);	// thr 8
    run2.push_back(5051);	// thr 9
    run2.push_back(5052);	// thr 10
    run2.push_back(5054);	// thr 11
    run2.push_back(5056);	// thr 12

    // energy 3 GeV:
    //run3.push_back(5066);	// thr 3
    run3.push_back(5065);	// thr 4
    run3.push_back(5065);	// thr 4
    run3.push_back(5064);	// thr 5
    run3.push_back(5063);	// thr 6
    run3.push_back(5062);	// thr 7
    run3.push_back(5061);	// thr 8
    run3.push_back(5060);	// thr 9
    run3.push_back(5059);	// thr 10
    run3.push_back(5058);	// thr 11
    run3.push_back(5057);	// thr 12

    // energy 4.4 GeV:
    run4.push_back(103);	// thr 3
    run4.push_back(101);	// thr 4
    run4.push_back(99);		// thr 5
    run4.push_back(98);		// thr 6
    run4.push_back(97);		// thr 7
    run4.push_back(96);		// thr 8
    run4.push_back(94);		// thr 9
    run4.push_back(92);		// thr 10
    run4.push_back(90);		// thr 11
    run4.push_back(88);		// thr 12

    // energy 5 GeV:
    run5.push_back(85);		// thr 3
    run5.push_back(83);		// thr 4
    run5.push_back(78);		// thr 5
    run5.push_back(37);		// thr 6
    run5.push_back(76);		// thr 7
    run5.push_back(75);		// thr 8
    run5.push_back(72);		// thr 9
    run5.push_back(71);		// thr 10
    run5.push_back(68);		// thr 11
    run5.push_back(64);		// thr 12

    // energy 2 GeV, 20mm Data:
    run2th.push_back(233);	// thr 3
    run2th.push_back(234);	// thr 4
    run2th.push_back(235);	// thr 5
    run2th.push_back(236);	// thr 6
    run2th.push_back(237);	// thr 7
    run2th.push_back(238);	// thr 8
    run2th.push_back(239);	// thr 9
    run2th.push_back(240);	// thr 10
    run2th.push_back(241);	// thr 11
    run2th.push_back(245);	// thr 12

    // energy 3 GeV, 20mm Data:
    run3th.push_back(246);	// thr 3
    run3th.push_back(247);	// thr 4
    run3th.push_back(248);	// thr 5
    run3th.push_back(249);	// thr 6
    run3th.push_back(250);	// thr 7
    run3th.push_back(251);	// thr 8
    run3th.push_back(252);	// thr 9
    run3th.push_back(253);	// thr 10
    run3th.push_back(254);	// thr 11
    run3th.push_back(255);	// thr 12

    // energy 5 GeV, 20mm Data:
    run5th.push_back(258);	// thr 3
    run5th.push_back(259);	// thr 4
    run5th.push_back(261);	// thr 5
    run5th.push_back(262);	// thr 6
    run5th.push_back(263);	// thr 7
    run5th.push_back(264);	// thr 8
    run5th.push_back(265);	// thr 9
    run5th.push_back(266);	// thr 10
    run5th.push_back(267);	// thr 11
    run5th.push_back(268);	// thr 12

    // energy 120 GeV CERN, 150mm Data:
    run120.push_back(752);	// thr 3
    run120.push_back(753);	// thr 4
    run120.push_back(754);	// thr 5
    run120.push_back(755);	// thr 6
    run120.push_back(756);	// thr 7
    run120.push_back(757);	// thr 8
    run120.push_back(758);	// thr 9
    run120.push_back(760);	// thr 10        alternative: 759
    run120.push_back(761);	// thr 11
    run120.push_back(762);	// thr 12

  }

  // Runmode: 0 for all, 1 for clustersize, 2 for threshold, 3 for threshold and clustersize, 4 for noise, 5 for E plot, 9 for testing...
  // 6 for clustersie only, 7 for pointing
  Int_t runmode = 9;

  if (runmode==9)
  {
    global_thresh=6;



    telescopebuild = "AnaTel.geom";
    planedistance = 150;
    for(int j=0;j<nplanes;j++)
      posx[j] = planedistance*j;

    fitter(98,5.4);
    fitter(37,6.);
    //fitter(5063,3.6);
    //fitter(5043,2.4);

    /*
       fitter(37,5.0);
       fitter(37,5.1);
       fitter(37,5.2);
       fitter(37,5.3);
       fitter(37,5.4);
       fitter(37,5.5);
       fitter(37,5.6);
       fitter(37,5.7);
       fitter(37,5.8);
       fitter(37,5.9);
       fitter(37,6.0);
     */
    /*    
	  fitter(5063,3.3);
	  fitter(5063,3.4);
	  fitter(5063,3.5);
	  fitter(5063,3.6);
	  fitter(5063,3.7);
	  fitter(5063,3.8);
	  fitter(5063,3.9);
	  fitter(5063,4.0);
	  fitter(5063,4.1);
	  fitter(5063,4.2);
	  fitter(5063,4.3);
     */
    /*
       fitter(37,6.1);
       fitter(37,6.2);
       fitter(98,5.5);
       fitter(98,5.6);

     */

    telescopebuild = "AnaTel_thin.geom";
    planedistance = 20;
    for(int j=0;j<nplanes;j++)
      posx[j] = planedistance*j;

    fitter(236,2.4);
    fitter(262,6.0);


    /*   fitter(296,12.5);
	 fitter(297,12.5);
	 fitter(298,12.5);
	 fitter(299,12.5);
	 fitter(300,12.5);*/
    /*
       fitter(testing[1],12.5);
       fitter(testing[2],12.5);
       fitter(testing[3],12.5);
       fitter(testing[4],12.5);
       fitter(testing[5],12.5);
       fitter(testing[6],12.5);
       fitter(testing[7],12.5);*/
    /*
       telescopebuild = "AnaTel_thin.geom";
       planedistance = 20;
       for(int j=0;j<nplanes;j++)
       posx[j] = planedistance*j;


       fitter(testing[1],5.0);
     */
  }





  // Simply run over everything
  if(runmode == 0)
  {
    cout << " " << endl;
    cout << "Mode 0" << endl;
    cout << " " << endl;
    cout << "Running over all runs" << endl;
    cout << " " << endl;
    cout << "Wide Geometry" << endl;
    cout << " " << endl;

    telescopebuild = "AnaTel.geom";
    planedistance = 150;
    for(int j=0;j<nplanes;j++)
      posx[j] = 150.0*j;

    for(int i=0;i<run2.size();i++)
    {
      global_thresh = i+3;
      fitter( run2[i], 2.0 );
    }
    for(int i=0;i<run3.size();i++)
    {
      global_thresh = i+3;
      fitter( run3[i], 3.0 );
    }
    for(int i=0;i<run4.size();i++)
    {
      global_thresh = i+3;
      fitter( run4[i], 4.4 );
    }
    for(int i=0;i<run5.size();i++)
    {
      global_thresh = i+3;
      fitter( run5[i], 5.0 );
    }

    for(int i=0;i<run120.size();i++)
    {
      global_thresh = i+3;
      fitter( run120[i], 120.0 );
    }

    cout << " " << endl;
    cout << "Thin Geometry" << endl;
    cout << " " << endl;

    telescopebuild = "AnaTel_thin.geom";
    planedistance = 20;
    for(int j=0;j<nplanes;j++)
      posx[j] = 20.0*j;

    for(int i=0;i<run2th.size();i++)
    {
      global_thresh = i+3;
      fitter( run2th[i], 2.0 );
    }
    for(int i=0;i<run3th.size();i++)
    {
      global_thresh = i+3;
      fitter( run3th[i], 3.0 );
    }
    for(int i=0;i<run5th.size();i++)
    {
      global_thresh = i+3;
      fitter( run5th[i], 5.0 );
    }
  }

  // Run according to clustersize, this should be at one threshold -> to be done -> Geometry needs implementing
  if(runmode == 1)
  {
    cout << "Running over clustersize" << endl;
    // Smallest clustersize is 1, largest is clustercount -1
    const Int_t clustercount = 5;
    Double_t clusterresult[clustercount];
    Double_t clustererror[clustercount];
    Double_t x[clustercount];
    for (Int_t j=1;j<clustercount;j++)
    {
      x[j] = j+0.0;
      ostringstream convert;
      convert << j;
      submask = convert.str();

      // this needs to be done:
      for(int i=0;i<run3.size();i++)
	fitter( run3[i], 3 );
      clusterresult[j] = m26_resolution*1000.0;
      clustererror[j] = m26_res_error*1000.0;
    }

    // As a comparison: put "clustersize = 0" as the normal uncut run.
    submask = "";
    for(int i=0;i<run3.size();i++)
      fitter( run3[i], 3 );
    clusterresult[0] = m26_resolution*1000.0;
    clustererror[0] = m26_res_error*1000.0;

    TCanvas *c1 = new TCanvas("c1","Resolution vs. Clustersize",10,10,800,600);
    c1->SetFillColor(0);
    c1->SetGrid();
    TGraph *gr = new TGraph(clustercount,x,clusterresult);
    gr->SetLineColor(kBlack);
    gr->SetLineWidth(2);
    gr->SetMarkerColor(kRed);
    gr->SetMarkerStyle(21);
    gr->SetTitle("Resolution vs. Clustersize");
    gr->GetXaxis()->SetTitle("Clustersize");
    gr->GetYaxis()->SetTitle("#sigma_{M26} in #mum");
    gr->Draw("ACP");
    c1->Update();
    //   c1->GetFrame()->SetFillColor(0);
    //    c1->GetFrame()->SetBorderSize(0);
    c1->Modified();
    c1->Print("pics/clustersize_somename.eps");
    _outputFile->cd();
    c1->Write();
    c1->Close();
  }

  // Run over thresholds
  if(runmode == 2)
  {
    cout << "Threshold mode" << endl;
    const Int_t threshcount = 12;
    Double_t threshresult2[threshcount];
    Double_t threshresult3[threshcount];
    Double_t threshresult4[threshcount];
    Double_t threshresult5[threshcount];
    Double_t threshresult120[threshcount];
    Double_t thresherror2[threshcount];
    Double_t thresherror3[threshcount];
    Double_t thresherror4[threshcount];
    Double_t thresherror5[threshcount];
    Double_t thresherror120[threshcount];

    Double_t threshresult2th[threshcount];
    Double_t threshresult3th[threshcount];
    Double_t threshresult5th[threshcount];
    Double_t thresherror2th[threshcount];
    Double_t thresherror3th[threshcount];
    Double_t thresherror5th[threshcount];

    Double_t x[threshcount];
    Double_t xerrorthresh[threshcount];
    for (Int_t j=0;j<(threshcount-2);j++)
    {

      telescopebuild = "AnaTel.geom";
      planedistance = 150;
      for(int jj=0;jj<nplanes;jj++)
	posx[jj] = planedistance*jj;


      x[j] = j+3;
      global_thresh = x[j];
      xerrorthresh[j] = 0.0;
      submask = "";


      fitter( run2[j], 2.2 ); 
      threshresult2[j] = m26_resolution*1000.0;
      thresherror2[j] = global_plot_error;



      fitter( run3[j], 3.3 );
      threshresult3[j] = m26_resolution*1000.0;
      thresherror3[j] = global_plot_error;



      fitter( run4[j], 4.84 );
      threshresult4[j] = m26_resolution*1000.0;
      thresherror4[j] = global_plot_error;
      fitter( run5[j], 5.5 );
      threshresult5[j] = m26_resolution*1000.0;
      thresherror5[j] = global_plot_error;

      /*
	 fitter( run120[j], 120.0 );
	 threshresult120[j] = m26_resolution*1000.0;
	 thresherror120[j] = global_plot_error;





         telescopebuild = "AnaTel_thin.geom";
	 planedistance = 20;
	 for(int jj=0;jj<nplanes;jj++)
	 posx[jj] = planedistance*jj;



	 fitter( run2th[j], 2.0 );
	 threshresult2th[j] = m26_resolution*1000.0;
	 thresherror2th[j] = global_plot_error;

	 fitter( run3th[j], 3.0 );
	 threshresult3th[j] = m26_resolution*1000.0;
	 thresherror3th[j] = global_plot_error;

	 fitter( run5th[j], 5.0 );
	 threshresult5th[j] = m26_resolution*1000.0;
	 thresherror5th[j] = global_plot_error;

       */
    }

    TGraphErrors *gr2 = new TGraphErrors((threshcount-2),x,threshresult2,xerrorthresh,thresherror2);
    TGraphErrors *gr3 = new TGraphErrors((threshcount-2),x,threshresult3,xerrorthresh,thresherror3);
    TGraphErrors *gr4 = new TGraphErrors((threshcount-2),x,threshresult4,xerrorthresh,thresherror4);
    TGraphErrors *gr5 = new TGraphErrors((threshcount-2),x,threshresult5,xerrorthresh,thresherror5);
    TGraphErrors *gr120 = new TGraphErrors((threshcount-2),x,threshresult120,xerrorthresh,thresherror120);

    TGraphErrors *gr2th = new TGraphErrors((threshcount-2),x,threshresult2th,xerrorthresh,thresherror2th);
    TGraphErrors *gr3th = new TGraphErrors((threshcount-2),x,threshresult3th,xerrorthresh,thresherror3th);
    TGraphErrors *gr5th = new TGraphErrors((threshcount-2),x,threshresult5th,xerrorthresh,thresherror5th);

    TH1D *h_axis = new TH1D("th_axis","th_axis",1, 2.0, 13.0);
    TCanvas *threshold = new TCanvas("threshold","threshold",10,10,800,600);
    gStyle->SetPadBorderMode(0);
    gStyle->SetOptStat(0);
    threshold->SetFillColor(0);
    threshold->Divide(1,1);
    threshold->cd(1);
    gStyle->SetErrorX(0);
    gPad->SetLogx(0);
    gPad->SetLogy(0);
    gPad->SetGridx();
    gPad->SetGridy();

    gr2->SetMarkerStyle(20);
    gr2->SetMarkerColor(kBlack);
    gr2->SetMarkerSize(3);
    gr3->SetMarkerStyle(21);
    gr3->SetMarkerColor(kGreen);
    gr3->SetMarkerSize(3);
    gr4->SetMarkerStyle(22);
    gr4->SetMarkerColor(kRed);
    gr4->SetMarkerSize(3);
    gr5->SetMarkerStyle(23);
    gr5->SetMarkerColor(kBlue);
    gr5->SetMarkerSize(3);
    gr120->SetMarkerStyle(23);
    gr120->SetMarkerColor(kOrange);
    gr120->SetMarkerSize(3);

    gr2th->SetMarkerStyle(24);
    gr2th->SetMarkerColor(kBlack);
    gr2th->SetMarkerSize(3);
    gr3th->SetMarkerStyle(25);
    gr3th->SetMarkerColor(kGreen);
    gr3th->SetMarkerSize(3);
    gr5th->SetMarkerStyle(26);
    gr5th->SetMarkerColor(kBlue);
    gr5th->SetMarkerSize(3);

    histo_cfg(h_axis, "Threshold (s/n)","#sigma_{M26} (#mum)","");
    h_axis->SetMinimum(0.0);
    h_axis->SetMaximum(5.0);
    h_axis->Draw("hist");
    gr2->Draw("P");
    gr3->Draw("P");
    gr4->Draw("P");
    gr5->Draw("P");
    gr120->Draw("P");

    gr2th->Draw("P");
    gr3th->Draw("P");
    gr5th->Draw("P");


    TLegend *leg = new TLegend(0.59,0.55,0.90,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Energy:");
    leg->AddEntry(gr2,"p = 2 GeV, #Delta_{z} = 150 mm","p");
    leg->AddEntry(gr3,"p = 3 GeV, #Delta_{z} = 150 mm","p");
    leg->AddEntry(gr4,"p = 4.4 GeV, #Delta_{z} = 150 mm","p");
    leg->AddEntry(gr5,"p = 5 GeV, #Delta_{z} = 150 mm","p");
    leg->AddEntry(gr120,"p = 120 GeV, #Delta_{z} = 150 mm","p");

    leg->AddEntry(gr2th,"p = 2 GeV, #Delta_{z} = 20 mm","p");
    leg->AddEntry(gr3th,"p = 3 GeV, #Delta_{z} = 20 mm","p");
    leg->AddEntry(gr5th,"p = 5 GeV, #Delta_{z} = 20 mm","p");

    leg->Draw();

    // Output
    threshold->Print("pics/threshold_all.eps");
    _outputFile->cd();
    threshold->Write();
    threshold->Close();
  }

  // Run over thresholds (x) and clustersize (y), res is z, geometry needs reset and thin implementing
  if(runmode == 3)
  {
    cout << "Threshold and clustersize mode" << endl;
    const Int_t threshcount = 12;
    const Int_t clustercount = 5;
    Double_t threshresult2[threshcount][clustercount];
    Double_t threshresult3[threshcount][clustercount];
    Double_t threshresult4[threshcount][clustercount];
    Double_t threshresult5[threshcount][clustercount];
    Double_t threshresulttotal[threshcount][clustercount];
    Double_t thresherror2[threshcount][clustercount];
    Double_t thresherror3[threshcount][clustercount];
    Double_t thresherror4[threshcount][clustercount];
    Double_t thresherror5[threshcount][clustercount];
    Double_t thresherrortotal[threshcount][clustercount];
    Double_t x[threshcount][clustercount];
    Double_t xerrorthresh[threshcount][clustercount];
    Double_t y[threshcount][clustercount];
    Double_t yerrorcluster[threshcount][clustercount];

    for (Int_t i=0;i<clustercount;i++)
    {
      for (Int_t j=0;j<(threshcount-2);j++)
      {
	cout << "I am now at clustersize " << i << " and threshold " << j+3 << " !" << endl;
	y[j][i] = i;
	x[j][i] = j+3;
	global_thresh = x[j][i];
	xerrorthresh[j][i] = 0.0;
	ostringstream convert;
	convert << i;
	submask = convert.str();
	if (i == 0)
	  submask = "";
	/*	fitter( run2[j], 2.0 );
		threshresult2[j][i] = m26_resolution*1000.0;
		thresherror2[j][i] = m26_res_error*1000.0;
		fitter( run3[j], 3.0 );
		threshresult3[j][i] = m26_resolution*1000.0;
		thresherror3[j][i] = m26_res_error*1000.0; */
	fitter( run4[j], 4.4 );
	threshresult4[j][i] = m26_resolution*1000.0;
	thresherror4[j][i] = m26_res_error*1000.0;
	/*	fitter( run5[j], 5.0 );
		threshresult5[j][i] = m26_resolution*1000.0;
		thresherror5[j][i] = m26_res_error*1000.0; */
      }
    }

    /*   TCanvas *c_thr_vs_clu_2 = new TCanvas("thr_vs_clu_2", "Threshold vs. Clustersize", 600, 400);
	 c_thr_vs_clu_2->cd();
	 TH2D *hist2D_thr_vs_clu_2 = new TH2D("hist2D_thr_vs_clu_2", "Histo_thr_vs_clu_2", 10, 3., 13., 5, 0., 5.);
	 for(Int_t i=0;i<clustercount;i++)
	 {
	 for(Int_t j=0;j<(threshcount-2);j++)
	 {
	 hist2D_thr_vs_clu_2->Fill((j+3),i,threshresult2[j][i]);
	 }
	 }
	 hist2D_thr_vs_clu_2->GetXaxis()->SetTitle("Threshold");
	 hist2D_thr_vs_clu_2->GetYaxis()->SetTitle("Clustersize");
	 hist2D_thr_vs_clu_2->GetZaxis()->SetTitle("#sigma_{M26}");
	 hist2D_thr_vs_clu_2->Draw("LEGO2");
	 c_thr_vs_clu_2->Print("pics/thr_vs_clu_2.eps");
	 c_thr_vs_clu_2->Write();

	 TCanvas *c_thr_vs_clu_3 = new TCanvas("thr_vs_clu_3", "Threshold vs. Clustersize", 600, 400);
	 c_thr_vs_clu_3->cd();
	 TH2D *hist2D_thr_vs_clu_3 = new TH2D("hist2D_thr_vs_clu_3", "Histo_thr_vs_clu_3", 10, 3., 13., 5, 0., 5.);
	 for(Int_t i=0;i<clustercount;i++)
	 {
	 for(Int_t j=0;j<(threshcount-2);j++)
	 {
	 hist2D_thr_vs_clu_3->Fill((j+3),i,threshresult3[j][i]);
	 }
	 }
	 hist2D_thr_vs_clu_3->GetXaxis()->SetTitle("Threshold");
	 hist2D_thr_vs_clu_3->GetYaxis()->SetTitle("Clustersize");
	 hist2D_thr_vs_clu_3->GetZaxis()->SetTitle("#sigma_{M26}");
	 hist2D_thr_vs_clu_3->Draw("LEGO2");
	 c_thr_vs_clu_3->Print("pics/thr_vs_clu_3.eps");
	 c_thr_vs_clu_3->Write();
     */
    TCanvas *c_thr_vs_clu_4 = new TCanvas("thr_vs_clu_4", "Threshold vs. Clustersize", 600, 400);
    c_thr_vs_clu_4->cd();
    TH2D *hist2D_thr_vs_clu_4 = new TH2D("hist2D_thr_vs_clu_4", "Histo_thr_vs_clu_4", 10, 3., 13., 5, 0., 5.);
    for(Int_t i=0;i<clustercount;i++)
    {
      for(Int_t j=0;j<(threshcount-2);j++)
      {
	hist2D_thr_vs_clu_4->Fill((j+3),i,threshresult4[j][i]);
      }
    }
    hist2D_thr_vs_clu_4->GetXaxis()->SetTitle("Threshold");
    hist2D_thr_vs_clu_4->GetYaxis()->SetTitle("Clustersize");
    hist2D_thr_vs_clu_4->GetZaxis()->SetTitle("#sigma_{M26}");
    hist2D_thr_vs_clu_4->Draw("LEGO2");
    c_thr_vs_clu_4->Print("pics/thr_vs_clu_4.eps");
    c_thr_vs_clu_4->Write();

    /*  TCanvas *c_thr_vs_clu_5 = new TCanvas("thr_vs_clu_5", "Threshold vs. Clustersize", 600, 400);
	c_thr_vs_clu_5->cd();
	TH2D *hist2D_thr_vs_clu_5 = new TH2D("hist2D_thr_vs_clu_5", "Histo_thr_vs_clu_5", 10, 3., 13., 5, 0., 5.);
	for(Int_t i=0;i<clustercount;i++)
	{
	for(Int_t j=0;j<(threshcount-2);j++)
	{
	hist2D_thr_vs_clu_5->Fill((j+3),i,threshresult5[j][i]);
	}
	}
	hist2D_thr_vs_clu_5->GetXaxis()->SetTitle("Threshold");
	hist2D_thr_vs_clu_5->GetYaxis()->SetTitle("Clustersize");
	hist2D_thr_vs_clu_5->GetZaxis()->SetTitle("#sigma_{M26}");
	hist2D_thr_vs_clu_5->Draw("LEGO2");
	c_thr_vs_clu_5->Print("pics/thr_vs_clu_5.eps");
	c_thr_vs_clu_5->Write();



	TCanvas *c_thr_vs_clu_total = new TCanvas("thr_vs_clu_total", "Threshold vs. Clustersize", 600, 400);
	c_thr_vs_clu_total->cd();
	TH2D *hist2D_thr_vs_clu_total = new TH2D("hist2D_thr_vs_clu_total", "Histo_thr_vs_clu_total", 10, 3., 13., 5, 0., 5.);
	for(Int_t i=0;i<clustercount;i++)
	{
	for(Int_t j=0;j<(threshcount-2);j++)
	{
	threshresulttotal[j][i]=((threshresult2[j][i]+threshresult3[j][i]+threshresult4[j][i]+threshresult5[j][i])/4.);
	hist2D_thr_vs_clu_total->Fill((j+3),i,threshresulttotal[j][i]);
	}
	}
	hist2D_thr_vs_clu_total->GetXaxis()->SetTitle("Threshold");
	hist2D_thr_vs_clu_total->GetYaxis()->SetTitle("Clustersize");
	hist2D_thr_vs_clu_total->GetZaxis()->SetTitle("#sigma_{M26}");
	hist2D_thr_vs_clu_total->Draw("LEGO2");
	c_thr_vs_clu_total->Print("pics/thr_vs_clu_total.eps");
	c_thr_vs_clu_total->Write();

     */
    /*

       TGraphErrors *gr2 = new TGraphErrors((threshcount-2),x,threshresult2,xerrorthresh,thresherror2);
       TGraphErrors *gr3 = new TGraphErrors((threshcount-2),x,threshresult3,xerrorthresh,thresherror3);
       TGraphErrors *gr4 = new TGraphErrors((threshcount-2),x,threshresult4,xerrorthresh,thresherror4);
       TGraphErrors *gr5 = new TGraphErrors((threshcount-2),x,threshresult5,xerrorthresh,thresherror5);
       TH1D *h_axis = new TH1D("th_axis","th_axis",1, 2.0, 13.0);
       TCanvas *threshold = new TCanvas("threshold","threshold",10,10,800,600);
       gStyle->SetPadBorderMode(0);
       gStyle->SetOptStat(0);
       threshold->SetFillColor(0);
       threshold->Divide(1,1);
       threshold->cd(1);
       gStyle->SetErrorX(0);
       gPad->SetLogx(0);
       gPad->SetLogy(0);
       gr2->SetMarkerStyle(20);
       gr2->SetMarkerColor(kBlack);
       gr2->SetMarkerSize(2);
       gr3->SetMarkerStyle(21);
       gr3->SetMarkerColor(kGreen);
       gr3->SetMarkerSize(2);
       gr4->SetMarkerStyle(22);
       gr4->SetMarkerColor(kRed);
       gr4->SetMarkerSize(2);
       gr5->SetMarkerStyle(23);
       gr5->SetMarkerColor(kBlue);
       gr5->SetMarkerSize(2);
       histo_cfg(h_axis, "Threshold (s/n)","#sigma_{M26} (#mum)","");
       th_axis->SetMinimum(0.0);
       th_axis->SetMaximum(5.0);
       th_axis->Draw("hist");
       gr2->Draw("P");
       gr3->Draw("P");
       gr4->Draw("P");
       gr5->Draw("P");

       TLegend *leg = new TLegend(0.59,0.55,0.90,0.85);
       leg->SetBorderSize(0);
       leg->SetFillColor(0);
       leg->SetFillStyle(0);
       leg->SetHeader("Energy:");
       leg->AddEntry(gr2,"p = 2 GeV","p");
       leg->AddEntry(gr3,"p = 3 GeV","p");
       leg->AddEntry(gr4,"p = 4.4 GeV","p");
       leg->AddEntry(gr5,"p = 5 GeV","p");
       leg->Draw();

    // Output
    threshold->Print("pics/threshold_all.eps"); */
  }

  // Run over threshold to get efficiency and noise -> 120gev missing
  if(runmode == 4)
  {

    // Results go in here
    Double_t threshnoise2[12];
    Double_t thresheffi2[12];
    Double_t threshnoise3[12];
    Double_t thresheffi3[12];
    Double_t threshnoise4[12];
    Double_t thresheffi4[12];
    Double_t threshnoise5[12];
    Double_t thresheffi5[12];
    Double_t threshnoise2_error[12];
    Double_t thresheffi2_error[12];
    Double_t threshnoise3_error[12];
    Double_t thresheffi3_error[12];
    Double_t threshnoise4_error[12];
    Double_t thresheffi4_error[12];
    Double_t threshnoise5_error[12];
    Double_t thresheffi5_error[12];
    Double_t thin_threshnoise2[12];
    Double_t thin_thresheffi2[12];
    Double_t thin_threshnoise3[12];
    Double_t thin_thresheffi3[12];
    Double_t thin_threshnoise5[12];
    Double_t thin_thresheffi5[12];
    Double_t thin_threshnoise2_error[12];
    Double_t thin_thresheffi2_error[12];
    Double_t thin_threshnoise3_error[12];
    Double_t thin_thresheffi3_error[12];
    Double_t thin_threshnoise5_error[12];
    Double_t thin_thresheffi5_error[12];
    Double_t x[12];
    Double_t xerror[12] = {0.0};

    cout << " " << endl;
    cout << "Mode 4" << endl;
    cout << " " << endl;
    cout << "Running over all runs - efficiency and noise" << endl;
    cout << " " << endl;
    cout << "Wide Geometry" << endl;
    cout << " " << endl;

    telescopebuild = "AnaTel.geom";
    planedistance = 150;
    for(int j=0;j<nplanes;j++)
      posx[j] = 150.0*j;

    for(int i=0;i<run2.size();i++)
    {
      noise( run2[i] );
      x[i] = i+3-0.15;
      threshnoise2[i] = avgnoise;
      thresheffi2[i] = avgeffi;
      threshnoise2_error[i] = avgnoise_error;
      thresheffi2_error[i] = avgeffi_error;
    }

    for(int i=0;i<run3.size();i++)
    {
      noise( run3[i] );
      x[i] = i+3-0.1;
      threshnoise3[i] = avgnoise;
      thresheffi3[i] = avgeffi;
      threshnoise3_error[i] = avgnoise_error;
      thresheffi3_error[i] = avgeffi_error;
    }

    for(int i=0;i<run4.size();i++)
    {
      noise( run4[i] );
      x[i] = i+3-0.05;
      threshnoise4[i] = avgnoise;
      thresheffi4[i] = avgeffi;
      threshnoise4_error[i] = avgnoise_error;
      thresheffi4_error[i] = avgeffi_error;
    }

    for(int i=0;i<run5.size();i++)
    {
      noise( run5[i] );
      x[i] = i+3;
      threshnoise5[i] = avgnoise;
      thresheffi5[i] = avgeffi;
      threshnoise5_error[i] = avgnoise_error;
      thresheffi5_error[i] = avgeffi_error;
    }

    cout << " " << endl;
    cout << "Thin Geometry" << endl;
    cout << " " << endl;

    telescopebuild = "AnaTel_thin.geom";
    planedistance = 20;

    for(int j=0;j<nplanes;j++)
      posx[j] = 20.0*j;

    for(int i=0;i<run2th.size();i++)
    {
      noise( run2th[i] );
      x[i] = i+3+0.05;
      thin_threshnoise2[i] = avgnoise;
      thin_thresheffi2[i] = avgeffi;
      thin_threshnoise2_error[i] = avgnoise_error;
      thin_thresheffi2_error[i] = avgeffi_error;
    }

    for(int i=0;i<run3th.size();i++)
    {
      noise( run3th[i] );
      x[i] = i+3+0.1;
      thin_threshnoise3[i] = avgnoise;
      thin_thresheffi3[i] = avgeffi;
      thin_threshnoise3_error[i] = avgnoise_error;
      thin_thresheffi3_error[i] = avgeffi_error;
    }

    for(int i=0;i<run5th.size();i++)
    {
      noise( run5th[i] );
      x[i] = i+3+0.15;
      thin_threshnoise5[i] = avgnoise;
      thin_thresheffi5[i] = avgeffi;
      thin_threshnoise5_error[i] = avgnoise_error;
      thin_thresheffi5_error[i] = avgeffi_error;
    }

    // Create graphs with the information
    TGraphErrors *gr2n = new TGraphErrors(10,x,threshnoise2,xerror,threshnoise2_error);
    TGraphErrors *gr2e = new TGraphErrors(10,x,thresheffi2,xerror,thresheffi2_error);
    TGraphErrors *gr3n = new TGraphErrors(10,x,threshnoise3,xerror,threshnoise3_error);
    TGraphErrors *gr3e = new TGraphErrors(10,x,thresheffi3,xerror,thresheffi3_error);
    TGraphErrors *gr4n = new TGraphErrors(10,x,threshnoise4,xerror,threshnoise4_error);
    TGraphErrors *gr4e = new TGraphErrors(10,x,thresheffi4,xerror,thresheffi4_error);
    TGraphErrors *gr5n = new TGraphErrors(10,x,threshnoise5,xerror,threshnoise5_error);
    TGraphErrors *gr5e = new TGraphErrors(10,x,thresheffi5,xerror,thresheffi5_error);

    TGraphErrors *gr2nth = new TGraphErrors(10,x,thin_threshnoise2,xerror,thin_threshnoise2_error);
    TGraphErrors *gr2eth = new TGraphErrors(10,x,thin_thresheffi2,xerror,thin_thresheffi2_error);
    TGraphErrors *gr3nth = new TGraphErrors(10,x,thin_threshnoise3,xerror,thin_threshnoise3_error);
    TGraphErrors *gr3eth = new TGraphErrors(10,x,thin_thresheffi3,xerror,thin_thresheffi3_error);
    TGraphErrors *gr5nth = new TGraphErrors(10,x,thin_threshnoise5,xerror,thin_threshnoise5_error);
    TGraphErrors *gr5eth = new TGraphErrors(10,x,thin_thresheffi5,xerror,thin_thresheffi5_error);

    // Let's plot this
    TH1D *h_axis = new TH1D("th_axis","th_axis",1, 2.0, 13.0);
    TCanvas *threshold = new TCanvas("threshold","threshold",10,10,800,600);
    gStyle->SetPadBorderMode(0);
    gStyle->SetOptStat(0);
    threshold->SetFillColor(0);
    threshold->Divide(1,1);
    threshold->cd(1);
    gStyle->SetErrorX(0);
    gPad->SetLogx(0);
    gPad->SetLogy(0);

    // Set apperance
    gr2n->SetMarkerStyle(5);
    gr2n->SetMarkerColor(kBlack);
    gr2n->SetMarkerSize(3);
    gr2e->SetMarkerStyle(20);
    gr2e->SetMarkerColor(kBlack);
    gr2e->SetMarkerSize(3);
    gr3n->SetMarkerStyle(5);
    gr3n->SetMarkerColor(kRed);
    gr3n->SetMarkerSize(3);
    gr3e->SetMarkerStyle(20);
    gr3e->SetMarkerColor(kRed);
    gr3e->SetMarkerSize(3);
    gr4n->SetMarkerStyle(5);
    gr4n->SetMarkerColor(kGreen);
    gr4n->SetMarkerSize(3);
    gr4e->SetMarkerStyle(20);
    gr4e->SetMarkerColor(kGreen);
    gr4e->SetMarkerSize(3);
    gr5n->SetMarkerStyle(5);
    gr5n->SetMarkerColor(kBlue);
    gr5n->SetMarkerSize(3);
    gr5e->SetMarkerStyle(20);
    gr5e->SetMarkerColor(kBlue);
    gr5e->SetMarkerSize(3);

    gr2nth->SetMarkerStyle(4);
    gr2nth->SetMarkerColor(kBlack);
    gr2nth->SetMarkerSize(3);
    gr2eth->SetMarkerStyle(24);
    gr2eth->SetMarkerColor(kBlack);
    gr2eth->SetMarkerSize(3);
    gr3nth->SetMarkerStyle(4);
    gr3nth->SetMarkerColor(kRed);
    gr3nth->SetMarkerSize(3);
    gr3eth->SetMarkerStyle(24);
    gr3eth->SetMarkerColor(kRed);
    gr3eth->SetMarkerSize(3);
    gr5nth->SetMarkerStyle(4);
    gr5nth->SetMarkerColor(kBlue);
    gr5nth->SetMarkerSize(3);
    gr5eth->SetMarkerStyle(24);
    gr5eth->SetMarkerColor(kBlue);
    gr5eth->SetMarkerSize(2);

    histo_cfg(h_axis, "Threshold (s/n)","N","");
    h_axis->SetMinimum(0.0);
    h_axis->SetMaximum(1.0);
    h_axis->Draw("hist");

    gr2n->Draw("P");
    gr2e->Draw("P");
    gr3n->Draw("P");
    gr3e->Draw("P");
    gr4n->Draw("P");
    gr4e->Draw("P");
    gr5n->Draw("P");
    gr5e->Draw("P");

    gr2nth->Draw("P");
    gr2eth->Draw("P");
    gr3nth->Draw("P");
    gr3eth->Draw("P");
    gr5nth->Draw("P");
    gr5eth->Draw("P");

    // The legend
    TLegend *leg = new TLegend(0.59,0.55,0.90,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Performance:");

    leg->AddEntry(gr2n,"2 GeV avg. noise","p");
    leg->AddEntry(gr2e,"2 GeV avg. efficiency","p");
    leg->AddEntry(gr3n,"3 GeV avg. noise","p");
    leg->AddEntry(gr3e,"3 GeV avg. efficiency","p");
    leg->AddEntry(gr4n,"4.4 GeV avg. noise","p");
    leg->AddEntry(gr4e,"4.4 GeV avg. efficiency","p");
    leg->AddEntry(gr5n,"5 GeV avg. noise","p");
    leg->AddEntry(gr5e,"5 GeV avg. efficiency","p");

    leg->AddEntry(gr2nth,"2 GeV avg. noise 20mm","p");
    leg->AddEntry(gr2eth,"2 GeV avg. efficiency 20mm","p");
    leg->AddEntry(gr3nth,"3 GeV avg. noise 20mm","p");
    leg->AddEntry(gr3eth,"3 GeV avg. efficiency 20mm","p");
    leg->AddEntry(gr5nth,"5 GeV avg. noise 20mm","p");
    leg->AddEntry(gr5eth,"5 GeV avg. efficiency 20mm","p");

    leg->Draw();

    // Output
    threshold->Print("pics/noiseeffi.eps");
    _outputFile->cd();
    threshold->Write();
    threshold->Close();

  }

  // Do a measured resolutions vs E plot
  if(runmode == 5)
  {

    cout << " " << endl;
    cout << "Mode 5" << endl;
    cout << " " << endl;

    Double_t resolution[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Double_t resolutionerror[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    Double_t energy[7] = { 2, 3, 4.4, 5, 2, 3, 5};
    Double_t energyerror[7] = { 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01};

    std::vector<int> energyruns;

    energyruns.push_back(5043);
    energyruns.push_back(5063);
    energyruns.push_back(98);
    energyruns.push_back(37);
    energyruns.push_back(236);
    energyruns.push_back(249);
    energyruns.push_back(262);

    for (int j=0;j<7;j++)
    {

      int runnumber = energyruns.at(j);

      TString filename("../filteredhistos/run00");
      if (runnumber <= 999)
	filename+="0";
      if (runnumber <= 99)
	filename+="0";
      filename+=runnumber;
      filename+="-fitter.root";
      TFile *f6 = new TFile(filename);
      //TH1D *h_m26[12];
      f6->cd();


      TH1D *h_m26[2];

      // Load file

      h_m26[0] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftX"))->Clone();
      h_m26[1] = (TH1D*) (f6->Get("DUTHisto03/DUTshiftY"))->Clone();

      float temp1 = 0.0;
      float temp2 = 0.0;

      for (int i=0;i<2;i++)
      {
	TF1 *f1 = new TF1("f1","gaus");
	h_m26[i]->Fit(f1,"EMQI0","");
	Double_t sigma_1 = f1->GetParameter(2);
	Double_t error_1 = f1->GetParError(2);

	temp1 += sigma_1;
	temp2 += error_1;

	cout << "mean 1 " << sigma_1 << endl;
      }

      temp1 = temp1 / 2.0;
      temp2 = temp2 / 2.0;

      resolution[j] = temp1*1000;
      resolutionerror[j] = temp2*1000;

      cout << "resolution is " << resolution[j] << " at j " << j << endl;

    }
    // Create graphs with the information
    TGraphErrors *gr = new TGraphErrors(7,energy,resolution,energyerror,resolutionerror);

    // Let's plot this
    TH1D *h_axis = new TH1D("th_axis","th_axis",1, 0.0, 10.0);
    TCanvas *energyplot = new TCanvas("energyplot","energyplot",10,10,800,600);
    gStyle->SetPadBorderMode(0);
    gStyle->SetOptStat(0);
    energyplot->SetFillColor(0);
    energyplot->Divide(1,1);
    energyplot->cd(1);
    gStyle->SetErrorX(0);
    gPad->SetLogx(0);
    gPad->SetLogy(0);

    histo_cfg(h_axis, "Threshold (s/n)","N","");
    h_axis->SetMinimum(0.0);
    h_axis->SetMaximum(10.0);
    h_axis->Draw("hist");

    gr->Draw("P");

    energyplot->Print("pics/energyplot.eps");
    _outputFile->cd();
    energyplot->Write();
    energyplot->Close();

  }


  // plot clustersize
  if(runmode == 6)
  {

    // Results go in here
    Double_t threshcluster2[12];
    Double_t threshcluster3[12];
    Double_t threshcluster4[12];
    Double_t threshcluster5[12];
    Double_t threshcluster120[12];
    Double_t threshcluster2_error[12];
    Double_t threshcluster3_error[12];
    Double_t threshcluster4_error[12];
    Double_t threshcluster5_error[12];
    Double_t threshcluster120_error[12];
    Double_t thin_threshcluster2[12];
    Double_t thin_threshcluster3[12];
    Double_t thin_threshcluster5[12];
    Double_t thin_threshcluster2_error[12];
    Double_t thin_threshcluster3_error[12];
    Double_t thin_threshcluster5_error[12];
    Double_t x[12];
    Double_t xerror[12] = {0.0};

    cout << " " << endl;
    cout << "Mode 6" << endl;
    cout << " " << endl;
    cout << "Running over all runs - cluster size" << endl;
    cout << " " << endl;
    cout << "Wide Geometry" << endl;
    cout << " " << endl;

    telescopebuild = "AnaTel.geom";
    planedistance = 150;
    for(int j=0;j<nplanes;j++)
      posx[j] = 150.0*j;

    for(int i=0;i<run2.size();i++)
    {
      getclusize( run2[i] );
      x[i] = i+3;
      threshcluster2[i] = avgclustersize;
      threshcluster2_error[i] = avgclustersize_error;
    }

    for(int i=0;i<run3.size();i++)
    {
      getclusize( run3[i] );
      x[i] = i+3;
      threshcluster3[i] = avgclustersize;
      threshcluster3_error[i] = avgclustersize_error;
    }

    for(int i=0;i<run4.size();i++)
    {
      getclusize( run4[i] );
      x[i] = i+3;
      threshcluster4[i] = avgclustersize;
      threshcluster4_error[i] = avgclustersize_error;
    }

    for(int i=0;i<run5.size();i++)
    {
      getclusize( run5[i] );
      x[i] = i+3;
      threshcluster5[i] = avgclustersize;
      threshcluster5_error[i] = avgclustersize_error;
    }

    for(int i=0;i<run120.size();i++)
    {
      getclusize( run120[i] );
      x[i] = i+3;
      threshcluster120[i] = avgclustersize;
      threshcluster120_error[i] = avgclustersize_error;
    }

    cout << " " << endl;
    cout << "Thin Geometry" << endl;
    cout << " " << endl;

    telescopebuild = "AnaTel_thin.geom";
    planedistance = 20;

    for(int j=0;j<nplanes;j++)
      posx[j] = 20.0*j;

    for(int i=0;i<run2th.size();i++)
    {
      getclusize( run2th[i] );
      x[i] = i+3;
      thin_threshcluster2[i] = avgclustersize;
      thin_threshcluster2_error[i] = avgclustersize_error;
    }

    for(int i=0;i<run3th.size();i++)
    {
      getclusize( run3th[i] );
      x[i] = i+3;
      thin_threshcluster3[i] = avgclustersize;
      thin_threshcluster3_error[i] = avgclustersize_error;
    }

    for(int i=0;i<run5th.size();i++)
    {
      getclusize( run5th[i] );
      x[i] = i+3;
      thin_threshcluster5[i] = avgclustersize;
      thin_threshcluster5_error[i] = avgclustersize_error;
    }

    // Create graphs with the information
    TGraphErrors *gr2n = new TGraphErrors(10,x,threshcluster2,xerror,threshcluster2_error);
    TGraphErrors *gr3n = new TGraphErrors(10,x,threshcluster3,xerror,threshcluster3_error);
    TGraphErrors *gr4n = new TGraphErrors(10,x,threshcluster4,xerror,threshcluster4_error);
    TGraphErrors *gr5n = new TGraphErrors(10,x,threshcluster5,xerror,threshcluster5_error);
    TGraphErrors *gr120n = new TGraphErrors(10,x,threshcluster120,xerror,threshcluster120_error);

    TGraphErrors *gr2nth = new TGraphErrors(10,x,thin_threshcluster2,xerror,thin_threshcluster2_error);
    TGraphErrors *gr3nth = new TGraphErrors(10,x,thin_threshcluster3,xerror,thin_threshcluster3_error);
    TGraphErrors *gr5nth = new TGraphErrors(10,x,thin_threshcluster5,xerror,thin_threshcluster5_error);

    // Let's plot this
    TH1D *h_axis = new TH1D("th_axis","th_axis",1, 2.0, 13.0);
    TCanvas *threshold = new TCanvas("threshold","threshold",10,10,800,600);
    gStyle->SetPadBorderMode(0);
    gStyle->SetOptStat(0);
    threshold->SetFillColor(0);
    threshold->Divide(1,1);
    threshold->cd(1);
    gStyle->SetErrorX(0);
    gPad->SetLogx(0);
    gPad->SetLogy(0);

    // Set apperance
    gr2n->SetMarkerStyle(22);
    gr2n->SetMarkerColor(kRed);
    gr2n->SetMarkerSize(2);
    gr3n->SetMarkerStyle(22);
    gr3n->SetMarkerColor(kBlue);
    gr3n->SetMarkerSize(2);
    gr4n->SetMarkerStyle(22);
    gr4n->SetMarkerColor(kGreen);
    gr4n->SetMarkerSize(2);
    gr5n->SetMarkerStyle(22);
    gr5n->SetMarkerColor(kBlack);
    gr5n->SetMarkerSize(2);
    gr120n->SetMarkerStyle(22);
    gr120n->SetMarkerColor(kOrange);
    gr120n->SetMarkerSize(2);

    gr2nth->SetMarkerStyle(34);
    gr2nth->SetMarkerColor(kRed);
    gr2nth->SetMarkerSize(2);
    gr3nth->SetMarkerStyle(34);
    gr3nth->SetMarkerColor(kBlue);
    gr3nth->SetMarkerSize(2);
    gr5nth->SetMarkerStyle(34);
    gr5nth->SetMarkerColor(kBlack);
    gr5nth->SetMarkerSize(2);

    histo_cfg(h_axis, "Threshold (s/n)","N","");
    h_axis->SetMinimum(0.0);
    h_axis->SetMaximum(3.0);
    h_axis->Draw("hist");

    gr2n->Draw("P");
    gr3n->Draw("P");
    gr4n->Draw("P");
    gr5n->Draw("P");
    gr120n->Draw("P");

    gr2nth->Draw("P");
    gr3nth->Draw("P");
    gr5nth->Draw("P");

    // The legend
    TLegend *leg = new TLegend(0.59,0.55,0.90,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    leg->SetHeader("Performance:");

    leg->AddEntry(gr2n,"2 GeV avg. cluster size","p");
    leg->AddEntry(gr3n,"3 GeV avg. cluster size","p");
    leg->AddEntry(gr4n,"4.4 GeV avg. cluster size","p");
    leg->AddEntry(gr5n,"5 GeV avg. cluster size","p");
    leg->AddEntry(gr120n,"120 GeV avg. cluster size","p");

    leg->AddEntry(gr2nth,"2 GeV avg. cluster size 20mm","p");
    leg->AddEntry(gr3nth,"3 GeV avg. cluster size 20mm","p");
    leg->AddEntry(gr5nth,"5 GeV avg. cluster size 20mm","p");

    leg->Draw();

    // Output
    threshold->Print("pics/clusize.eps");
    _outputFile->cd();
    threshold->Write();
    threshold->Close();

  }


  // plot pointing
  if(runmode == 7)
  {

    // Results go in here
    Double_t threshcluster2[12];
    Double_t threshcluster2_error[12];
    Double_t threshcluster2w[12];
    Double_t threshcluster2w_error[12];
    Double_t x[12];
    Double_t xerror[12] = {0.0};
    Double_t xw[12];
    Double_t xwerror[12] = {0.0};

    cout << " " << endl;
    cout << "Mode 7" << endl;
    cout << " " << endl;

    double thresh = 6-3;
    // x -3 with x the threshold

    getpointing( run2[thresh] , 3.415567);
    x[0] = 2;
    xerror[0] = 0.2;
    threshcluster2[0] = avgmeas;
    threshcluster2_error[0] = avgmeas_error;

    getpointing( run3[thresh] , 3.468273);
    x[1] = 3;
    xerror[1] = 0.3;
    threshcluster2[1] = avgmeas;
    threshcluster2_error[1] = avgmeas_error;

    getpointing( run4[thresh] , 3.530031);
    x[2] = 4.4;
    xerror[2] = 0.44;
    threshcluster2[2] = avgmeas;
    threshcluster2_error[2] = avgmeas_error;

    getpointing( run5[thresh] , 3.444235);
    x[3] = 5;
    xerror[3] = 0.5;
    threshcluster2[3] = avgmeas;
    threshcluster2_error[3] = avgmeas_error;


    // 3.194782 thr 6
    // 3.17949 thr 7
    // 3.31353 thr 8
    // 3.4976 thr 9 
    getpointing( run120[thresh+2] , 3.31353);
    x[4] = 120;
    xerror[4] = 12;
    threshcluster2[4] = avgmeas;
    threshcluster2_error[4] = avgmeas_error;

    getpointing( run2th[thresh] , 3.490106);
    xw[0] = 2;
    xwerror[0] = 0.2;
    threshcluster2w[0] = avgmeas;
    threshcluster2w_error[0] = avgmeas_error;

    getpointing( run3th[thresh] , 3.445178);
    xw[1] = 3;
    xwerror[1] = 0.3;
    threshcluster2w[1] = avgmeas;
    threshcluster2w_error[1] = avgmeas_error;

    getpointing( run5th[thresh] , 3.420506);
    xw[2] = 5;
    xwerror[2] = 0.5;
    threshcluster2w[2] = avgmeas;
    threshcluster2w_error[2] = avgmeas_error;


    //getpointing( testing[0] , 3.194782);

    //getpointing( 293 , 3.01298);

    //getpointing( 294 , 3.17294);


    //getpointing( 295 , 3.3148);
    //getpointing( 296 , 3.60078);

    // all fixed but m26

    getpointing( 293 , 2.51878);
    getpointing( 294 , 2.59982);
    getpointing( 295 , 2.74087);
    getpointing( 296 , 2.96215);
    getpointing( 297 , 3.17257);
    getpointing( 298 , 3.33224);
    getpointing( 299 , 3.51072);
    getpointing( 300 , 3.60078);

    // all fixed but m26, energy
    /*
       getpointing( 293 , 3.08702);
       getpointing( 294 , 3.17234);
       getpointing( 295 , 3.3148);
       getpointing( 296 , 3.52033);
       getpointing( 297 , 3.72182);
       getpointing( 298 , 3.88004);
       getpointing( 299 , 4.04436);
       getpointing( 300 , 4.14077);
     */

    // all fixed but m26, energy, thick
    /*
       getpointing( 294 , 3.16115);
       getpointing( 295 , 3.32987);
       getpointing( 296 , 3.52034);
       getpointing( 297 , 3.74096);
       getpointing( 298 , 3.88069);
       getpointing( 299 , 4.02916);
       getpointing( 300 , 4.1261);
     */
    x[5] = 12.5;
    xerror[5] = 1.25;
    threshcluster2[5] = avgmeas;
    threshcluster2_error[5] = avgmeas_error;



    // Create graphs with the information
    TGraphErrors *gr2n = new TGraphErrors(6,x,threshcluster2,xerror,threshcluster2_error);
    TGraphErrors *gr2w = new TGraphErrors(3,xw,threshcluster2w,xwerror,threshcluster2w_error);


    // Let's plot this
    TH1D *h_axis = new TH1D("th_axis","th_axis",100, 0.0, 130.0);
    TCanvas *threshold = new TCanvas("threshold","threshold",10,10,800,600);
    gStyle->SetPadBorderMode(0);
    gStyle->SetOptStat(0);
    threshold->SetFillColor(0);
    threshold->Divide(1,1);
    threshold->cd(1);
    gStyle->SetErrorX(0);
    gPad->SetLogx(0);
    gPad->SetLogy(0);

    // Set apperance
    gr2n->SetMarkerStyle(22);
    gr2n->SetMarkerColor(kRed);
    gr2n->SetMarkerSize(3);
    gr2w->SetMarkerStyle(22);
    gr2w->SetMarkerColor(kBlue);
    gr2w->SetMarkerSize(3);


    histo_cfg(h_axis, "Beam Energy in GeV","#sigma_{Point}","");
    h_axis->SetMinimum(0.0);
    h_axis->SetMaximum(10.0);
    h_axis->Draw("hist");

    gr2n->Draw("P");
    gr2w->Draw("P");


    // The legend
    TLegend *leg = new TLegend(0.59,0.55,0.90,0.85);
    leg->SetBorderSize(0);
    leg->SetFillColor(0);
    leg->SetFillStyle(0);
    // leg->SetHeader("Performance:");

    leg->AddEntry(gr2n,"#Delta_{z} = 150mm","p");
    leg->AddEntry(gr2w,"#Delta_{z} = 20mm","p");
    /*    leg->AddEntry(gr4n,"4.4 GeV avg. cluster size","p");
	  leg->AddEntry(gr5n,"5 GeV avg. cluster size","p");
	  leg->AddEntry(gr120n,"120 GeV avg. cluster size","p");

	  leg->AddEntry(gr2nth,"2 GeV avg. cluster size 20mm","p");
	  leg->AddEntry(gr3nth,"3 GeV avg. cluster size 20mm","p");
	  leg->AddEntry(gr5nth,"5 GeV avg. cluster size 20mm","p");
     */
    leg->Draw();

    // Output
    threshold->Print("pics/pointing.eps");
    _outputFile->cd();
    threshold->Write();
    threshold->Close();

  }

  if(runmode == 8)
  {
    cout << "test mode" << endl;

    // histoplot(37);

    for (int i = 0;i<10;i++)
    {
      histoplot (run2[i]);
      histoplot (run3[i]);
      histoplot (run4[i]);
      histoplot (run5[i]);
      histoplot (run120[i]);
      histoplot (run5th[i]);
      histoplot (run3th[i]);
      histoplot (run2th[i]);

    }


    _outputFile->cd();
    delta0->Write();
    delta1->Write();
    delta2->Write();
    delta3->Write();
    delta4->Write();
    delta5->Write();
    delta6->Write();
    delta7->Write();
    delta8->Write();
    delta9->Write();
    delta10->Write();
    delta11->Write();


  } 


  _outputFile->Close();

  // And we're done
  cout << endl;
  cout << "Done :-)" << endl;
  cout << endl;

}
