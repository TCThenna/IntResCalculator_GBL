/*    ROOT class for analytical EUDET telescope error estimate 
 *
 *    Author: A.F.Zarnecki, University of Warsaw 
 *    email: zarnecki@fuw.edu.pl
 *    Date: April 2008
 *
 *   This source code was developed with the EUDET project and is
 *   based on the algorithm implemented in the Eutelescope package of Marlin.
 *   You are free to use this source file for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 *   For results based on this code please refer to:  
 *      A.F.Zarnecki, P.Niezurawski, 
 *      'EUDET Telescope Geometry and Resolution Studies',
 *      EUDET-REPORT-2007-01, Mar 2007. 
 */


// Edited July 2013 Thomas Eichhorn

#ifndef AnaTel_Included
#define AnaTel_Included

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
//#include "TF1.h"

class AnaTel {

public:
  
  virtual ~AnaTel();

  AnaTel(char * GeomFile);
  // Constructor: create telescope from geometry file

  AnaTel(Int_t Npl);
  // Constructor: create telescope "by hand"

  void SetPlane(Int_t Ipl, Double_t Position,  Double_t Thickness,  
	   Double_t X0, Double_t Resolution); 

  void  SetResolution(Int_t Ipl, Double_t Resolution);

  void  SetResolution(Double_t Resolution);

  void  SetResolution(Double_t * ResArray);

  void  SetThickness(Double_t Thickness);

  void SetBeam(Double_t Energy, Double_t Spread=0.);

  void SetDUT(Int_t idut, Bool_t UseInFit);

  Int_t GetNplanes();

  Double_t * GetPosition();

  Double_t * GetThickness();

  Double_t * GetResolution();

  void GetPlane(Int_t Ipl, Double_t  * Position, Double_t  * Thickness,  
	   Double_t  * X0, Double_t  * Resolution);

  Double_t GetError(Int_t Ipl, Bool_t UseInFit);

  Double_t GetWidth(Int_t Ipl, Bool_t UseInFit);

  Double_t GetDUTError();

  Double_t GetDUTWidth();
  
private:

    int _nTelPlanes;
    double _eBeam ;
    int _iDUT;
    bool   _useDUT;

    double * _planePosition;
    double * _planeThickness;
    double * _planeX0;
    double * _planeResolution;
    bool   * _isActive;

    double * _planeDist ;
    double * _planeScat ;
    double * _tempScat1 ;
    double * _tempScat2 ;
    double * _tempScat3 ;

    double * _fitX  ;
    double * _fitArray ;
    double * _nominalError ;

    bool   _useBeamConstraint ;
    double _beamSpread;

    int DoAnalFit(double * , double *);
    int GaussjSolve(double * ,double * ,int );


  ClassDef(AnaTel,1);

};


AnaTel::AnaTel(char * GeomFile)
{
  ifstream geometryFile;

  geometryFile.open(GeomFile);

  //cout << "Reading telescope geometry description from " << GeomFile << endl;

  geometryFile >> _nTelPlanes  >> _iDUT ;

  // !!! Important: DUT position in geometry file numbered 0...N-1 !!!
  //     _iDUT = -1 denotes no DUT in the Setup

  _planePosition   = new double[_nTelPlanes];
  _planeThickness  = new double[_nTelPlanes];
  _planeX0         = new double[_nTelPlanes];
  _planeResolution = new double[_nTelPlanes];
  _isActive        = new bool[_nTelPlanes];

  // Planes have to be ordered in position along the beam line !
  // This is not checked !!!

  for(int ipl=0; ipl < _nTelPlanes; ipl++)
    {
      int iActive;
      double resolution;

      // All dimensions should be given in mm !!!

      geometryFile >> _planePosition[ipl]
                   >> _planeThickness[ipl]
                   >> _planeX0[ipl]
                   >> iActive 
                   >> resolution ;

      if(iActive)
	{
	  _isActive[ipl] = true ;
	  _planeResolution[ipl]=resolution;
	}
      else
	{
          _isActive[ipl] = false ;
          _planeResolution[ipl]=0.;
	}
    }

  // Print out geometry information
  //cout << "Telescope configuration with " << _nTelPlanes << " planes" << endl;


  for(int ipl=0; ipl < _nTelPlanes; ipl++)
    {
      stringstream ss ; 

      if ( ipl == _iDUT)
	{
        if(_isActive[ipl])
	  ss << ipl << " : active  DUT   at" ;
        else
	  ss << ipl << " : passive  DUT  at" ; 
	}
      else
	{
        if(_isActive[ipl])
	  ss << ipl << " : active  plane at" ;
        else
	  ss << ipl << " : passive plane at" ; 
	}

      
      ss << "  Z [mm] = " << _planePosition[ipl] 
	 << " dZ [um] = " << _planeThickness[ipl]*1000. ;
      
      if(_isActive[ipl])
	ss << "  Res [um] = " << _planeResolution[ipl]*1000. ;
      
     //cout << ss.str()  << endl;
    }
  // Allocate arrays for track fitting

  _planeDist = new double[_nTelPlanes];
  _planeScat = new double[_nTelPlanes];
  _tempScat1 = new double[_nTelPlanes];
  _tempScat2 = new double[_nTelPlanes];
  _tempScat3 = new double[_nTelPlanes];

  _fitX  = new double[_nTelPlanes];
  _nominalError  = new double[_nTelPlanes];

  int arrayDim = _nTelPlanes * _nTelPlanes;

  _fitArray = new double[arrayDim];


  if(_iDUT>=0)
    _useDUT=true;
  else
    _useDUT=false;

  _eBeam=0.;
  _useBeamConstraint=false ;
}


AnaTel::AnaTel(Int_t Npl)
{
  _nTelPlanes=Npl;

  _planePosition   = new double[_nTelPlanes];
  _planeThickness  = new double[_nTelPlanes];
  _planeX0         = new double[_nTelPlanes];
  _planeResolution = new double[_nTelPlanes];
  _isActive        = new bool[_nTelPlanes];

  _planeDist = new double[_nTelPlanes];
  _planeScat = new double[_nTelPlanes];
  _tempScat1 = new double[_nTelPlanes];
  _tempScat2 = new double[_nTelPlanes];
  _tempScat3 = new double[_nTelPlanes];

  _fitX  = new double[_nTelPlanes];
  _nominalError  = new double[_nTelPlanes];

  int arrayDim = _nTelPlanes * _nTelPlanes;

  _fitArray = new double[arrayDim];
  
  _iDUT=-1;
  _useDUT=false;

  _useBeamConstraint=false ;
}


void AnaTel::SetPlane(Int_t Ipl, Double_t Position,  Double_t Thickness,  
	   Double_t X0, Double_t Resolution)
{

  _planePosition[Ipl]=Position;
  _planeThickness[Ipl]=Thickness;
  _planeX0[Ipl]=X0;
 
  if(Resolution>0)
	{
	  _isActive[Ipl] = true ;
	  _planeResolution[Ipl]=Resolution;
	}
    else
	{
          _isActive[Ipl] = false ;
          _planeResolution[Ipl]=0.;
	}


      stringstream ss ; 
      if(_isActive[Ipl])
	  ss << "Active  plane at" ;
      else
	  ss << "Passive plane at" ; 
      
      ss << "  Z [mm] = " << _planePosition[Ipl] 
	 << " dZ [um] = " << _planeThickness[Ipl]*1000. ;
      
      if(_isActive[Ipl])
	ss << "  Res [um] = " << _planeResolution[Ipl]*1000. ;
      
     cout << ss.str()  << endl;
} 


void  AnaTel::SetDUT(Int_t idut, Bool_t UseInFit)
{
    _iDUT=idut;
    _useDUT=UseInFit;

   if (_iDUT >=0 && _iDUT <  _nTelPlanes) 
     {
       if(_isActive[_iDUT] && UseInFit)
	  cout << "Active plane " << idut << "  set as DUT, " << "used in fit" << endl;
       else
	 if(_isActive[_iDUT])
	  cout << "Active plane " << idut << "  set as DUT, " << "not used in fit" << endl;
	 else
          cout  << "Passive plane " << idut << "  set as DUT "  << endl;

 
       if(!_isActive[_iDUT] && UseInFit)
         {
	 cout << "Warning: passive plane can not be used in fit "  << endl;
         _useDUT=false;
         }

     }
   else
	  cout << "Setup without DUT "  << endl;

}

void  AnaTel::SetBeam(Double_t Energy, Double_t Spread)
{
  _eBeam = Energy;

  if(Spread>0.)
    {
     _useBeamConstraint= true ;
     _beamSpread = Spread;
    }
  else
    {
     _useBeamConstraint= false ;
     _beamSpread = 0.;
    }

  //   if(_useBeamConstraint)
  // cout <<  "Assuming beam energy of " << _eBeam 
  //      << " GeV and beam spread of " << _beamSpread << " rad" << endl ;
  // else
  // cout <<  "Assuming beam energy of " << _eBeam 
  //      << " GeV and no beam constraint" << endl ;

}

void AnaTel::SetResolution(Int_t Ipl, Double_t Resolution)
{
_planeResolution[Ipl]= Resolution;
}


void AnaTel::SetThickness(Double_t thickness)
{
  for(int ipl=0;ipl<_nTelPlanes;ipl++)
       _planeThickness[ipl]= thickness;
}


void AnaTel::SetResolution(Double_t Resolution)
{
  for(int ipl=0;ipl<_nTelPlanes;ipl++)
       _planeResolution[ipl]= Resolution;
}


void AnaTel::SetResolution(Double_t * Resolution)
{
  for(int ipl=0;ipl<_nTelPlanes;ipl++)
       _planeResolution[ipl]= Resolution[ipl];
}


Int_t AnaTel::GetNplanes()
{
return _nTelPlanes;
}

Double_t * AnaTel::GetPosition()
{
return _planePosition;
}

Double_t * AnaTel::GetThickness()
{
return _planeThickness;
}

Double_t * AnaTel::GetResolution()
{
return _planeResolution;
}

void AnaTel::GetPlane(Int_t Ipl, Double_t * Position,  Double_t * Thickness,  
	   Double_t * X0, Double_t * Resolution)
{
  *Position=_planePosition[Ipl];
  *Thickness=_planeThickness[Ipl];
  *X0=_planeX0[Ipl];
  *Resolution=_planeResolution[Ipl];
} 

/*
Double_t AnaTel::GetError(Int_t Ipl, Bool_t UseInFit)
{
  if(_eBeam <=0.)
    {
    cerr << " Beam energy not set ! Use SetBeam() method" << endl ;
    return 0.;
    }

  for(int ipl=0; ipl<_nTelPlanes ; ipl++)
    {
      if(ipl>0)
	_planeDist[ipl-1]=1./(_planePosition[ipl] - _planePosition[ipl-1]) ;

     _tempScat1[ipl]= (0.0136/_eBeam * sqrt(_planeThickness[ipl]/_planeX0[ipl])*(1.+0.038*std::log(_planeThickness[ipl]/_planeX0[ipl])));// + (0.0136/_eBeam * sqrt(150.0/304200.0) * (1.+0.038*std::log(150.0/304200.0)));
      _tempScat2[ipl] = (0.0136/_eBeam * sqrt(150.0/304200.0) * (1.+0.038*std::log(150.0/304200.0)));
	_tempScat3[ipl] = (0.0136/_eBeam * sqrt(0.005/286.0) * (1.+0.038*std::log(0.005/286.0))) ;
      //_tempScat3[ipl] = 1.0;


      if(ipl==0 && _useBeamConstraint)
      {
	_planeScat[ipl]= 1./(_tempScat1[ipl]*_tempScat1[ipl]+ _beamSpread*_beamSpread);
//	_planeScat[ipl]+= 1./(_tempScat2[ipl]*_tempScat2[ipl]+ _beamSpread*_beamSpread);
//	_planeScat[ipl]+= 1./(_tempScat3[ipl]*_tempScat3[ipl]+ _beamSpread*_beamSpread);
      }
      else
      {
	_planeScat[ipl]= 1./(_tempScat1[ipl] * _tempScat1[ipl]);
	_planeScat[ipl]+= 1./(_tempScat2[ipl] * _tempScat2[ipl]);
	_planeScat[ipl]+= 1./(_tempScat3[ipl] * _tempScat3[ipl]);
      }
      
//      if(ipl==0 || ipl==_nTelPlanes)
//	_planeScat[ipl]= 1./(_tempScat1[ipl] * _tempScat1[ipl]);	


      _fitX[ipl] = 0. ;
      _nominalError[ipl]= _planeResolution[ipl];
    }

  // Fit with nominal parameters

  if(!UseInFit) _nominalError[Ipl]=0;

  int status = DoAnalFit(_fitX,_nominalError);

  if(status)
    cerr << "\n Fit failed !?!" << endl ;

  return _nominalError[Ipl];
}*/

Double_t AnaTel::GetError(Int_t Ipl, Bool_t UseInFit)
{
  if(_eBeam <=0.)
    {
    cerr << " Beam energy not set ! Use SetBeam() method" << endl ;
    return 0.;
    }


  float oldscat = 1.;
  for(int ipl=0; ipl<_nTelPlanes ; ipl++)
    {
    if(ipl>0)
      _planeDist[ipl-1]=1./(_planePosition[ipl] - _planePosition[ipl-1]) ;

    //_planeScat[ipl]= 0.0136/_eBeam * sqrt(_planeThickness[ipl]/_planeX0[ipl] + 150.0/304200.0 + 0.1/286.0) * (1.+0.038*std::log(_planeThickness[ipl]/_planeX0[ipl] + 150.0/304200.0 + 0.05/286.0)) ;

    if(ipl>0)
    {
      Double_t tempdist = (_planePosition[ipl] - _planePosition[ipl-1])*1.0;
      //cout << "tempdist is " << tempdist << endl;
      //  Double_t tempdist = 0;
      _planeScat[ipl]= 0.0136/_eBeam * sqrt(_planeThickness[ipl]/_planeX0[ipl] + tempdist/304200.0 + 0.05/286.0) * (1.+0.038*std::log(_planeThickness[ipl]/_planeX0[ipl] + tempdist/304200.0 + 0.05/286.0)) ;
      
      //if ( abs((_planeScat[ipl] - oldscat)/oldscat) > 0.01) { 
      //  cout << "scatter is " << ipl << " : " << _planeScat[ipl] << "   thick = " << _planeThickness[ipl] << endl;
//	oldscat = _planeScat[ipl];
      //}
    }
    if(ipl==0)
    {
      _planeScat[ipl]= 0.0136/_eBeam * sqrt(_planeThickness[ipl]/_planeX0[ipl] + 0.05/286.0) * (1.+0.038*std::log(_planeThickness[ipl]/_planeX0[ipl] + 0.05/286.0)) ;
      // why 0.1 ? is this a wrong estimate of the capton thickness??
    }
    if(ipl==0 && _useBeamConstraint)
      _planeScat[ipl]= 1./(_planeScat[ipl]*_planeScat[ipl]+ _beamSpread*_beamSpread) ; 
    else
      _planeScat[ipl]= 1./(_planeScat[ipl] * _planeScat[ipl]) ; 

    _fitX[ipl] = 1. ;
    _nominalError[ipl]= _planeResolution[ipl];
    }

  // Fit with nominal parameters

  if(!UseInFit) {
    _nominalError[Ipl]=0;
    _fitX[Ipl] = 0. ;
  }
  
  int status = DoAnalFit(_fitX,_nominalError);
  
  //cout << "After fit for plane " << Ipl << " is " << _nominalError[Ipl] << endl;

  if(status)
    cerr << "\n Fit failed !?!" << endl ;

 
  return _nominalError[Ipl];
}


Double_t AnaTel::GetWidth(Int_t Ipl, Bool_t UseInFit)
{
  
  

  if(_planeResolution[Ipl] <=0.)
    {
    cerr << " Can not calculate residual width for passive plane !" << endl ;
    return 0.;
    }


  Double_t FitError=GetError( Ipl, UseInFit);

  Double_t Width;

  if(UseInFit)
      Width=sqrt(_planeResolution[Ipl]*_planeResolution[Ipl]-FitError*FitError);
  else
      Width=sqrt(_planeResolution[Ipl]*_planeResolution[Ipl]+FitError*FitError);

  //if(Ipl == 0) cout << "sigma_i " << Ipl << " = " << _planeResolution[Ipl]*1000.0 << "  pointing res = " << FitError*1000.0 << endl;
  //if(Ipl == 3) cout << "sigma_i " << Ipl << " = " << _planeResolution[Ipl]*1000.0 << "  pointing res = " << FitError*1000.0 << endl;
  //cout << " " << endl;
  
  return Width;
}


Double_t AnaTel::GetDUTError()
{
  Double_t width=0.;

  if (_iDUT >=0 && _iDUT <  _nTelPlanes) width = GetError(_iDUT,_useDUT);

  return width;
}


Double_t AnaTel::GetDUTWidth()
{
  Double_t width=0.;

  if (_iDUT >=0 && _iDUT <  _nTelPlanes) width = GetWidth(_iDUT,_useDUT);

  return width;
}


// ==============================================================
//
// Private functions
//

int AnaTel::DoAnalFit(double * pos, double *err)
{
  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      //cout << pos[ipl] << " ";
      if(_isActive[ipl] && err[ipl]>0)
        err[ipl]=1./err[ipl]/err[ipl] ;
      else
        err[ipl] = 0. ;

      pos[ipl]*=err[ipl];
    }
    //cout << endl;

  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      //cout << pos[ipl] << " ";
    }
    //cout << endl;

  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    for(int jpl=0; jpl<_nTelPlanes;jpl++)
      {
	int imx=ipl+jpl*_nTelPlanes;

	_fitArray[imx] = 0.;

	if(jpl==ipl-2)
	  _fitArray[imx] += _planeDist[ipl-2]*_planeDist[ipl-1]*_planeScat[ipl-1] ;

	if(jpl==ipl+2)
	  _fitArray[imx] += _planeDist[ipl]*_planeDist[ipl+1]*_planeScat[ipl+1] ;

	if(jpl==ipl-1)
	  {
	    if(ipl>0 &&  ipl < _nTelPlanes-1)
              _fitArray[imx] -= _planeDist[ipl-1]*(_planeDist[ipl]+_planeDist[ipl-1])*_planeScat[ipl] ;
	    if(ipl>1)
              _fitArray[imx] -= _planeDist[ipl-1]*(_planeDist[ipl-1]+_planeDist[ipl-2])*_planeScat[ipl-1] ;
	  }

	if(jpl==ipl+1)
	  {
	    if(ipl>0 && ipl < _nTelPlanes-1)
	      _fitArray[imx] -= _planeDist[ipl]*(_planeDist[ipl]+_planeDist[ipl-1])*_planeScat[ipl] ;
	    if(ipl < _nTelPlanes-2)
	      _fitArray[imx] -= _planeDist[ipl]*(_planeDist[ipl+1]+_planeDist[ipl])*_planeScat[ipl+1] ;
	  }

	if(jpl==ipl)
	  {
            _fitArray[imx] += err[ipl] ;

            if(ipl>0 && ipl<_nTelPlanes-1)
	      _fitArray[imx] += _planeScat[ipl]*(_planeDist[ipl]+_planeDist[ipl-1])*(_planeDist[ipl]+_planeDist[ipl-1]) ;

            if(ipl > 1 )
	      _fitArray[imx] += _planeScat[ipl-1]*_planeDist[ipl-1]*_planeDist[ipl-1] ;

            if(ipl < _nTelPlanes-2)
	      _fitArray[imx] += _planeScat[ipl+1]*_planeDist[ipl]*_planeDist[ipl] ;
	  }

	// For beam constraint

	if(ipl==jpl && ipl<2 && _useBeamConstraint)
	  _fitArray[imx] += _planeScat[0]*_planeDist[0]*_planeDist[0] ;

	if(ipl+jpl==1 && _useBeamConstraint) 
	  _fitArray[imx] -= _planeScat[0]*_planeDist[0]*_planeDist[0] ;


      }

  int status=GaussjSolve(_fitArray,pos,_nTelPlanes) ;

  if(status)
    {
      cerr << "Singular matrix in track fitting algorithm ! " << endl;
      for(int ipl=0;ipl<_nTelPlanes;ipl++)
	err[ipl]=0. ;
    }
  else
    for(int ipl=0;ipl<_nTelPlanes;ipl++)
      err[ipl]=sqrt(_fitArray[ipl+ipl*_nTelPlanes]);

  return status ;
}




int AnaTel::GaussjSolve(double *alfa,double *beta,int n)
{
  int *ipiv;
  int *indxr;
  int *indxc;
  int i,j,k;
  int irow=0;
  int icol=0;
  double abs,big,help,pivinv;

  ipiv = new int[n];
  indxr = new int[n];
  indxc = new int[n];

  for(i=0;i<n;ipiv[i++]=0);

  for(i=0;i<n;i++)
    {
      big=0.;
      for(j=0;j<n;j++)
	{
	  if(ipiv[j]==1)continue;
	  for(k=0;k<n;k++)
	    {
	      if(ipiv[k]!=0)continue;
	      abs=fabs(alfa[n*j+k]);
	      if(abs>big)
		{
		  big=abs;
		  irow=j;
		  icol=k;
		}
	    }
	}
      ipiv[icol]++;

      if(ipiv[icol]>1)
	return 1;

      if(irow!=icol)
	{
	  help=beta[irow];
	  beta[irow]=beta[icol];
	  beta[icol]=help;
	  for(j=0;j<n;j++)
	    {
	      help=alfa[n*irow+j];
	      alfa[n*irow+j]=alfa[n*icol+j];
	      alfa[n*icol+j]=help;
	    }
	}
      indxr[i]=irow;
      indxc[i]=icol;

      if(alfa[n*icol+icol]==0.)
	return 1;

      help=alfa[n*icol+icol];
      pivinv=1./help;
      alfa[n*icol+icol]=1.;
      for(j=0;j<n;alfa[n*icol+(j++)]*=pivinv);
      beta[icol]*=pivinv;
      //cout << " beta[icol] = " << beta[icol] << endl;

      for(j=0;j<n;j++)
	{
	  if(j==icol)continue;
	  help=alfa[n*j+icol];
	  alfa[n*j+icol]=0.;
	  for(k=0;k<n;k++)
	    alfa[n*j+k]-=alfa[n*icol+k]*help;
	  beta[j]-=beta[icol]*help;
	}
    }

  for(i=n-1;i>=0;i--)
    {
      if(indxr[i]==indxc[i])continue;
      for(j=0;j<n;j++)
	{
	  help=alfa[n*j+indxr[i]];
	  alfa[n*j+indxr[i]]=alfa[n*j+indxc[i]];
	  alfa[n*j+indxc[i]]=help;
	}
    }

  delete [] ipiv;
  delete [] indxr;
  delete [] indxc;

  return 0;
}


#endif

