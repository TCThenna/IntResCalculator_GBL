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

/*   Exchange pointing resolution estimate by GBL calculation
     
     Author: H. Jansen
     email hendrik.jansen@desy.de
     Date: 3.8.15
*/

#include "GblTrajectory.h"
#include "GblPoint.h"
#include "GblData.h"


#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "AnaTel.h"
//#include "TF1.h"


//GblPoint getPoint(double dz, double res, TVectorD wscat, bool has_meas = true) {
//
//  // Propagate:
//  TMatrixD jacPointToPoint = Jac5(dz);
//  GblPoint point(jacPointToPoint);
//
//  // Add scatterer:
//  TVectorD scat(2);
//  scat.Zero(); // mean is zero
//  point.addScatterer(scat, wscat);
//
//  // Add measurement if requested:
//  if(has_meas) {
//    // measurement = residual
//    TVectorD meas(2);
//    meas.Zero(); // ideal
//
//    // Precision = 1/resolution^2
//    TVectorD measPrec(2);
//    measPrec[0] = 1.0 / res / res;
//    measPrec[1] = 1.0 / res / res;
//
//    TMatrixD proL2m(2,2);
//    proL2m.UnitMatrix();
// 
//    point.addMeasurement(proL2m, meas, measPrec);
//  }
//
//  return point;
//}
//
//GblPoint getPoint(double dz, TVectorD wscat) {
//  // This does not add a measurement - no resultion is given!
//  return getPoint(dz,0.0,wscat,false);
//}




void AnaTel::SetPlane(Int_t Ipl, Double_t Position,  Double_t Thickness,   Double_t X0, Double_t Resolution)
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


      std::stringstream ss ; 
      if(_isActive[Ipl])
	  ss << "Active  plane at" ;
      else
	  ss << "Passive plane at" ; 
      
      ss << "  Z [mm] = " << _planePosition[Ipl] 
	 << " dZ [um] = " << _planeThickness[Ipl]*1000. ;
      
      if(_isActive[Ipl])
	ss << "  Res [um] = " << _planeResolution[Ipl]*1000. ;
      
     std::cout << ss.str()  << std::endl;
} 


void  AnaTel::SetDUT(Int_t idut, Bool_t UseInFit)
{
    _iDUT=idut;
    _useDUT=UseInFit;

   if (_iDUT >=0 && _iDUT <  _nTelPlanes) 
     {
       if(_isActive[_iDUT] && UseInFit)
	  std::cout << "Active plane " << idut << "  set as DUT, " << "used in fit" << std::endl;
       else
	 if(_isActive[_iDUT])
	  std::cout << "Active plane " << idut << "  set as DUT, " << "not used in fit" << std::endl;
	 else
          std::cout  << "Passive plane " << idut << "  set as DUT "  << std::endl;

 
       if(!_isActive[_iDUT] && UseInFit)
         {
	 std::cout << "Warning: passive plane can not be used in fit "  << std::endl;
         _useDUT=false;
         }

     }
   else
	  std::cout << "Setup without DUT "  << std::endl;

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
  // std::cout <<  "Assuming beam energy of " << _eBeam 
  //      << " GeV and beam spread of " << _beamSpread << " rad" << std::endl ;
  // else
  // std::cout <<  "Assuming beam energy of " << _eBeam 
  //      << " GeV and no beam constraint" << std::endl ;

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



Double_t AnaTel::GetError(Int_t Ipl, Bool_t UseInFit)
{
  if(_eBeam <=0.)
    {
    std::cerr << " Beam energy not set ! Use SetBeam() method" << std::endl ;
    return 0.;
    }


  float oldscat = 1.;
  for(int ipl=0; ipl<_nTelPlanes ; ipl++)
    {
    if(ipl>0)
    {
      _planeDist[ipl-1]=1./(_planePosition[ipl] - _planePosition[ipl-1]) ;
      Double_t tempdist = (_planePosition[ipl] - _planePosition[ipl-1])*1.0;
      //std::cout << "tempdist is " << tempdist << std::endl;
      _planeScat[ipl]= 0.0136/_eBeam * sqrt(_planeThickness[ipl]/_planeX0[ipl] + tempdist/304200.0 + 0.05/286.0) * (1.+0.038*std::log(_planeThickness[ipl]/_planeX0[ipl] + tempdist/304200.0 + 0.05/286.0)) ;
      
    }
    if(ipl==0)
    {
      _planeScat[ipl]= 0.0136/_eBeam * sqrt(_planeThickness[ipl]/_planeX0[ipl] + 0.05/286.0) * (1.+0.038*std::log(_planeThickness[ipl]/_planeX0[ipl] + 0.05/286.0)) ;
    }
    if(ipl==0 && _useBeamConstraint)
      _planeScat[ipl]= 1./(_planeScat[ipl]*_planeScat[ipl]+ _beamSpread*_beamSpread) ; 
    else
      _planeScat[ipl]= 1./(_planeScat[ipl] * _planeScat[ipl]) ; 

    _fitX[ipl] = 0. ;
    _nominalError[ipl]= _planeResolution[ipl];
    }

  // Fit with nominal parameters

  if(!UseInFit) {
    _nominalError[Ipl]=0;
    _fitX[Ipl] = 0. ;
  }
  
  int status = DoAnalFit(_fitX,_nominalError);
  
  //std::cout << "After fit for plane " << Ipl << " is " << _nominalError[Ipl] << std::endl;

  if(status)
    std::cerr << "\n Fit failed !?!" << std::endl ;

 
  return _nominalError[Ipl];
}


Double_t AnaTel::GetWidth(Int_t Ipl, Bool_t UseInFit)
{
  
  

  if(_planeResolution[Ipl] <=0.)
    {
    std::cerr << " Can not calculate residual width for passive plane !" << std::endl ;
    return 0.;
    }


  Double_t FitError=GetError( Ipl, UseInFit);

  Double_t Width;

  if(UseInFit)
      Width=sqrt(_planeResolution[Ipl]*_planeResolution[Ipl]-FitError*FitError);
  else
      Width=sqrt(_planeResolution[Ipl]*_planeResolution[Ipl]+FitError*FitError);

  //if(Ipl == 0) std::cout << "sigma_i " << Ipl << " = " << _planeResolution[Ipl]*1000.0 << "  pointing res = " << FitError*1000.0 << std::endl;
  //if(Ipl == 3) std::cout << "sigma_i " << Ipl << " = " << _planeResolution[Ipl]*1000.0 << "  pointing res = " << FitError*1000.0 << std::endl;
  //std::cout << " " << std::endl;
  
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
      //std::cout << pos[ipl] << " ";
      if(_isActive[ipl] && err[ipl]>0)
        err[ipl]=1./err[ipl]/err[ipl] ;
      else
        err[ipl] = 0. ;

      pos[ipl]*=err[ipl];
    }
    //std::cout << std::endl;

  for(int ipl=0; ipl<_nTelPlanes;ipl++)
    {
      //std::cout << pos[ipl] << " ";
    }
    //std::cout << std::endl;

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
      std::cerr << "Singular matrix in track fitting algorithm ! " << std::endl;
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
      //std::cout << " beta[icol] = " << beta[icol] << std::endl;

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



