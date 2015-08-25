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
#include "TMatrixD.h"
//#include "TF1.h"

TMatrixD Jac5( double ds )
{
  /*
     straight line, no B-field
     track = 
     q/p, x', y', x, y
     0,   1,  2,  3, 4
     */
  TMatrixD jac(5, 5);
  jac.UnitMatrix();
  jac[3][1] = ds; // x = xp * ds
  jac[4][2] = ds; // y = yp * ds
  return jac;
}

gbl::GblPoint AnaTel::getPoint(double step, double res, TVectorD wscat, bool IsPlane, bool has_meas) {

  // Propagate:
  TMatrixD jacPointToPoint = Jac5(step);
  gbl::GblPoint point(jacPointToPoint);

  // Add scatterer:
  TVectorD scat(2);
  scat.Zero(); // mean is zero
  point.addScatterer(scat, wscat);

  // Add measurement if requested:
  if(has_meas) {
    // measurement = residual
    TVectorD meas(2);
    meas.Zero(); // ideal

    // Precision = 1/resolution^2
    TVectorD measPrec(2);
    measPrec[0] = 1.0 / res / res;
    measPrec[1] = 1.0 / res / res;

    TMatrixD proL2m(2,2);
    proL2m.UnitMatrix();

    point.addMeasurement(proL2m, meas, measPrec);

  }

  if(IsPlane) _ID.push_back(_sPoint.size());

  _s += step;
  _sPoint.push_back(_s);

  return point;
}

gbl::GblPoint AnaTel::getPoint(double step, TVectorD wscat, bool IsPlane) {
  // This does not add a measurement - no resultion is given!
  return AnaTel::getPoint(step, 0.0, wscat, IsPlane, false);
}

AnaTel::AnaTel(std::string GeomFile, double _eBeam)
{

  AnaTel::SetBeam(_eBeam, 0.0);
  ifstream geometryFile;

  std::string help = "../geometries/";
  help += GeomFile;
  std::cout << " Reading geometry file from: " << help << std::endl;
  geometryFile.open(help.c_str());

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

  std::cout << " nTelPlanes = " << _nTelPlanes << std::endl;

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
    std::stringstream ss ; 

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
  _pointingResolution  = new double[_nTelPlanes];

  int arrayDim = _nTelPlanes * _nTelPlanes;

  _fitArray = new double[arrayDim];


  if(_iDUT>=0)
    _useDUT=true;
  else
    _useDUT=false;

  //_eBeam=0.; // not needed anymore as constructor requires to pass the energy already
  _useBeamConstraint=false ;

  // GBL

  _s = 0;
  _sPoint.clear();
  _listOfPoints.clear();

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
  _pointingResolution  = new double[_nTelPlanes];

  int arrayDim = _nTelPlanes * _nTelPlanes;

  _fitArray = new double[arrayDim];

  _iDUT=-1;
  _useDUT=false;

  _useBeamConstraint=false ;
}




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

}

void AnaTel::SetResolution(Int_t Ipl, Double_t res)
{
  _planeResolution[Ipl] = res;
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

void AnaTel::SetThickness(Double_t thickness)
{
  for(int ipl=0;ipl<_nTelPlanes;ipl++)
    _planeThickness[ipl]= thickness;
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

Double_t AnaTel::GetPointingRes(Int_t Ipl, Bool_t UseInFit)
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
    _pointingResolution[ipl]= _planeResolution[ipl]; // first set to the intrinsic resolution, overwritten by DoAnalFit
  }

  // Fit with nominal parameters

  if(!UseInFit) {
    _pointingResolution[Ipl]=0;
    _fitX[Ipl] = 0. ;
  }

  int status = DoAnalFit(_fitX,_pointingResolution);

  //std::cout << "After fit for plane " << Ipl << " is " << _pointingResolution[Ipl] << std::endl;

  if(status)
    std::cerr << "\n Fit failed !?!" << std::endl ;


  return _pointingResolution[Ipl];
}

Double_t AnaTel::GetPointingResGBL(Int_t Ipl, Bool_t UseInFit)
{

  std::vector<double> _pointingResolutionLoc;
  _s = 0;
  _sPoint.clear();
  _listOfPoints.clear();

  std::cout << "Building the GBL t'scope " << std::endl;
  _listOfPoints.reserve(3*6);

  unsigned int ipl = 0;
  double step = 0;
  _s = 0;
  double tetSi, tetDUT, tetAir;

  TVectorD wscatSi(2); // FIXME can these be in the header ?
  TVectorD wscatDUT(2);
  TVectorD wscatAir(2);

  std::cout << "E = " << _eBeam << std::endl;
  while(1)
  {
    _X0_Si_frac = _planeThickness[ipl] / _X0_Si + 50.e-3 / _X0_Kapton; // Si + Kapton per plane
    //if(ipl == 0) std::cout << "X0 fraction of sensors = " << _X0_Si_frac << std::endl; // FIXME add verbose level

    tetSi  = 0.0136 * sqrt(_X0_Si_frac) / _eBeam * ( 1 + 0.038*log(_X0_Si_frac) );
    //tetDUT = 0.0136 * sqrt(_X0_DUT_frac) / _eBeam * ( 1 + 0.038*log(_X0_DUT_frac) ); // FIXME for now w/o DUT

    wscatSi[0] = 1.0 / ( tetSi * tetSi ); // weight
    wscatSi[1] = 1.0 / ( tetSi * tetSi );

    //wscatDUT[0] = 1.0 / ( tetDUT * tetDUT ); // weight // FIXME for now w/o DUT
    //wscatDUT[1] = 1.0 / ( tetDUT * tetDUT );



    if(UseInFit || ipl != Ipl) {
      _listOfPoints.push_back(AnaTel::getPoint(step, _planeResolution[ipl], wscatSi, true)); // add plane
    }
    else  _listOfPoints.push_back(AnaTel::getPoint(step, wscatSi, true)); // add plane w/o meas


    if((ipl+1) == _nTelPlanes) break;
    _dz = _planePosition[ipl+1] - _planePosition[ipl];
    // Air between planes; Factor 0.5 as the air is devided into two scatterers
    _X0_Air_frac =   0.5*_dz  / _X0_Air; 
    tetAir = 0.0136 * sqrt(_X0_Air_frac) / _eBeam * ( 1 + 0.038*log(_X0_Air_frac) );

    wscatAir[0] = 1.0 / ( tetAir * tetAir ); // weight
    wscatAir[1] = 1.0 / ( tetAir * tetAir );


    step = 0.21*_dz;
    _listOfPoints.push_back(AnaTel::getPoint(step, wscatAir, false)); // add air, not added after last plane
    step = 0.58*_dz;
    _listOfPoints.push_back(AnaTel::getPoint(step, wscatAir, false)); // add air, not added after last plane
    // Set step for distance to next plane
    step = 0.21*_dz; 

    ipl++;

  }

  // fit trajectory:

  gbl::GblTrajectory traj( _listOfPoints, 0 );
  //std::cout << " Is traj valid? " << traj.isValid() << std::endl;
  //traj.printPoints();
  //traj.printTrajectory();
  //traj.printData();

  double Chi2;
  int Ndf;
  double lostWeight;

  traj.fit( Chi2, Ndf, lostWeight );

  TVectorD aCorrection(5);
  TMatrixDSym aCovariance(5);

  for(int iipl = 0; iipl < _nTelPlanes; iipl++){
    traj.getResults( _ID[iipl], aCorrection, aCovariance );
    _pointingResolutionLoc.push_back(sqrt(aCovariance(3,3)));
    //std::cout << iipl << "  " << _ID[iipl] <<" pointingRes = " << _pointingResolutionLoc.back() << ",  ";
  }
  std::cout << std::endl;

  return _pointingResolutionLoc.at(Ipl);
}

Double_t AnaTel::GetWidth(Int_t Ipl, Bool_t UseInFit)
{

  if(_planeResolution[Ipl] <=0.)
  {
    std::cerr << " Can not calculate residual width for passive plane !" << std::endl ;
    return 0.;
  }

  Double_t PointingRes = GetPointingRes( Ipl, UseInFit);

  Double_t Width;

  //std::cout << " standard PR for plane " << Ipl << " = " << PointingRes << std::endl;
  if(UseInFit)
    Width=sqrt(_planeResolution[Ipl]*_planeResolution[Ipl]-PointingRes*PointingRes);
  else
    Width=sqrt(_planeResolution[Ipl]*_planeResolution[Ipl]+PointingRes*PointingRes);

  return Width;
}

Double_t AnaTel::GetWidthGBL(Int_t Ipl, Bool_t UseInFit)
{

  if(_planeResolution[Ipl] <=0.)
  {
    std::cerr << " Can not calculate residual width for passive plane !" << std::endl ;
    return 0.;
  }

  Double_t PointingRes = GetPointingResGBL( Ipl, UseInFit);
  //std::cout << "      GBL PR for plane " << Ipl << " = " << PointingRes << std::endl;

  Double_t Width;

  //std::cout << " Current intrinsic resolution used in GetWidthGBL(): " << _planeResolution[Ipl] << std::endl;


  if(UseInFit)
    Width=sqrt(_planeResolution[Ipl]*_planeResolution[Ipl]-PointingRes*PointingRes);
  else
    Width=sqrt(_planeResolution[Ipl]*_planeResolution[Ipl]+PointingRes*PointingRes);

  return Width;
}

Double_t AnaTel::GetDUTError()
{
  Double_t width=0.;

  if (_iDUT >=0 && _iDUT <  _nTelPlanes) width = GetPointingRes(_iDUT,_useDUT);

  return width;
}


Double_t AnaTel::GetDUTWidth()
{
  Double_t width=0.;

  if (_iDUT >=0 && _iDUT <  _nTelPlanes) width = GetWidth(_iDUT,_useDUT);

  return width;
}

void AnaTel::PrintPlanes()
{
  std::cout << "ntelPlanes = " << _nTelPlanes << std::endl;
  std::cout << "E = " << _eBeam << std::endl;
  std::cout << "s = " << _s << std::endl;
  std::cout << "dz = " << _dz << std::endl;
  std::cout << " size _sPoint = " << _sPoint.size() << std::endl;
  std::cout << " size _LoP = " << _listOfPoints.size() << std::endl;

  for(unsigned int ipl = 0; ipl < _nTelPlanes; ipl++)
  {
    std::cout << ipl << " thickness = " << _planeThickness[ipl] 
                   << " position  = " << _planePosition[ipl]
		   << " planeX0   = " << _planeX0[ipl]
		   << " planeReso = " << _planeResolution[ipl] 
		   << std::endl;
  }
  
}

// ==============================================================
//
// Private functions // DEPRECATED, USE GBL
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



