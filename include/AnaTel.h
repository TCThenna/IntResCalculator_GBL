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

#ifndef AnaTel_Included
#define AnaTel_Included

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
//#include "TF1.h"

//using namespace gbl;

class AnaTel {

public:
  
  //virtual ~AnaTel();

  AnaTel(char * GeomFile);
  // Constructor: create telescope from geometry file

  AnaTel(Int_t Npl);
  // Constructor: create telescope "by hand"

  // Dectructor
  ~AnaTel() = default;

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

    // GBL additions
    gbl::GblPoint* point;
    std::vector<gbl::GblPoint> listOfPoints;

  //ClassDef(AnaTel,1);

};

AnaTel::AnaTel(char * GeomFile)
{
  // GBL
  //listOfPoints.reserve(3*6);
  //TVectorD wscatSi(2);
  //    double tetSi = 0.1;
  //    wscatSi[0] = 1.0 / ( tetSi * tetSi ); // weight
  //    wscatSi[1] = 1.0 / ( tetSi * tetSi );
  //point->getPoint(150, wscat);

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


#endif

