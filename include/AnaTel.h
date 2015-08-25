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

    // default constructor
    AnaTel() = default;

    // Constructor: create telescope from geometry file
    AnaTel(std::string GeomFile, double energy);

    // Constructor: create telescope from geometry file, additional integer value to turn of certain plane
    AnaTel(std::string GeomFile, unsigned int iPlaneOff);

    // Constructor: create telescope from geometry file, additional integer value to turn of certain plane
    AnaTel(std::string GeomFile, unsigned int * PlanesOff);

    // Constructor: create telescope "by hand"
    AnaTel(Int_t Npl);

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

    Double_t GetPointingRes(Int_t Ipl, Bool_t UseInFit);
    Double_t GetPointingResGBL(Int_t Ipl, Bool_t UseInFit);

    Double_t GetWidth(Int_t Ipl, Bool_t UseInFit);
    Double_t GetWidthGBL(Int_t Ipl, Bool_t UseInFit);

    Double_t GetDUTError();

    Double_t GetDUTWidth();

    gbl::GblPoint getPoint(double step, double res, TVectorD wscat, bool IsPlane, bool has_meas = true);
    gbl::GblPoint getPoint(double step, TVectorD wscat, bool IsPlane);

    void PrintID()
    {
      for (int i = 0; i < _nTelPlanes; i++)
        std::cout << _ID[i] << " ";
      std::cout << std::endl;
    }

    void PrintPlanes();

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

    double * _planeDist;
    double * _planeScat;
    double * _tempScat1;
    double * _tempScat2;
    double * _tempScat3;

    double * _fitX;
    double * _fitArray;
    double * _pointingResolution;

    bool   _useBeamConstraint;
    double _beamSpread;

    int DoAnalFit(double * , double *);
    int GaussjSolve(double * ,double * ,int );

    // GBL additions
    //gbl::GblPoint* point;
    std::vector<gbl::GblPoint> _listOfPoints;
    std::vector<double> _sPoint;

    double _dz;
    double _s; // arc

    //TVectorD _wscatSi(2);
    //TVectorD _wscatDUT(2);
    //TVectorD _wscatAir(2);

    double _X0_Si_frac;
    double _X0_DUT_frac;
    double _X0_Air_frac;

    static constexpr double _X0_Si = 93.65; //  mm
    static constexpr double _X0_Air = 304200; //  [mm] from PDG   
    static constexpr double _X0_Kapton =  285.6; // mm

    std::vector<int> _ID; // contains the position of the M26 planes in listOfPoints


};


#endif

