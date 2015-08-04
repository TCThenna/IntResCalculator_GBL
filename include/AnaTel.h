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

    AnaTel(std::string GeomFile);
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
    //gbl::GblPoint* point;
    std::vector<gbl::GblPoint> listOfPoints;

    //ClassDef(AnaTel,1);

};


#endif

