#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GaussIntegrator.h"



int energyshift(){

  cout << "Energy shift calculator" << endl;

  double central = 6.;
  double d_up = 0.1; 
  double d_down = 1.4; 
  double lowerlimit = central - d_down;
  double upperlimit = central + d_up;
  double L = d_up + d_down;;

  TF1 f("1 over x", "1/x", lowerlimit, upperlimit);
 
  double integral = f.Integral(lowerlimit,upperlimit);


 
  //double integral = ig.Integral(lowerlimit, upper) << endl;
  double inv_int = 1./integral;

  cout << "L = " << L << endl;
  cout << "int f = " << integral << endl;
  cout << "int f / L = " << integral / L << endl;
  cout << "inv normalint = " << L / integral << endl;
  cout << "\% of central = " << L / integral / central << endl;

  return 0;

}
