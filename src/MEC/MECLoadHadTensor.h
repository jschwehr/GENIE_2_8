//____________________________________________________________________________
/*!

\class    genie::INukeHadroData

\brief    Singleton class to load & save mec hadron tensor tables - heavily
          based on MECLoadXSecFiles (which is based on INukeHadroData.cxx)

\author   Jackie Schwehr

\created  September 12, 2014

\notes    


*/
//_____________________________________________________________________

#ifndef _MEC_LOAD_HAD_TENSOR_H_
#define _MEC_LOAD_HAD_TENSOR_H_

#include "Numerical/BLI2D.h"
#include "GHEP/GHepParticle.h"

namespace genie {

class MECLoadHadTensor
{
 public:
  static MECLoadHadTensor * Instance (void);
  
  // public functions
  double XSec(int, double, double, double); // pdg, E, T, Costheta;
  double MaxXSec(double);// Enu

  // Public Variables
  const vector <BLI2DNonUnifGrid *> HadTensorC12 (void) const {return HadTensor_C12_2DGrids;}
  const vector <double> MaxXSecC12 (void) const {return MaxXSecvect;}
  const vector <double> EnuValues (void) const {return Enuvect;}

 private:

  // private offical stuff
  MECLoadHadTensor();
  MECLoadHadTensor(const MECLoadHadTensor & shx);
  ~MECLoadHadTensor();
  static MECLoadHadTensor * fInstance;

  // private functions
  void LoadTensorTables(void);
  void ReadHadTensorqzq0File(string filename, int nwpoints, int nqzpoints, int nq0points, double hadtensor_w_array[][57600]);
  void MakeMaxXSecTables();

  // private varables
  vector <BLI2DNonUnifGrid * >  HadTensor_C12_2DGrids;
  vector <double> MaxXSecvect;
  vector <double> Enuvect;

  // singleton cleaner
  struct Cleaner {
    void DummyMethodAndSilentCompiler(){}
    ~Cleaner(){
      if (MECLoadHadTensor::fInstance !=0){
	delete MECLoadHadTensor::fInstance;
	MECLoadHadTensor::fInstance = 0;
      }
    }
  };
  friend struct Cleaner;
};

} // genie namespace

#endif // _MEC_LOAD_HAD_TENSOR_H_
