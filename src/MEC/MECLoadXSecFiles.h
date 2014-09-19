//____________________________________________________________________________
/*!

\class    genie::INukeHadroData

\brief    Singleton class to load & save mec x-section tables - heavily
          based on INukeHadroData.cxx

\author   Jackie Schwehr

\created  March 12, 2014


\notes    


*/
//_____________________________________________________________________

#ifndef _MEC_LOAD_XSEC_FILES_H_
#define _MEC_LOAD_XSEC_FILES_H_

#include "Numerical/BLI2D.h"
#include "GHEP/GHepParticle.h"

namespace genie {

class MECLoadXSecFiles
{
 public:
  static MECLoadXSecFiles * Instance (void);
  
  // public functions
  void testMECNievesLoadXSecFiles(vector <BLI2DNonUnifGrid*>  Nieves2DGrids, int num, string filename);
  int EtoIndex(double);
  double IndextoE(int);
  double linearinterp(double, double, double, double, double);
  double MaxXSec(double);
  double XSec(double, double, double);

  // Public Variables?
  const vector <BLI2DNonUnifGrid *> Nieves14C12 (void) const {return Nieves_14_C12_2DGrids; }
  const int number (void) const {return 7;}



  
 private:
  MECLoadXSecFiles();
  // what is this one?
  MECLoadXSecFiles(const MECLoadXSecFiles & shx);
  ~MECLoadXSecFiles();

  // private functions
  void LoadXSecTables(void);
  void ReadNievesTCosthFile( string filename, int ncosthpoints, int ntpoints, double * costh_array, double * t_array, double * xsec_array, double enu);

  static MECLoadXSecFiles * fInstance;

  // private varables
  vector <BLI2DNonUnifGrid * >  Nieves_14_C12_2DGrids;
  
  // singleton cleaner
  struct Cleaner {
    void DummyMethodAndSilentCompiler(){}
    ~Cleaner(){
      if (MECLoadXSecFiles::fInstance !=0){
	delete MECLoadXSecFiles::fInstance;
	MECLoadXSecFiles::fInstance = 0;
      }
    }
  };
  friend struct Cleaner;
};

} // genie namespace

#endif // _MEC_LOAD_XSEC_FILES_H_
