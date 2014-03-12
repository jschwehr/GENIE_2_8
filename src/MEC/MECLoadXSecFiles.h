//____________________________________________________________________________
/*!

\class    genie::INukeHadroData

\brief    Singleton class to load & serve mec x-section tables - heavily
          based on INukeHadroData.cxx

\author   Jackie Schwehr

\created  March 12, 2014

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
  void testMECNievesLoadXSecFiles(vector <BLI2DNonUnifGrid*>  Nieves2DGrids);

  // elements worth returning!
  const vector <BLI2DNonUnifGrid> * const Nieves14C12 (void) const {return Nieves_14_C12_Graphs; }
  
 private:
  MECLoadXSecFiles();
  // what is this one?
  MECLoadXSecFiles(const MECLoadXSecFiles & shx);
  ~MECLoadXSecFiles();

  // private functions
  void LoadXSecTables(void);
  void ReadNievesTCosthFile( string filename, int ncosthpoints, int ntpoints, double * costh_array, double * t_array, double * xsec_array);

  static MECLoadXSecFiles * fInstance;

  // private varables
  vector <BLI2DNonUnifGrid> * Nieves_14_C12_Graphs;
  
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
