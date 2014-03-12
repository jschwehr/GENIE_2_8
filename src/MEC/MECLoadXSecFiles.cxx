//_______________________________________________________________
/*



*/
//_______________________________________________________________


// --- Includes --- //

// 
#include <cmath>
#include <cstdlib>
#include <sstream>
#include <string>
#include <fstream>
#include <cassert>
#include <string>
#include <iostream>
#include <vector>

#include "Conventions/Constants.h"
#include "HadronTransport/INukeHadroData.h"
#include "Messenger/Messenger.h"
#include "Numerical/Spline.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"
#include "Numerical/BLI2D.h"
#include "MEC/MECLoadXSecFiles.h"

#include <TSystem.h>
#include <TFile.h>
#include <TGraph2D.h>
#include <TH2.h>

using std::ostringstream;
using std::istream;
using std::ios;
using std::endl;

using namespace genie;
using namespace genie::constants;

//________________________________________________________________
// instance
MECLoadXSecFiles * MECLoadXSecFiles::fInstance = 0;

//_________________________________________________________________
// define variables

//________________________________________________________________
// keep it simple
MECLoadXSecFiles::MECLoadXSecFiles()
{
  this->LoadXSecTables();
  fInstance = 0;
}
//________________________________________________________________
// destructor
MECLoadXSecFiles::~MECLoadXSecFiles()
{
  // xsec arrays
  delete Nieves_14_C12_Graphs;

}
//________________________________________________________________
// instance
MECLoadXSecFiles * MECLoadXSecFiles::Instance()
{
  if(fInstance == 0){
    //LOG("MECLoadXSecFiles", pINFO) << "MECLoadXSecFiles late initialization";
    static MECLoadXSecFiles::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new MECLoadXSecFiles;
  }
  return fInstance;
}
//_________________________________________________________________
// load xsec tables
void MECLoadXSecFiles::LoadXSecTables(void)
{

  // define directory
  string data_dir = string("/home/jackie/work/neut/5.3.2/crsdat");
  std::cout << "loading files from: " << data_dir << std::endl;  

  // define vector of grids with entries = num energies (92)
  vector <BLI2DNonUnifGrid*>  Nieves_14_C12_2DGrids;

  // loop over all 92 files for each configuration
  for (int i=1; i<93; i++){

    std::cout << "pulling file " << i << " of 92" << std::endl;

    // build filenames
    ostringstream datafile_14c;
    datafile_14c << data_dir << "/Nieves_14_C12_" << i << ".dat";
  
    // make sure data files are available
    assert (! gSystem->AccessPathName(datafile_14c.str().c_str()));
  
    // define dimensions of data in data files
    int ntpoints = 31;
    int ncosthpoints = 31;

    // define arrays to fill from data files
    double nieves_costh_array[ncosthpoints]; 
    double nieves_t_array[ntpoints]; 
    double nieves_xsec_array[ntpoints*ncosthpoints];
  
    // read data files
    ReadNievesTCosthFile(datafile_14c.str(), ncosthpoints, ntpoints, nieves_costh_array, nieves_t_array, nieves_xsec_array);

    // "save" xsecs into non uniform grid
    BLI2DNonUnifGrid *NievesGrid = new BLI2DNonUnifGrid(ncosthpoints, ntpoints, nieves_costh_array, nieves_t_array, nieves_xsec_array);
   
    // add grid to the array of grids
    Nieves_14_C12_2DGrids.push_back(NievesGrid);

  }
  
  // check the array of grids
  //std::cout << "array size: " << Nieves_14_C12_2DGrids.size() << std::endl;
  //std::cout << "evaluate a random point in one of the grids: " << Nieves_14_C12_2DGrids[31]->Evaluate(.5,.003) << std::endl;
  
  // test program...
  //testMECNievesLoadXSecFiles( Nieves_14_C12_2DGrids );

}
//___________________________________________________________________
void MECLoadXSecFiles::ReadNievesTCosthFile( 
  string filename, int ncosthpoints, int ntpoints, 
  double * costh_array, double * t_array, double * xsec_array)
{

  // ------ open and check file ---------- //

  // open file
  std::ifstream nieves_stream(filename.c_str(), ios::in);

  // check file exists
  if(!nieves_stream.good()){
    // log or output file read error //
    //std::cout << "bad file name: " << filename << std::endl;
    return;
  }
  
  //std::cout << "opening filename: " << filename << std::endl;

  // ---- TO DO ---- 
  // check if file is of correct format (enu, costh, t, xsec)
  // ?correct energy?
  // correct number of columns?
  // correct number of rows? (ntpoints (31) * ncosthpoints (31) = 961)
  // (and arrays provided are the correct size...?)
    
  //--------- read file ---------- //

  // define columns
  double enu = 0, costh = 0, t = 0, xsec = 0 ;
  
  // loop over rows
  for (int i = 0; i < (ntpoints*ncosthpoints); i++){
  
    // read row
    nieves_stream >> enu >> costh >> t >> xsec ; 
    
    // fill arrays
    //    only fill costh once per costh value
    if (i%ncosthpoints==0) costh_array[(i/ncosthpoints)]=costh;
    //    only fill t for the first set of t values
    if (i < ntpoints) t_array[i]=t;
    //    fill xsec every time.
    xsec_array[i]=xsec;

  }// end for i
}
//_________________________________________________________________
// evaluate an element of vector and make 2d root hist.
void testMECNievesLoadXSecFiles( vector <BLI2DNonUnifGrid*>  Nieves2DGrids){

  //could loop over all elements in vector and make a hist for each...

  double xmin=Nieves2DGrids[61]->XMin();
  double xmax=Nieves2DGrids[61]->XMax();
  double xstep = (xmax - xmin)/100.;
  double ymin=Nieves2DGrids[61]->YMin();
  double ymax=Nieves2DGrids[61]->YMax();
  double ystep = (ymax - ymin)/100.;

  TH2F *gridhist1 = new TH2F ( "gridhist1", "grid 1", 100, ymin, ymax, 100, xmin, xmax);

  for (int i = 0 ; i < 100 ; i++ ){
    for (int j = 0 ; j < 100 ; j++) {

      gridhist1->SetBinContent(j,i,
			      Nieves2DGrids[61]->Evaluate( xmin+(i*xstep), 
							   ymin+(j*ystep) ));

    }
  }  

  TFile f("testhist.root","recreate");
  gridhist1->Write();
  f.Close();
}
//___________________________________________________________________
// future plans/potentially usefull functions?
  
//_______________________________________//
// interpolate target //
//  given a target, and the targets available, find way to scale everything to the new target //

//______________________________________//
// interpolate E //
//  given a neutrino energy, and the energies available, find the way to scale everything to the new energy //
//  also identify the two bracketing energies for future reference?

//______________________________________//
// interpolate xsec //


//_____________________________________//
// get maximum xsec //


//______________________________________//
// (get minimum xsec?) //
  

