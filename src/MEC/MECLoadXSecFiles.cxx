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
  //delete Nieves_14_C12_2DGrids;
  // need to delete all the pointers to the grids that fill the vector

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
  //std::cout << "loading files from: " << data_dir << std::endl;  

  // define vector of grids with entries = num energies (92)
  //vector <BLI2DNonUnifGrid*>  Nieves_14_C12_2DGrids;

  // define dimensions of data in data files
  int nepoints = 92;
  int ntpoints = 31;
  int ncosthpoints = 31;
    
  // define arrays to fill from data files
  double enu;
  double nieves_e_array[nepoints];
  double nieves_costh_array[ncosthpoints]; 
  double nieves_t_array[ntpoints]; 
  double nieves_xsec_array[ntpoints*ncosthpoints];


  // loop over all 92 files for each configuration
  for (int i=0; i<(nepoints); i++){

    std::cout << "pulling file " << i << " of 92" << std::endl;

    // build filenames
    ostringstream datafile_14c;
    datafile_14c << data_dir << "/Nieves_14_C12_" << i+1 << ".dat";
  
    // make sure data files are available
    assert (! gSystem->AccessPathName(datafile_14c.str().c_str()));
  
    // read data files
    ReadNievesTCosthFile(datafile_14c.str(), ncosthpoints, ntpoints, nieves_costh_array, nieves_t_array, nieves_xsec_array, enu);

    nieves_e_array[i] = enu;

    // "save" xsecs into non uniform grid
    BLI2DNonUnifGrid *NievesGrid = new BLI2DNonUnifGrid(ncosthpoints, ntpoints, nieves_costh_array, nieves_t_array, nieves_xsec_array);
   
    // add grid to the array of grids
    Nieves_14_C12_2DGrids.push_back(NievesGrid);
    
  }
  
  // check the array of grids
  std::cout << "array size: " << Nieves_14_C12_2DGrids.size() << std::endl;
  std::cout << "evaluate a random point in one of the grids: " << Nieves_14_C12_2DGrids[31]->Evaluate(.5,.003) << std::endl;
  
  // test program...
  //testMECNievesLoadXSecFiles( Nieves_14_C12_2DGrids );

  std::cout << "Energies 0 - 1 - 91: " <<  nieves_e_array[0] << " - " <<  nieves_e_array[1] << " - " << nieves_e_array[91] << std::endl;

}
//___________________________________________________________________
void MECLoadXSecFiles::ReadNievesTCosthFile( 
  string filename, int ncosthpoints, int ntpoints, 
  double * costh_array, double * t_array, double * xsec_array, double enu)
{

  // ------ open and check file ---------- //

  // open file
  std::ifstream nieves_stream(filename.c_str(), ios::in);

  // check file exists
  if(!nieves_stream.good()){
    // log or output file read error //
    //std::cout << "bad file name: " << filename << std::endl;
    return ;
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
  double costh = 0, t = 0, xsec = 0 ;
  
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

//_____________________________________//
int MECLoadXSecFiles::EtoIndex(double ein)
{
  // right now, just rounding down.  Later real interpolation will be
  // implemented

  if (ein <= 0.1) return 0;
  else if (ein >= 30.) return 91;
  else {
    
    double e[] = 
      {0.1, 0.1025, 0.1050, 0.1075, 0.11,
       0.1125, 0.1150, 0.1175, 0.12, 0.1225,
       0.1250, 0.1275, 0.13, 0.1325, 0.1350,
       0.1375, 0.14, 0.1425, 0.145, 0.1475,
       0.15, 0.1525, 0.155, 0.1575, 0.16,
       0.1625, 0.165, 0.17, 0.175, 0.18,
       0.185, 0.19, 0.2, 0.205, 0.210, 
       0.22, 0.23, 0.24, 0.25, 0.26,
       0.275, 0.3, 0.325, 0.35, 0.375,
       0.4, 0.425, 0.45, 0.475, 0.5,
       0.525, 0.55, 0.6, 0.65, 0.7,
       0.75, 0.8, 0.85, 0.9, 0.95,
       1, 1.05, 1.1, 1.15, 1.2,
       1.25, 1.3, 1.4, 1.49, 1.4999999,
       1.5000001, 1.502, 1.5075, 1.513, 1.515,
       1.52, 1.525, 1.55, 1.575, 1.6,
       1.65, 1.8, 1.9, 2, 2.5, 
       3, 4, 6, 7.5, 10, 20, 30};

    for (int i=0; i< 92; i++){
      if (ein < e[i]) return i-1;
    }
  }
  return -1;
}

//_________________________________________________________________
double MECLoadXSecFiles::IndextoE(int index)
{
  // right now, just rounding down.  Later real interpolation will be
  // implemented

  if (index >= 91 ) return 30.;
  else {
    
    double e[] = 
      {0.1, 0.1025, 0.1050, 0.1075, 0.11,
       0.1125, 0.1150, 0.1175, 0.12, 0.1225,
       0.1250, 0.1275, 0.13, 0.1325, 0.1350,
       0.1375, 0.14, 0.1425, 0.145, 0.1475,
       0.15, 0.1525, 0.155, 0.1575, 0.16,
       0.1625, 0.165, 0.17, 0.175, 0.18,
       0.185, 0.19, 0.2, 0.205, 0.210, 
       0.22, 0.23, 0.24, 0.25, 0.26,
       0.275, 0.3, 0.325, 0.35, 0.375,
       0.4, 0.425, 0.45, 0.475, 0.5,
       0.525, 0.55, 0.6, 0.65, 0.7,
       0.75, 0.8, 0.85, 0.9, 0.95,
       1, 1.05, 1.1, 1.15, 1.2,
       1.25, 1.3, 1.4, 1.49, 1.4999999,
       1.5000001, 1.502, 1.5075, 1.513, 1.515,
       1.52, 1.525, 1.55, 1.575, 1.6,
       1.65, 1.8, 1.9, 2, 2.5, 
       3, 4, 6, 7.5, 10, 20, 30};

    return e[index];

  }
  return -1;
}



//_________________________________________________________________
// evaluate an element of vector and make 2d root hist.
void testMECNievesLoadXSecFiles( vector <BLI2DNonUnifGrid*>  Nieves2DGrids, int num, string filename){

  //could loop over all elements in vector and make a hist for each...

  double xmin=Nieves2DGrids[num]->XMin();
  double xmax=Nieves2DGrids[num]->XMax();
  double xstep = (xmax - xmin)/100.;
  double ymin=Nieves2DGrids[num]->YMin();
  double ymax=Nieves2DGrids[num]->YMax();
  double ystep = (ymax - ymin)/100.;

  TH2F *gridhist1 = new TH2F ( "gridhist1", "grid 1",20, ymin, ymax, 20, xmin, xmax);

  for (int i = 0 ; i < 20 ; i++ ){
    for (int j = 0 ; j < 20 ; j++) {

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
// Linear Interpolate - General (also extrapolates...)
double MECLoadXSecFiles::linearinterp(double x, double x1, double x2, double y1, double y2)
{
  return ( ((y1 - y2)/(x1 - x2)) * x ) + ( ((y2 * x1) - (y1 * x2))/(x1 - x2) );
}

//____________________________________________________________________
// XsecMax - linear interp max
double MECLoadXSecFiles::MaxXSec(double Enu)//(target,neuflavor)
{

  // get E index high and low
  int eindLow = EtoIndex(Enu);
  int eindHigh = eindLow+1;
  
  // Get E high and low
  double ELow = IndextoE(eindLow);
  double EHigh = IndextoE(eindHigh);

  // get maxxsec high and low
  double maxxsecLow = Nieves14C12()[eindLow]->ZMax();
  double maxxsecHigh = Nieves14C12()[eindHigh]->ZMax();

  // Linear interp for maxxsec

  return linearinterp(Enu, ELow, EHigh, maxxsecLow, maxxsecHigh);
}


//____________________________________________________________________
// XSec - linear interp
// get xsec(e)
double MECLoadXSecFiles::XSec(double Enu, double Costh, double T)//(target,neuflavor)
{

  // get E index high and low
  int eindLow = EtoIndex(Enu);
  int eindHigh = eindLow+1;
  
  // Get E high and low
  double ELow = IndextoE(eindLow);
  double EHigh = IndextoE(eindHigh);

  // get maxxsec high and low
  double xsecLow = Nieves14C12()[eindLow]->Evaluate(Costh,T);
  double xsecHigh = Nieves14C12()[eindHigh]->Evaluate(Costh,T);

  // Linear interp for maxxsec

  return linearinterp(Enu, ELow, EHigh, xsecLow, xsecHigh);
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
// interpolate xsec (E, t, costh)//


//_____________________________________//
// get maximum xsec (on a table)//
//  there is a zmax value that can be accessed?


// get maximum xsec (for any energy)


  

