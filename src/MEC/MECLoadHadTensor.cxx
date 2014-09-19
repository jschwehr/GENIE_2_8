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
#include "MEC/MECLoadHadTensor.h"

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
MECLoadHadTensor * MECLoadHadTensor::fInstance = 0;

//_________________________________________________________________
// define variables

//________________________________________________________________
// keep it simple
MECLoadHadTensor::MECLoadHadTensor()
{
  this->LoadTensorTables();
  this->MakeMaxXSecTables();
  fInstance = 0;
}
//________________________________________________________________
// destructor
MECLoadHadTensor::~MECLoadHadTensor()
{

}
//________________________________________________________________
// instance
MECLoadHadTensor * MECLoadHadTensor::Instance()
{
  if(fInstance == 0){
    //LOG("MECLoadHadTensor", pINFO) << "MECLoadHadTensor late initialization";
    static MECLoadHadTensor::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new MECLoadHadTensor;
  }
  return fInstance;
}
//_________________________________________________________________
// load xsec tables
void MECLoadHadTensor::LoadTensorTables(void)
{

  // define directory
  string data_dir = string("/home/jackie/work/genie/rikcode/mec-20140813");

  // define dimensions of data in data files
  int nwpoints = 5;
  int nq0points = 240;
  int nqzpoints = 240;

  // define arrays to fill from data files
  double hadtensor_q0_array[nq0points];
  double hadtensor_qz_array[nqzpoints];
  double hadtensor_w_array[nwpoints][57600];

  
  // fill q0 array (in GeV)
  // 240 5MeV bins
  for (int a = 0; a < 240; a++){
    hadtensor_q0_array[a]=a*0.005;
  }

  // fill qz array (in GeV)
  // 240 5MeV bins
  for (int a = 0; a < 240; a++){
    hadtensor_qz_array[a]=a*0.005;
  }
 
  // build filenames
  ostringstream datafile_c12;
  datafile_c12 << data_dir << "/HadTensor240-C12Fullpn.dat";
  
  // make sure data files are available
  assert (! gSystem->AccessPathName(datafile_c12.str().c_str()));
  
  // read data file
  ReadHadTensorqzq0File(datafile_c12.str(), nwpoints, nqzpoints, nq0points,hadtensor_w_array);
  
  //loop over all 5 tensors 
  for (int i=0; i<(nwpoints); i++){
   
    // "save" xsecs into non uniform grid
    BLI2DNonUnifGrid *HadTensorGrid = new BLI2DNonUnifGrid(nq0points, nqzpoints, hadtensor_q0_array, hadtensor_qz_array, hadtensor_w_array[i]);
    
    // add grid to the array of grids
    HadTensor_C12_2DGrids.push_back(HadTensorGrid);
  }
  
}
//___________________________________________________________________
void MECLoadHadTensor::ReadHadTensorqzq0File( 
					     string filename, int nwpoints, int nqzpoints, int nq0points, double hadtensor_w_array[][57600])
{
  // ------ open and check file ---------- //

  // open file
  std::ifstream tensor_stream(filename.c_str(), ios::in);

  // check file exists
  if(!tensor_stream.good()){
    // log or output file read error //
    //std::cout << "bad file name: " << filename << std::endl;
    return ;
  }
  
  // -- to do -- //
  // check if file is of correct format/size (nwpoints, nq0points, nqzpoints)

  //--------- read file ---------- //

  double temp;
  
  for (int ij = 0; ij < (nqzpoints*nq0points); ij++){  
    for (int k = 0; k < nwpoints; k++){
      tensor_stream >> temp;
      hadtensor_w_array[k][ij]=temp;
    }}

}

//____________________________________________________________________
// Xsec
double MECLoadHadTensor::XSec(int pdgnu, double Enu, double Tmu, double Costheta)
{

  //--- variables and constants ---//
  double xsec;
  double Mmu = 0.105658357;
  double pi=3.141592653589793;
  TLorentzVector v4Mu;
  TLorentzVector v4Nu = TLorentzVector(0,0,Enu,Enu);   // assuming traveling along z:
  TLorentzVector v4q;
  double q0nucleus;
  double Gfermi=1.1664e-5;
  double facconv = std::pow(0.19733,2)*1e15;

  double Qvalue;
  if (pdgnu > 0) Qvalue = 16.827/1000.; // GeV neutrino
  else Qvalue = 13.880/1000.; // GeV anti-neutrino

  //--- kinetics ---//
  //angles
  double Sintheta = std::sqrt( 1. - std::pow(Costheta,2));
  double Cosh = std::cos(std::acos(Costheta)/2.);
  double Sinh = std::sin(std::acos(Costheta)/2.);
  //muon
  v4Mu.SetE( Tmu + Mmu );
  //q
  double q0 = v4Nu.E() - v4Mu.E();
  q0nucleus = q0 - Qvalue;

  // Define some calculation placeholders
  double part1, part2;
  double modkprime ;
  double w[5];
  double wtotd[5];

  //--- event selection loop ---//
  // if energy was transfered to the nucleus... commence!
  if (q0nucleus > 0){
    modkprime = std::sqrt(std::pow(v4Mu.E(),2)-std::pow(Mmu,2));
    v4Mu.SetX(modkprime*Sintheta); 
    v4Mu.SetY(0);
    v4Mu.SetZ(modkprime*Costheta);

    //q: v4q = v4Nu - v4Mu;
    v4q.SetE(q0nucleus);
    v4q.SetX(v4Nu.X() - v4Mu.X());
    v4q.SetY(v4Nu.Y() - v4Mu.Y());
    v4q.SetZ(v4Nu.Z() - v4Mu.Z());
    //remember: qm = v4q.Vect().Mag()
  
    // pull elements from tensor tables
    for (int i=0 ; i < 5; i++){
      wtotd[i]=HadTensor_C12_2DGrids[i]->Evaluate(v4q.Vect().Mag(),v4q.E());
    }

    // calculate!
    w[0]=wtotd[2]/2.;
    w[1]=(wtotd[0]+wtotd[2]+(std::pow(v4q.E(),2)/std::pow(v4q.Vect().Mag(),2)
			     *(wtotd[4]-wtotd[2]))
	  -(2.*v4q.E()/v4q.Vect().Mag()*wtotd[1]))/2.;
    w[2]=wtotd[3]/v4q.Vect().Mag();
    w[3]=(wtotd[4]-wtotd[2])/(2.*std::pow(v4q.Vect().Mag(),2));
    w[4]=(wtotd[1]-v4q.E()/v4q.Vect().Mag()*(wtotd[4]-wtotd[2]))/v4q.Vect().Mag();
    // adjust for anti neutrinos
    if (pdgnu < 0) w[3] = -1. * w[3];

    // calculate cross section, in parts

    part1 = w[0]*Costheta - w[1]/2.*Costheta 
      + w[2]/2.*(v4Nu.E()+modkprime-(v4Mu.E()+v4Nu.E())*Costheta)
      + w[3]/2.*(std::pow(Mmu,2)*Costheta+2.*v4Mu.E()
		 *(v4Mu.E()+modkprime)*std::pow(Sinh,2))
      - w[4]*(v4Mu.E()+modkprime)/2.;
    part2 = 2.*w[0]*std::pow(Sinh,2) + w[1]*std::pow(Cosh,2)
      - w[2]*(v4Mu.E()+v4Nu.E())*std::pow(Sinh,2)
      + std::pow(Mmu,2)/(v4Mu.E()*(v4Mu.E()+modkprime))*part1;
    xsec = modkprime*v4Mu.E()*std::pow(Gfermi,2)*2./pi*part2*facconv;
  }

  else {//nothing transfered to nucleus thus no mec
    xsec = 0.;
  }
  
  return xsec;
}

//____________________________________________________________________
// make max xsec table
//  -- this should be made once and saved as a text file to read in - takes too long!
void MECLoadHadTensor::MakeMaxXSecTables(){

  // save as 2 vectors

  // -- variables -- //
  // x->Tmu, y->CosTheta, z->E
  double xbins=1000;
  double ybins=100;
  double zbins=100;
  //double xbins=100;
  //double ybins=10;
  //double zbins=10;

  double xmin = 0;
  double xmax = 5;
  double ymin = -1;
  double ymax = 1;
  double zmin = 0;
  double zmax = 5;

  TH2F *temphist = new TH2F ( "temphist", "temphist", xbins, xmin, xmax, ybins, ymin, ymax);

  int apdgnu = 14;
  double aEnu;
  double aTmu; 
  double aCostheta;
  double xsec;

  // -- make sure the vectors are empty before filling -- //
  MaxXSecvect.clear();
  Enuvect.clear();
  
  // -- loop over E, Tmu, Costheta and pull xsec -- //
  for (int z = 0; z< zbins; z++){
    aEnu = (zmin + ((zmax - zmin)/zbins)/2. )+ z * (zmax - zmin)/zbins ;
    
    for (int i = 0; i < xbins; i++){
      aTmu = (xmin + ((xmax - xmin)/xbins)/2. )+ i * (xmax - xmin)/xbins ;
      
      for (int j = 0; j < ybins; j++){
	aCostheta = (ymin + ((ymax - ymin)/ybins)/2. )+ j * (ymax - ymin)/ybins;

	// get xsec and store in 2d histogram
	xsec =  XSec(apdgnu, aEnu, aTmu, aCostheta);
	temphist->SetBinContent( i, j, xsec);
	
      }}
    
    // fill xsec and enu vectors
    MaxXSecvect.push_back(temphist->GetBinContent(temphist->GetMaximumBin()));
    Enuvect.push_back(aEnu);
  }
   
}

//________________________________________________________________
// get max xsec - interpolate/extrapolate.
double MECLoadHadTensor::MaxXSec(double Enu)
{

  // -- get the right vector bins -- //
  // get bins:
  int ebinhi=-1;
  int iterate = 0;
  int ebinlo;
  double elo;
  double xseclo;

  while (ebinhi == -1){
    if(Enuvect[iterate]>Enu) ebinhi=iterate;
    else if (iterate == (int)Enuvect.size()) ebinhi=0;
    else iterate++;
  }

  if (ebinhi == 0) {
    ebinlo = 0;  
    elo = Enuvect[ebinlo]; 
    xseclo = MaxXSecvect[ebinlo];
  }
  else {
    ebinlo = ebinhi-1.;
    elo = Enuvect[ebinlo]; 
    xseclo = MaxXSecvect[ebinlo];
  }
  
  double ehi = Enuvect[ebinhi];
  double xsechi = MaxXSecvect[ebinhi];

  // Linear inter/extra-polation of xsec at ENEU and find maximum value of xsec
  // Formula: y = (((y1-y2)/(x1-x2)) * x ) + (((y2*x1)-(y1*x2)) / (x1-x2))
  // where y = extra/inter-polated xsec, Enu = x
  // y1, y2 = xsecs from tables at Enu = x1, x2.
  

  return ( ( xseclo - xsechi )/( elo - ehi ) * Enu ) + ( ( xsechi*elo - xseclo*ehi ) / ( elo - ehi ) ); 

}

