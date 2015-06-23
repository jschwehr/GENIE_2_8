//_______________________________________________________________
/*

January 14th, 2015 - All neutrino flavors, Carbon and Oxygen targets

Hardcoded Inputs:
   Hadron tensor file directory: Line 89
   ** creates MaxXSec files in directory in which it is run, and then looks for these files later.  To speed up processing, do not delete/move them.  (to do - put these files in the directory with the hadron tensor files?)

   ** 587 - max xsec files hard linked location - should use global variables

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
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/CmdLnArgParser.h"
#include "Numerical/BLI2D.h"
#include "MECTensor/MECLoadHadTensor.h"

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
// 
MECLoadHadTensor * MECLoadHadTensor::fInstance = 0;
//________________________________________________________________
//
MECLoadHadTensor::MECLoadHadTensor(int targetpdg, int nupdg)
{
  this->LoadTensorTables(targetpdg);
  //this->WriteMaxXSecTables(targetpdg, nupdg);
  //this->ReadMaxXSecTables(targetpdg, nupdg);
  //std::cout << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Worked!~~~~~" << endl;
  fInstance = 0;
}
//________________________________________________________________
// 
MECLoadHadTensor::~MECLoadHadTensor()
{
}
//________________________________________________________________
//
MECLoadHadTensor * MECLoadHadTensor::Instance(int targetpdg, int nupdg)
{
  if(fInstance == 0){
    //LOG("MECLoadHadTensor", pINFO) << "MECLoadHadTensor late initialization";
    static MECLoadHadTensor::Cleaner cleaner;
    cleaner.DummyMethodAndSilentCompiler();
    fInstance = new MECLoadHadTensor(targetpdg,nupdg);
  }
  return fInstance;
}
//_________________________________________________________________
// Create the hadron tensor tables
void MECLoadHadTensor::LoadTensorTables(int targetpdg)
{

  // define directory of hadron tensor files
  // Put this in the xml configuration
  // string data_dir = string("/home/jackie/work/genie/GENIE_2_8/data/tensors");
  string data_dir = string("/home/jackie/work/genie/GENIE_2_8/data/evgen/mectensor/nieves");
  //string data_dir = string("/usr/local/physics/hep/genie-2.8.6-20150204/src/MEC");
  //string data_dir = string("./");

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
    hadtensor_q0_array[a]=double(a+1)*0.005;
  }

  // fill qz array (in GeV)
  // 240 5MeV bins
  for (int a = 0; a < 240; a++){
    hadtensor_qz_array[a]=double(a+1)*0.005;
  }

  for (int tables = 0; tables < 4; tables++){ 
  // build filenames
  ostringstream datafile;
  if (targetpdg == 1000060120){//carbon
    if(tables == 0) datafile << data_dir << "/HadTensor240-C12FullAll-20150210.dat";
    if(tables == 1) datafile << data_dir << "/HadTensor240-C12Fullpn-20150210.dat";
    if(tables == 2) datafile << data_dir << "/HadTensor240-C12DeltaAll-20150210.dat";
    if(tables == 3) datafile << data_dir << "/HadTensor240-C12Deltapn-20150210.dat";
  }
  else if (targetpdg == 1000080160) {//oxygen
    if(tables == 0) datafile << data_dir << "/HadTensor240-O16FullAll-20150210.dat";
    if(tables == 1) datafile << data_dir << "/HadTensor240-O16Fullpn-20150210.dat";
    if(tables == 2) datafile << data_dir << "/HadTensor240-O16DeltaAll-20150210.dat";
    if(tables == 3) datafile << data_dir << "/HadTensor240-O16Deltapn-20150210.dat";
  }
  else if (targetpdg == 1000200400) {//calcium 40 isoscalar similar to argon 40
    if(tables == 0) datafile << data_dir << "/HadTensor240-Ca40FullAll-20150210.dat";
    if(tables == 1) datafile << data_dir << "/HadTensor240-Ca40Fullpn-20150210.dat";
    if(tables == 2) datafile << data_dir << "/HadTensor240-Ca40DeltaAll-20150210.dat";
    if(tables == 3) datafile << data_dir << "/HadTensor240-Ca40Deltapn-20150210.dat";
  }
  else { 
    std::cout << "bad target" << std::endl;
    assert(false);
  }

  // make sure data files are available
  assert (! gSystem->AccessPathName(datafile.str().c_str()));
  
  // read data file
  ReadHadTensorqzq0File(datafile.str(), nwpoints, nqzpoints, nq0points,hadtensor_w_array);
  
  //loop over all 5 tensors 
  for (int i=0; i<(nwpoints); i++){
   
    // "save" xsecs into non uniform grid
    genie::BLI2DNonUnifGrid *HadTensorGrid = new genie::BLI2DNonUnifGrid(nq0points, nqzpoints, hadtensor_q0_array, hadtensor_qz_array, hadtensor_w_array[i]);
    
    // add grid to the array of grids
    if (targetpdg == 1000060120){//carbon
      if(tables == 0)       HadTensorFullAll_C12_2DGrids.push_back(HadTensorGrid);
      if(tables == 1)       HadTensorFullpn_C12_2DGrids.push_back(HadTensorGrid);
      if(tables == 2)       HadTensorDeltaAll_C12_2DGrids.push_back(HadTensorGrid);
      if(tables == 3)       HadTensorDeltapn_C12_2DGrids.push_back(HadTensorGrid);

    }
    else if (targetpdg ==1000080160){//oxygen
      if(tables == 0)       HadTensorFullAll_O16_2DGrids.push_back(HadTensorGrid);
      if(tables == 1)       HadTensorFullpn_O16_2DGrids.push_back(HadTensorGrid);
      if(tables == 2)       HadTensorDeltaAll_O16_2DGrids.push_back(HadTensorGrid);
      if(tables == 3)       HadTensorDeltapn_O16_2DGrids.push_back(HadTensorGrid);
    }
    else if (targetpdg ==1000200400){//oxygen
      if(tables == 0)       HadTensorFullAll_Ca40_2DGrids.push_back(HadTensorGrid);
      if(tables == 1)       HadTensorFullpn_Ca40_2DGrids.push_back(HadTensorGrid);
      if(tables == 2)       HadTensorDeltaAll_Ca40_2DGrids.push_back(HadTensorGrid);
      if(tables == 3)       HadTensorDeltapn_Ca40_2DGrids.push_back(HadTensorGrid);
    }
   
  }
  
  
  }// end loop over tables (probably need to clear things here...)  
}
//___________________________________________________________________
void MECLoadHadTensor::ReadHadTensorqzq0File( 
					     string filename, int nwpoints, int nqzpoints, int nq0points, double hadtensor_w_array[][57600]){
  // ------ open and check file ---------- //

  // open file
  std::ifstream tensor_stream(filename.c_str(), ios::in);

  // check file exists
  if(!tensor_stream.good()){
    // log or output file read error //
    //std::cout << "bad file name: " << filename << std::endl;
    return ;
  }
  
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
double MECLoadHadTensor::XSec(int targetpdg, int nupdg, double Enu, double Tmu, double Costheta,  vector <genie::BLI2DNonUnifGrid *> HadTensor)
{
  // These units should deliver 10^{41} cm^2 / GeV for d2sigma/(dTmu dcos_mu)

  
  //--- variables and constants ---//

  // lepton:
  double Mlep = 0.0;
  if( nupdg == 12 || nupdg == -12 ){
    Mlep = 0.000510998928;
  }
  else if( nupdg == 14 || nupdg == -14 ){
    Mlep = 0.105658357;
  }
  else if( nupdg == 16 || nupdg == -16){
    Mlep = 1.77682;
  }
  double xsec;
  double pi=3.141592653589793;
  TLorentzVector v4lep;
  TLorentzVector v4Nu = TLorentzVector(0,0,Enu,Enu);   // assuming traveling along z:
  TLorentzVector v4q;
  double q0nucleus;
  double Gfermi=1.1664e-5;
  double facconv = 0.0389391289e15; // std::pow(0.19733,2)*1e15;
  // These units should deliver 10^{41} cm^2 / GeV for d2sigma/(dTmu dcos_mu)

  double Qvalue;
  // From Table I in Nieves et al. PRC 70 055503 (2004)
  if(targetpdg == 1000060120){
    if (nupdg > 0) Qvalue = 16.827/1000.; // GeV neutrino
    else Qvalue = 13.880/1000.; // GeV anti-neutrino
  } else if(targetpdg == 1000080160){
    if (nupdg > 0) Qvalue = 14.906/1000.; // GeV neutrino
    else Qvalue = 10.931/1000.; // GeV anti-neutrino
  } else if(targetpdg == 1000200400){
    if (nupdg > 0) Qvalue = 13.809/1000.; // GeV neutrino
    else Qvalue = 1.822/1000.; // GeV anti-neutrino
  } else {
    LOG("MEC", pERROR) << "Bad nucleus " << targetpdg << " did you forget to implement other nuclei?";
    assert(false);
  }
  

  //--- kinetics ---//
  //angles
  double Sintheta = 1. - Costheta * Costheta;
  if(Sintheta < 0.0)Sintheta = 0.0;
  else Sintheta = TMath::Sqrt(Sintheta);   
  
  double Cosh = TMath::Cos(TMath::ACos(Costheta)/2.);
  double Sinh = TMath::Sin(TMath::ACos(Costheta)/2.);
  //lepton
  v4lep.SetE( Tmu + Mlep );
  //q
  double q0 = v4Nu.E() - v4lep.E();
  q0nucleus = q0 - Qvalue;

  // Define some calculation placeholders
  double part1, part2;
  double modkprime ;
  double w[5];
  double wtotd[5];

  //--- event selection loop ---//
  if (q0nucleus > 0){ // true if energy was transfered to the nucleus
    modkprime = std::sqrt(TMath::Power(v4lep.E(),2)-TMath::Power(Mlep,2));
    v4lep.SetX(modkprime*Sintheta); 
    v4lep.SetY(0);
    v4lep.SetZ(modkprime*Costheta);

    //q: v4q = v4Nu - v4lep;
    v4q.SetE(q0nucleus);
    v4q.SetX(v4Nu.X() - v4lep.X());
    v4q.SetY(v4Nu.Y() - v4lep.Y());
    v4q.SetZ(v4Nu.Z() - v4lep.Z());
    //remember: qm = v4q.Vect().Mag()
  
    /*
    if (targetpdg == 1000060120){
      for (int i=0 ; i < 5; i++){
	wtotd[i]=HadTensor[i]->Evaluate(v4q.Vect().Mag(),v4q.E());
      }
    }
    else if (targetpdg == 1000080160){
      for (int i=0 ; i < 5; i++){
	wtotd[i]=HadTensor[i]->Evaluate(v4q.Vect().Mag(),v4q.E());
      }
    }
    */
    for (int i=0 ; i < 5; i++){
      wtotd[i]=HadTensor[i]->Evaluate(v4q.Vect().Mag(),v4q.E());
    }

    // calculate hadron tensor components
    // these are footnote 2 of Nieves PRC 70 055503
    w[0]=wtotd[2]/2.;
    w[1]=(wtotd[0]+wtotd[2]+(TMath::Power(q0,2)/TMath::Power(v4q.Vect().Mag(),2)
			     *(wtotd[4]-wtotd[2]))
	  -(2.*q0/v4q.Vect().Mag()*wtotd[1]))/2.;
    w[2]=wtotd[3]/v4q.Vect().Mag();
    w[3]=(wtotd[4]-wtotd[2])/(2.*TMath::Power(v4q.Vect().Mag(),2));
    w[4]=(wtotd[1]-(q0/v4q.Vect().Mag()*(wtotd[4]-wtotd[2])))/v4q.Vect().Mag();

    // adjust for anti neutrinos
    if (nupdg < 0) w[3] = -1. * w[3];

    // calculate cross section, in parts
    part1 = w[0]*Costheta - w[1]/2.*Costheta 
      + w[2]/2.*(v4Nu.E()+modkprime-(v4lep.E()+v4Nu.E())*Costheta)
      + w[3]/2.*(TMath::Power(Mlep,2)*Costheta+2.*v4lep.E()
		 *(v4lep.E()+modkprime)*TMath::Power(Sinh,2))
      - w[4]*(v4lep.E()+modkprime)/2.;
    part2 = 2.*w[0]*TMath::Power(Sinh,2) + w[1]*TMath::Power(Cosh,2)
      - w[2]*(v4lep.E()+v4Nu.E())*TMath::Power(Sinh,2)
      + TMath::Power(Mlep,2)/(v4lep.E()*(v4lep.E()+modkprime))*part1;
    xsec = modkprime*v4lep.E()*TMath::Power(Gfermi,2)*2./pi*part2*facconv;

    if(xsec >= 0.0){
      //nothing, but this will return false for negative values AND nan.
    } else {
      //  Sometimes Costheta is just over 1.0 due to fluke or numerical precision, so Sintheta would be undefined.
      /*
      std::cout << "RIK is bogus xsec ? " << xsec << " " << part1 << " " << part2
		<< " w[0] " << w[0] << " " << w[1] << " " << w[2] << " " << w[3] << " " << w[4]
		<< " wtotd[] " << wtotd[0] << " " << wtotd[1] << " " << wtotd[2] << " " << wtotd[3] << " " << wtotd[4] << " "
		<< " vec " << v4q.Vect().Mag() << " " << q0  << " " << v4q.Px() << " " << v4q.Py() << " " << v4q.Pz()
		<< " v4qX " << v4Nu.X() << " " << v4lep.X() << " " << Costheta << " " << Sintheta << " " << modkprime
		<< std::endl;
      */
      xsec = 0.;      
    }
  }

  else {//nothing transfered to nucleus thus no mec
    xsec = 0.;
  }
  
  return xsec * (1.0E-41 * units::cm2);
}

//___________________________________________________________________________
double MECLoadHadTensor::TotalXsecAtE(int targetpdg, int nupdg, double Enu){

  // this does a numerical integral in q0 q3 space.

  // at 3.0 GeV neutrino Carbon, this delivers 1467.16e-41 cm2 
  // This is good to 0.5% with the 1200, 0.001 step size integration.
  
  double maxstep = 1200;
  double stepsize = 0.001;
  long double totalXsec = 0.0;

  int lpdg;
  if(TMath::Abs(nupdg)==16)lpdg=15;
  else if(TMath::Abs(nupdg)==14)lpdg = 13;
  else if(TMath::Abs(nupdg)==12)lpdg = 11;
  else {
    //else not a neutrino?
    LOG("MEC", pERROR) << "Trying to get XSec for not a neutrino " << lpdg;
    assert(false);
  }
  
  double lmass = PDGLibrary::Instance()->Find( lpdg )->Mass();

  for(int iq3=0; iq3<=maxstep; iq3++){
    for(int iq0=0; iq0<=iq3; iq0++){
      double dq3 = stepsize * (double)iq3;
      double dq0 = stepsize * (double)iq0;

      double tmu = 0.0;
      double cost = 0.0;
      double area = 0.0;

      //int targetpdg, int nupdg, double Enu, double Tmu, double Costheta,  vector <BLI2DNonUnifGrid *> HadTensor
      GetTmuCostFromq0q3(dq0,dq3,Enu,lmass,tmu,cost,area);

      if(tmu < 0.0)continue;
      if(area <= 0.0)continue;
      if(cost < -1.0)continue;

      
      double Xsec = this->XSecFullAll(targetpdg, nupdg, Enu, tmu, cost);
      //std::cout << "RIK total " << "tgt " << targetpdg << " " << nupdg << " " << dq0 << " " << dq3 << " " << Enu << " " << lmass << " " << tmu << " " << cost << " " << area << " Xsec " << Xsec << " " <<  totalXsec << std::endl;

      totalXsec += (long double) ( Xsec * area  * stepsize * stepsize);
      
    }
  }

  return totalXsec; //  don't do this, its done by the Xsec call * (1.0E-41 * units::cm2);   // is in E-41 cm2.
  
}

double MECLoadHadTensor::GetTmuCostFromq0q3(double dq0, double dq3, double Enu, double lmass, double &tmu, double &cost, double &area){

  tmu = Enu - dq0 - lmass;
  if(tmu < 0.0){
    cost = -999;
    tmu = -999;
    area = 0.0;
    return -999;
  }

  double thisE = tmu + lmass;
  double thisp = sqrt( thisE * thisE - lmass * lmass);
  double numerator =  Enu * Enu + thisp * thisp - dq3 * dq3; 
  double denominator = 2.0 * thisp * Enu;
  if(denominator <= 0.0 ){
    cost = 0.0;
    if(denominator < 0.0){
      //std::cout << "Error in tmu cost from q0 q3 " << denominator << std::endl;
      return -999;
    }
  }
  else cost = numerator / denominator;

  if(TMath::Abs(numerator) > TMath::Abs(denominator)){
    //std::cout << "Error in tmu cost from q0 q3, bad cost " << cost << " " << numerator << " " << denominator << std::endl;
    cost = -999;
    tmu = -999;
    area = 0.0;
    return -999;
  }
  
  // xCrossSect is not yet in the right units for this particular case.
  // need areaElement to go dsigma/dTmudcost to dsigma/dq0dq3
  // Recompute the area element jacobian
  // dT/dq0 x dC/dq3 - dT/dq3 x dC/dq0
  double areaElement = 0.0;
  //double veryCloseToZero = 0.000001;  // in GeV, this should be safe.
  double insqrt = 0.0;
  numerator = 0.0;
  insqrt = Enu * Enu - 2.0 * Enu * dq0 + dq0 * dq0 - lmass * lmass;
  numerator = dq3 / Enu;
  if(insqrt < 0.0)areaElement=0.0;
  else areaElement = numerator / TMath::Sqrt(insqrt);
  area = areaElement;
  
  return 0;
  
}

//__________________________________________________________
// specific XSec
double MECLoadHadTensor::XSecFullAll(int targetpdg, int nupdg, double Enu, double Tmu, double CosTheta) {
  if (targetpdg == 1000060120){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorFullAll_C12_2DGrids );
    }
  }
  else if (targetpdg == 1000080160){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorFullAll_O16_2DGrids );
    }
  }
    else if (targetpdg == 1000200400){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorFullAll_Ca40_2DGrids );
    }
    }
  //else

  return -1;
}
double MECLoadHadTensor::XSecFullpn(int targetpdg, int nupdg, double Enu, double Tmu, double CosTheta) {
  if (targetpdg == 1000060120){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorFullpn_C12_2DGrids );
    }
  }
  else if (targetpdg == 1000080160){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorFullpn_O16_2DGrids );
    }
  }
  else if (targetpdg == 1000200400){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorFullpn_Ca40_2DGrids );
    }
  }
  //else
  return -1;
}
double MECLoadHadTensor::XSecDeltaAll(int targetpdg, int nupdg, double Enu, double Tmu, double CosTheta) {
  if (targetpdg == 1000060120){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorDeltaAll_C12_2DGrids );
    }
  }
  else if (targetpdg == 1000080160){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorDeltaAll_O16_2DGrids );
    }
  }
  else if (targetpdg == 1000200400){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorDeltaAll_Ca40_2DGrids );
    }
  }
  //else

  return -1;
}
double MECLoadHadTensor::XSecDeltapn(int targetpdg, int nupdg, double Enu, double Tmu, double CosTheta) {
  if (targetpdg == 1000060120){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorDeltapn_C12_2DGrids );
    }
  }
  else if (targetpdg == 1000080160){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorDeltapn_O16_2DGrids );
    }
  }
  else if (targetpdg == 1000200400){
    for (int i=0 ; i < 5; i++){
      return XSec( targetpdg, nupdg, Enu, Tmu, CosTheta, HadTensorDeltapn_Ca40_2DGrids );
    }
  }
  //else

  return -1;
}



//____________________________________________________________________
// write max xsec table
//  -- this should be made once and saved as a text file to read in
void MECLoadHadTensor::WriteMaxXSecTables(int targetpdg, int nupdg){

  // make file name:
  char filename [50];
  int n;

  ofstream maxxsecfile;
  vector <double> MaxXSecvect;

  // x->Tmu, y->CosTheta, z->E
  double xbins=1000;
  double ybins=1000;
  int zbins=92;

  double xmin = 0;
  double xmax = 5;
  double ymin = -1;
  double ymax = 1;

  TH2F *temphist = new TH2F ( "temphist", "temphist", xbins, xmin, xmax, ybins, ymin, ymax);

  double bEnu;
  double aTmu; 
  double aCostheta;
  double xsec;

  double aEnu[92] = {0.10  , 0.1025, 0.1050, 0.1075, 0.11  ,   
		  0.1125, 0.1150, 0.1175, 0.12  , 0.1225, 
		  0.1250, 0.1275, 0.13  , 0.1325, 0.1350, 
		  0.1375, 0.14  , 0.1425, 0.145 , 0.1475, 
		  0.15  , 0.1525, 0.155 , 0.1575, 0.16  ,
		  0.1625, 0.165 , 0.17  , 0.175 , 0.18  , 
		  0.185 , 0.19  , 0.2   , 0.205 , 0.210 , 
		  0.22  , 0.23  , 0.24  , 0.25  , 0.26  , 
		  0.275 , 0.3   , 0.325 , 0.35  , 0.375 , 
		  0.4   , 0.425 , 0.45  , 0.475 , 0.5   , 
		  0.525 , 0.55  , 0.6   , 0.65  , 0.7   , 
		  0.75  , 0.8   , 0.85  , 0.9   , 0.95  , 
		  1     , 1.05  , 1.1   , 1.15  , 1.2   , 
		  1.25  , 1.3   , 1.4   , 1.49  , 1.4999999,
		  1.5000001, 1.502 , 1.5075, 1.513 , 1.515 , 
		  1.52  , 1.525 , 1.55  , 1.575 , 1.6   , 
		  1.65  , 1.8   , 1.9   , 2     , 2.5   , 
		  3     , 4     , 6     , 7.5   , 10    , 
			20    , 30};


  for (int tables = 0; tables < 2; tables++){
    if (tables ==0) n=sprintf(filename, "/home/jackie/work/genie/GENIE_2_8/data/evgen/mectensor/nieves/MaxXSecAll%iv%i.txt", targetpdg, nupdg);
    else if (tables == 1) n=sprintf(filename, "/home/jackie/work/genie/GENIE_2_8/data/evgen/mectensor/nieves/MaxXSecDelta%iv%i.txt", targetpdg, nupdg);
    
    // only run this if the file doesn't already exist:
    if (FILE *file = fopen(filename, "r")) {
      fclose(file);
      std::cout << "max xsec file " << filename << " exists, moving right along" << std::endl;
    } 
    else {
      // output a text file
      
      std::cout << "max xsec file " << filename << " does not yet exist, so make one" << std::endl;
      
      // -- variables -- //
      
      maxxsecfile.open(filename);
      
      
      // -- make sure the vectors are empty before filling -- //
      MaxXSecvect.clear();
      Enuvect.clear();
      
      // -- loop over E, Tmu, Costheta and pull xsec -- //
      for (int z = 0; z< zbins; z++){
	bEnu = aEnu[z];
	
	for (int i = 0; i < xbins; i++){
	  aTmu = (xmin + ((xmax - xmin)/xbins)/2. )+ i * (xmax - xmin)/xbins ;
	  
	  for (int j = 0; j < ybins; j++){
	    aCostheta = (ymin + ((ymax - ymin)/ybins)/2. )+ j * (ymax - ymin)/ybins;
	    
	    // get xsec and store in 2d histogram
	    if (tables == 0) xsec =  XSecFullAll(targetpdg, nupdg, bEnu, aTmu, aCostheta);
	    else if (tables == 1) xsec =  XSecDeltaAll(targetpdg, nupdg, bEnu, aTmu, aCostheta);
	    
	    temphist->SetBinContent( i, j, xsec);
	  }
	}
	
	// fill xsec and enu vectors
	MaxXSecvect.push_back(temphist->GetBinContent(temphist->GetMaximumBin()));
	Enuvect.push_back(bEnu);
      }
      
      // write enuvect to file
      for (int i = 0; i < 92; i++){
	maxxsecfile << Enuvect[i] << " " << MaxXSecvect[i] << std::endl;
      }
      maxxsecfile.close();
    }

  }// end loop over tables
}

//____________________________________________________________________
// read files into arrays
void MECLoadHadTensor::ReadMaxXSecTables(int targetpdg,int nupdg)
{
  // ------ open and check file ---------- //


  // make file name
  char filename1 [50];
  char filename2 [50];
  int n1;
  int n2;
  if (nupdg > 0){
  n1=sprintf(filename1, "MaxXSecAll%iv%i.txt", targetpdg, nupdg);
  n2=sprintf(filename2, "MaxXSecDelta%iv%i.txt", targetpdg, nupdg);
  }
  else{
  n1=sprintf(filename1, "MaxXSecAll%iv%i.txt", targetpdg, nupdg);
  n2=sprintf(filename2, "MaxXSecDelta%iv%i.txt", targetpdg, nupdg);
  }

  // open file
  std::ifstream file1(filename1, ios::in);
  std::ifstream file2(filename2, ios::in);

  // check file exists
  if(!file1.good()){
    // log or output file read error //
    //std::cout << "bad file name: " << filename << std::endl;
    return ;
  }
  if(!file2.good()){
    // log or output file read error //
    //std::cout << "bad file name: " << filename << std::endl;
    return ;
  }

  

  double temp1, temp2, temp3, temp4;
  vector <double> tempXSecVect1;
  vector <double> tempXSecVect2;
  tempXSecVect1.clear();
  tempXSecVect2.clear();
  Enuvect.clear();

  // Read file and fill vectors
  for (int i=0; i<92; i++){
    file1 >> temp1 >> temp2;
    file2 >> temp3 >> temp4;

    tempXSecVect1.push_back(temp2);
    tempXSecVect2.push_back(temp4);
    Enuvect.push_back(temp1);
  }

  // now set correct vector, given target and neutrino pdg...
  if(nupdg > 0 ){
    if (nupdg == 12){
      if (targetpdg == 1000060120){
	MaxXSecAllvectC12v12 = tempXSecVect1;
	MaxXSecDeltavectC12v12 = tempXSecVect2;
      }
      else if (targetpdg == 1000080160)	{
	MaxXSecAllvectO16v12 = tempXSecVect1;
	MaxXSecDeltavectO16v12 = tempXSecVect2;
      }
      else if (targetpdg == 1000200400)	{
	MaxXSecAllvectCa40v12 = tempXSecVect1;
	MaxXSecDeltavectCa40v12 = tempXSecVect2;
      }
    }
    else if (nupdg == 14){
      if (targetpdg == 1000060120){
	MaxXSecAllvectC12v14 = tempXSecVect1;
	MaxXSecDeltavectC12v14 = tempXSecVect2;
      }
      else if (targetpdg == 1000080160)	{
	MaxXSecAllvectO16v14 = tempXSecVect1;
	MaxXSecDeltavectO16v14 = tempXSecVect2;
      }
      else if (targetpdg == 1000200400)	{
	MaxXSecAllvectCa40v14 = tempXSecVect1;
	MaxXSecDeltavectCa40v14 = tempXSecVect2;
      }
    } 
    else if (nupdg == 16){
      if (targetpdg == 1000060120) {
	MaxXSecAllvectC12v16 = tempXSecVect1;
	MaxXSecDeltavectC12v16 = tempXSecVect2;
      }
      else if (targetpdg == 1000080160)	{
	MaxXSecAllvectO16v16 = tempXSecVect1;
	MaxXSecDeltavectO16v16 = tempXSecVect2;
      }
      else if (targetpdg == 1000200400)	{
	MaxXSecAllvectCa40v16 = tempXSecVect1;
	MaxXSecDeltavectCa40v16 = tempXSecVect2;
      }
    } 
  }

  else{
    if (nupdg == -12){
      if (targetpdg == 1000060120){
	MaxXSecAllvectC12av12 = tempXSecVect1;
	MaxXSecDeltavectC12av12 = tempXSecVect2;
      }
      else if (targetpdg == 1000080160)	{
	MaxXSecAllvectO16av12 = tempXSecVect1;
	MaxXSecDeltavectO16av12 = tempXSecVect2;
      }
      else if (targetpdg == 1000200400)	{
	MaxXSecAllvectCa40av12 = tempXSecVect1;
	MaxXSecDeltavectCa40av12 = tempXSecVect2;
      }
    }
    else if (nupdg == -14){
      if (targetpdg == 1000060120){
	MaxXSecAllvectC12av14 = tempXSecVect1;
	MaxXSecDeltavectC12av14 = tempXSecVect2;
      }
      else if (targetpdg == 1000080160)	{
	MaxXSecAllvectO16av14 = tempXSecVect1;
	MaxXSecDeltavectO16av14 = tempXSecVect2;
      }
      else if (targetpdg == 1000200400)	{
	MaxXSecAllvectCa40av14 = tempXSecVect1;
	MaxXSecDeltavectCa40av14 = tempXSecVect2;
      }
    } 
    else if (nupdg == -16){
      if (targetpdg == 1000060120) {
	MaxXSecAllvectC12av16 = tempXSecVect1;
	MaxXSecDeltavectC12av16 = tempXSecVect2;
      }
      else if (targetpdg == 1000080160)	{
	MaxXSecAllvectO16av16 = tempXSecVect1;
	MaxXSecDeltavectO16av16 = tempXSecVect2;
      }
      else if (targetpdg == 1000200400)	{
	MaxXSecAllvectCa40av16 = tempXSecVect1;
	MaxXSecDeltavectCa40av16 = tempXSecVect2;
      }
    } 
  }

}


//________________________________________________________________
// get max xsec - interpolate/extrapolate.
double MECLoadHadTensor::MaxXSecAll(int targetpdg, int nupdg, double Enu)
{

  vector <double> MaxXSecvect;

  if(nupdg > 0 ){
    if (nupdg == 12){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecAllvectC12v12;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecAllvectO16v12;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecAllvectCa40v12;
    }
    else if (nupdg == 14){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecAllvectC12v14;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecAllvectO16v14;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecAllvectCa40v14;
    } 
    else if (nupdg == 16){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecAllvectC12v16;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecAllvectO16v16;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecAllvectCa40v16;
    } 
  }

  else if(nupdg < 0 ){
    if (nupdg == -12){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecAllvectC12av12;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecAllvectO16av12;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecAllvectCa40av12;
    }
    else if (nupdg == -14){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecAllvectC12av14;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecAllvectO16av14;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecAllvectCa40av14;
    } 
    else if (nupdg == -16){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecAllvectC12av16;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecAllvectO16av16;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecAllvectCa40av16;
    } 
  }
  return MaxXSec(targetpdg, nupdg, Enu, MaxXSecvect);
}



double MECLoadHadTensor::MaxXSecDelta(int targetpdg, int nupdg, double Enu)
{

  vector <double> MaxXSecvect;

  if(nupdg > 0 ){
    if (nupdg == 12){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecDeltavectC12v12;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecDeltavectO16v12;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecDeltavectCa40v12;
    }
    else if (nupdg == 14){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecDeltavectC12v14;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecDeltavectO16v14;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecDeltavectCa40v14;
    } 
    else if (nupdg == 16){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecDeltavectC12v16;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecDeltavectO16v16;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecDeltavectCa40v16;
    } 
  }

  else if(nupdg < 0 ){
    if (nupdg == -12){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecDeltavectC12av12;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecDeltavectO16av12;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecDeltavectCa40av12;
    }
    else if (nupdg == -14){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecDeltavectC12av14;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecDeltavectO16av14;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecDeltavectCa40av14;
    } 
    else if (nupdg == -16){
      if (targetpdg == 1000060120) MaxXSecvect = MaxXSecDeltavectC12av16;
      else if (targetpdg == 1000080160)	MaxXSecvect = MaxXSecDeltavectO16av16;
      else if (targetpdg == 1000200400)	MaxXSecvect = MaxXSecDeltavectCa40av16;
    } 
  }

  return MaxXSec(targetpdg, nupdg, Enu, MaxXSecvect);

}





  double MECLoadHadTensor::MaxXSec(int targetpdg, int nupdg, double Enu, vector<double> MaxXSecvect)
{

  // -- get the right vector bins -- //
  // get bins:
  int ebinhi=-1;
  int iterate = 0;
  int ebinlo;
  double elo;
  double xseclo;

  // select proper MaxXsecvect based on target and nu
  //  MaxXSecvect.clear();
  

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

