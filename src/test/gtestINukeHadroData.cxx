//____________________________________________________________________________
/*!

\program testINukeHadroData

 cannabalized to test my mec code ^_^;



*/
//____________________________________________________________________________

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



int main()
{
  std::cout << "Starting..." << std::endl;
  
  MECLoadHadTensor * tensorclassthingy = MECLoadHadTensor::Instance();
  
  const vector <BLI2DNonUnifGrid *> Grids2D = tensorclassthingy->HadTensorC12();
  
  int num = 2;

  double xmin=Grids2D[num]->XMin();
  double xmax=Grids2D[num]->XMax();
  double xstep = (xmax - xmin)/100.;
  double ymin=Grids2D[num]->YMin();
  double ymax=Grids2D[num]->YMax();
  double ystep = (ymax - ymin)/100.;

  TH2F *gridhist1 = new TH2F ( "gridhist1", "Tensor W[2]",20, ymin, ymax, 20, xmin, xmax);
  

  for (int i = 0 ; i < 20 ; i++ ){
    for (int j = 0 ; j < 20 ; j++) {

      gridhist1->SetBinContent(j,i,
			       Grids2D[2]->Evaluate( xmin+(i*xstep), 
						      ymin+(j*ystep) ));

    }
  }  

  // and now for the xsec. let's do it with Enu at 800 mev?
  double xbins=4;
  double ybins=4;
  xmin = 0;
  xmax = .8;
  ymin = -1;
  ymax = 1;

  TH2F *gridhist2 = new TH2F ( "gridhist2", "XSec (Enu = 800 MeV)", xbins, xmin, xmax, ybins, ymin, ymax);

  int pdgnu = 14;
  double Enu = 0.8;
  double Tmu, Costheta;
  double xsec;

  for (int i = 0; i < xbins; i++){
    for (int j = 0; j < ybins; j++){
      
      Tmu = (xmin + ((xmax - xmin)/xbins)/2. )+ i * (xmax - xmin)/xbins ;
      Costheta = (ymin + ((ymax - ymin)/ybins)/2. )+ j * (ymax - ymin)/ybins;
      xsec =  tensorclassthingy->XSec(pdgnu, Enu, Tmu, Costheta);
      std::cout << " ~~~ " << "e: " << Enu << ", Tmu: " << Tmu << ", costheta: " << Costheta << ", xsec: " << tensorclassthingy->XSec(pdgnu, Enu, Tmu, Costheta) << std::endl;
      gridhist2->SetBinContent( i, j, xsec);

    }}

  std::cout << "max from hist: " << gridhist2->GetBinContent(gridhist2->GetMaximumBin()) << std::endl;
  std::cout << "max from fn: " << tensorclassthingy->MaxXSec(.68) << std::endl;
  std::cout << "** test in gtest: " << tensorclassthingy->XSec(14, 0.8, 0.5, .75) << std::endl;

  TFile f("tensor.root","recreate");
  gridhist1->Write();
  gridhist2->Write();
  f.Close();


  /*  
  ////////
  MECLoadXSecFiles * mecfiles = MECLoadXSecFiles::Instance();
  //mecfiles->testMECNievesLoadXSecFiles(mecfiles->Nieves14C12(), 61,"shiggles");
  std:: cout << "Size of array, my fav number ... " << mecfiles->Nieves14C12().size() << endl;
  std:: cout << "Evaluate? ... " << mecfiles->Nieves14C12()[61]->Evaluate(.5,.5) << endl;

  std:: cout << "Max Xsec ... " << mecfiles->Nieves14C12()[61]->ZMax() << endl;
  std:: cout << "Min Xsec ... " << mecfiles->Nieves14C12()[61]->ZMin() << endl;

  std:: cout << "number? " << mecfiles->number() << endl;

  //std:: cout << "Evaluate? ... " << 
  
  */

  std::cout << "The End" << std::endl;
  return 0;
}
//____________________________________________________________________________
