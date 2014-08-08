//________________________________________________________
/*

generator for making NIEVES mec events

working version 
2014/07/10

*/
//_________________________________________________________

// -- includes -- //

#include <TMath.h>

#include "Base/XSecAlgorithmI.h"
#include "Conventions/Constants.h"
#include "Conventions/Controls.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGCore/RunningThreadInfo.h"
#include "EVGCore/EventGeneratorI.h"
#include "GHEP/GHepStatus.h"
#include "GHEP/GHepFlags.h"
#include "GHEP/GHepParticle.h"
#include "GHEP/GHepRecord.h"
#include "Messenger/Messenger.h"
#include "MEC/MECGenerator.h"
#include "Numerical/RandomGen.h"
#include "Nuclear/NuclearModelI.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"
#include "Utils/PrintUtils.h"
#include "MEC/MECLoadXSecFiles.h"
#include <iostream>

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;


//___________________________________________________________________________
MECGenerator::MECGenerator() :
EventRecordVisitorI("genie::MECGenerator")
{

}
//___________________________________________________________________________
MECGenerator::MECGenerator(string config) :
EventRecordVisitorI("genie::MECNievesGenerator", config)
{

}
//___________________________________________________________________________
MECGenerator::~MECGenerator()
{

}
//___________________________________________________________________________
void MECGenerator::ProcessEventRecord(GHepRecord * event) const
{

  // Start with Lepton //
  this -> SelectLeptonKinematics (event);
  this -> AddTargetRemnant       (event);
  this -> GenerateInitialHadrons (event);
  this -> RecoilNucleonCluster   (event);
  this -> DecayNucleonCluster    (event);

  /*
  this -> AddTargetRemnant      (event); /// shortly, this will be handled by the InitialStateAppender module
  this -> GenerateFermiMomentum (event);
  this -> SelectKinematics      (event);
  this -> AddFinalStateLepton   (event);

  this -> RecoilNucleonCluster  (event);
  this -> DecayNucleonCluster   (event);
  */

}
//___________________________________________________________________________
void MECGenerator::SelectLeptonKinematics (GHepRecord * event) const
{

  // -- Constants --------------------------------- //
  double Q3Max_data = 1.5;
  double Q3Max = 1.3;
  
  // Nieves' "Q" value - An energy related to the difference between
  //                     initial and final nuclear states

  double Q = 0.016827; // for now for carbon for neutrinos
 
 
  // -- Event Properties -----------------------------//
  Interaction * interaction = event->Summary();
  double Enu = interaction->InitState().ProbeE(kRfHitNucRest);
  double LepMass = interaction->FSPrimLepton()->Mass();

  // -- Lepton Kinematic Limits ----------------------------------------- //
 
  double Costh; // lepton angle
  double CosthMax=1.;
  double CosthMin;

  double T;  // lepton kinetic energy
  double TMax;
  double TMin;

  double Plep; // lepton 3 momentum
  double Elep; // lepton energy
  
  double Q0; // energy component of q four vector
  double Q3; // magnitude of transfered 3 momentum
  double Q2; // properly Q^2 (Q squared) - transfered 4 momentum.


  // -- load xsec tables --- //
  MECLoadXSecFiles * xsecfiles = MECLoadXSecFiles::Instance();



  TMax = Enu - LepMass - Q;
  //TMax = Enu - LepMass;
  
  if(Enu < Q3Max){
    TMin = 0 ;
    CosthMin = -1 ; 
  }
  else{
    TMin = TMath::Sqrt( TMath::Power(LepMass,2) + TMath::Power((Enu+Q3Max),2) ) - LepMass ;
    CosthMin = TMath::Sqrt( 1 - TMath::Power(( Q3Max / Enu ),2) ) ;
    
  }
  
  // -- Generate and Test the Kinematics----------------------------------//

  RandomGen * rnd = RandomGen::Instance();
  bool accept = false;
  unsigned int iter = 0;

  // loop over different (randomly) selected T and Costh
  while (!accept) {
    iter++;
    if(iter > kRjMaxIterations) {
      // error if try too many times
      LOG("MEC", pWARN)
           << "Couldn't select a valid W, Q^2 pair after " 
           << iter << " iterations";
        event->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
    }
    
    // generate random T and Costh
    T = TMin + (TMax-TMin)*rnd->RndKine().Rndm();
    Costh = CosthMin + (CosthMax-CosthMin)*rnd->RndKine().Rndm();
  
    // Calc Useful Values 
    Plep = TMath::Sqrt( T * (T + (2. * LepMass)));
    Q3 = TMath::Sqrt( TMath::Power(Plep,2) + TMath::Power(Enu,2) - (2. * Plep * Enu * Costh));
      
    // Check if allowed 3 momentum transfer Q3

    if (Q3 < Q3Max){

      // Accept/Reject

      // get max xsec
      double XSecMax = xsecfiles->MaxXSec(Enu);//(target,neuflavor)

      // get xsec(e)
      double XSec = xsecfiles->XSec(Enu, Costh, T);//(target,neuflavor)

      /*
      // get e index (and eindex+1)
      int  eindlow = xsecfiles->EtoIndex(Enu);
      int  eindhigh = eindlow+1;

      // get max xsec for given energy
      double XSecMaxlow = xsecfiles->Nieves14C12()[eindlow]->ZMax();
      double XSecMaxhigh = xsecfiles->Nieves14C12()[eindhigh]->ZMax();
      

      //accept = xsecfiles->Nieves14C12()[eindex]->Evaluate(T,Costh) > XSecMax*rnd->RndKine().Rndm();  
      accept = xsecfiles->Nieves14C12()[eindex]->Evaluate(Costh,T) > XSecMax*rnd->RndKine().Rndm();  
*/      
      accept = XSec > XSecMax*rnd->RndKine().Rndm();


    }// end if q3 test
  }// end while
  
  // -- finish lepton kinematics
  
  // define coordinate system wrt neutrino: z along neutrino, xy perp
  
  // Cos theta gives us z, the rest in xy:
  double PlepZ = Plep * Costh;
  double PlepXY = Plep * TMath::Sqrt(1. - TMath::Power(Costh,2));

  // random rotation about unit vector for phi direction
  double phi= 2 * kPi * rnd->RndLep().Rndm();
  // now fill x and y from PlepXY
  double PlepX = PlepXY * TMath::Cos(phi);
  double PlepY = PlepXY * TMath::Sin(phi);
  
  // Rotate lepton momentum vector from the reference frame (x'y'z') where 
  // {z':(neutrino direction), z'x':(theta plane)} to the LAB
  TVector3 unit_nudir = event->Probe()->P4()->Vect().Unit();
  TVector3 p3l(PlepX, PlepY, PlepZ);
  p3l.RotateUz(unit_nudir);

  // Lepton 4-momentum in LAB
  Elep = TMath::Sqrt(TMath::Power(LepMass,2)+TMath::Power(PlepX,2)+TMath::Power(PlepY,2)+TMath::Power(PlepZ,2));
  TLorentzVector p4l(p3l,Elep);

  // Figure out the final-state primary lepton PDG code
  int pdgc = interaction->FSPrimLepton()->PdgCode();

  // Lepton 4-position (= interacton vtx)
  TLorentzVector v4(*event->Probe()->X4());

  int momidx = event->ProbePosition();

  // -- Store Values ------------------------------------------//

  // -- Interaction: Q2
  Q0 = Enu - Elep;
  Q2 = TMath::Power(Q3,2) - TMath::Power(Q0,2);

  interaction->KinePtr()->SetQ2(Q2, true);

  // -- Lepton
  event->AddParticle( pdgc, kIStStableFinalState, momidx, -1, -1, -1, p4l, v4);

}
//____________________________________________________________________

//___________________________________________________________________________
// coppied verbatum //
void MECGenerator::AddTargetRemnant(GHepRecord * event) const
{
// Add the remnant nucleus (= initial nucleus - nucleon cluster) in the
// event record.

  GHepParticle * target  = event->TargetNucleus();
  GHepParticle * cluster = event->HitNucleon();

  int Z = target->Z();
  int A = target->A();

  if(cluster->Pdg() == kPdgClusterNN) { A-=2; ;     }
  if(cluster->Pdg() == kPdgClusterNP) { A-=2; Z-=1; }
  if(cluster->Pdg() == kPdgClusterPP) { A-=2; Z-=2; }

  int ipdgc = pdg::IonPdgCode(A, Z);

  const TLorentzVector & p4cluster = *(cluster->P4());
  const TLorentzVector & p4tgt     = *(target ->P4());

  const TLorentzVector p4 = p4tgt - p4cluster;
  const TLorentzVector v4(0.,0.,0., 0.);

  int momidx = event->TargetNucleusPosition();
  event->AddParticle(ipdgc,kIStStableFinalState, momidx,-1,-1,-1, p4,v4);  

  std::cout << "*** cluster p = " << p4.E() << std::endl; 

}
//_________________________________________________________________________
void MECGenerator::GenerateInitialHadrons  (GHepRecord * event) const
{

  // -- Inputs: Q4, Fermi Momentum --------------------------
  // get neutrino & its 4-momentum
  GHepParticle * neutrino = event->Probe();
  assert(neutrino);
  TLorentzVector p4nu(*neutrino->P4());

  // get final state primary lepton & its 4-momentum
  GHepParticle * fsl = event->FinalStatePrimaryLepton();
  assert(fsl);
  TLorentzVector p4l(*fsl->P4());

  // calculate 4-momentum transfer
  TLorentzVector Q4 = p4nu - p4l;

  std::cout << "*** Q0  =" << Q4.E() << std::endl; 

  GHepParticle * target_nucleus = event->TargetNucleus();
  assert(target_nucleus);
  GHepParticle * nucleon_cluster = event->HitNucleon();
  assert(nucleon_cluster);
  GHepParticle * remnant_nucleus = event->RemnantNucleus();
  assert(remnant_nucleus);

  // -- get nucleons from fermi gas
  Target tgt(target_nucleus->Pdg());
  PDGCodeList pdgv = this->NucleonClusterConstituents(nucleon_cluster->Pdg());
  assert(pdgv.size()==2);

  bool accept = false;
  TVector3 p31i;
  TVector3 p32i;
  unsigned int iter = 0;

  double M1i;
  double M2i;



  while(!accept){
    iter++;
    if(iter > kRjMaxIterations) {
      // error if try too many times
      LOG("MEC", pWARN)
           << "Couldn't select a valid W, Q^2 pair after " 
           << iter << " iterations";
        event->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
    }


    // generate nucleons
    tgt.SetHitNucPdg(pdgv[0]);
    fNuclModel->GenerateNucleon(tgt);
    p31i = fNuclModel->Momentum3();
    tgt.SetHitNucPdg(pdgv[1]);
    fNuclModel->GenerateNucleon(tgt);
    p32i = fNuclModel->Momentum3();

    // get properties
    // M1i, M2i, P1i, P2i, EBind1, EBind2
    M1i = PDGLibrary::Instance()->Find(pdgv[0])-> Mass(); // nucleon 1 mass
    M2i = PDGLibrary::Instance()->Find(pdgv[1])-> Mass(); // nucleon 2 mass
    double M12f;

  // m1f + m2f - for neutrino cases only (2112n, 2212p)
    if (pdgv[0] == pdgv[1]) M12f =  
			    PDGLibrary::Instance()->Find(2112)->Mass() +
			    PDGLibrary::Instance()->Find(2212)->Mass() ;
    else M12f = PDGLibrary::Instance()->Find(2212)->Mass() * 2.;



    double EBind1 = -.0250; // switching to gev... (mev ... ?)
    double EBind2 = -.0250;

    // Check if valid initial state.

    // will these initial states allow for possible final states?
    //  why wouldn't they? 
    //    - there is 'lost energy'
    //      - binding energy
    //      - mass change     
    //    - there is minimum momentum to escape
    //      - final mom > fermi momentum (for global, fermi = to initial mom?)
    //
    //  ** these are rough cuts - if it isn't "allowed" it will be
    //  ** identified later, these can just speed the code up by not
    //  ** allowing it to get that far if things are "obviously"
    //  ** unphysical.
    //
    // If use local fermi gas, then position matters. First check that position can give reasonable results
    // if position doesn't matter, just pick hadrons and see if they will work (going with this one...)?

    double P1i = p31i.Mag();
    double P2i = p32i.Mag();

    // calculate energy of hadronic system
    double W0 = Q4(3) + 
      TMath::Sqrt( TMath::Power(M1i,2) +  TMath::Power(P1i,2) ) + 
      TMath::Sqrt( TMath::Power(M2i,2) +  TMath::Power(P2i,2) ) + 
      EBind1 + EBind2;
     
    double W3min = TMath::Sqrt(TMath::Power(Q4(0),2) + TMath::Power(Q4(1),2) 
		     + TMath::Power(Q4(2),2)) - P1i - P2i;
    
    std::cout << "*** P1i, P2i = " << P1i << ", " << P2i << std::endl;

    if (W3min < 0) W3min = 0;
    
    double invMassMax = TMath::Sqrt(TMath::Power(W0,2) - TMath::Power(W3min,2));
    if (invMassMax < (M12f)){
      accept = false; // FAIL
    }

    //else if ( W0 < emin){  // WTF is EMIN?!?!
    //  accept = false; // FAIL
    //}
      
    else { // PASS
      
      // Pick random directions? psh, already got those!
      // calculate momentum of hadronic system!
      double W3x = Q4(0)+ p31i(0) + p32i(0);
      double W3y = Q4(1)+p31i(1)+p32i(1);
      double W3z = Q4(2)+p31i(2)+p32i(2);
      double W32 = TMath::Power(W3x,2) + TMath::Power(W3y,2) + TMath::Power(W3z,2);
      double invmass;

      if ((TMath::Power(W0,2) - W32) < 0) invmass = -1;
      else invmass = TMath::Sqrt((W0*W0)-W32);
      
      if (invmass > M12f) accept = true;
      else accept = false;
  
    }
    
  }//end while(!accept)
  // Now Physically allowed


  // set initial hadron and cluster states
  
  
  double Mi  = PDGLibrary::Instance()->Find(target_nucleus->Pdg() )-> Mass(); // initial nucleus mass
  double M2n = PDGLibrary::Instance()->Find(nucleon_cluster->Pdg())-> Mass(); // nucleon cluster mass
  
  // calcute nucleon cluster momentum
  TVector3 p3 = p31i + p32i;
    
  // nucleon cluster energy

  double EN = TMath::Sqrt(p3.Mag2() + M2n*M2n);

  // set the remnant nucleus and nucleon cluster 4-momenta at GHEP record

  TLorentzVector p4nclust   (   p3.Px(),    p3.Py(),    p3.Pz(),  EN   );
  TLorentzVector p4remnant   (-1*p3.Px(), -1*p3.Py(), -1*p3.Pz(), Mi-EN);
       

  std::cout << "*** cluster Pi = " << p4nclust.E() << std::endl;

  nucleon_cluster->SetMomentum(p4nclust);
  remnant_nucleus->SetMomentum(p4remnant);

  // set the nucleon cluster 4-momentum at the interaction summary 

  event->Summary()->InitStatePtr()->TgtPtr()->SetHitNucP4(p4nclust);
  //event->Summary()->InitStatePtr()->TgtPtr()

  
  TLorentzVector p41i(p31i,M1i);
  TLorentzVector p42i(p32i,M2i);

  TLorentzVector v41(* nucleon_cluster->X4());

  event->AddParticle( pdgv[0],kIStInitialState, 42, -1, -1, -1, p41i, v41);
  event->AddParticle( pdgv[1],kIStInitialState, 43, -1, -1, -1, p42i, v41);
  
}

//_________________________________________________________________________
void MECGenerator::RecoilNucleonCluster    (GHepRecord * event) const
{

   // get di-nucleon cluster & its 4-momentum
  GHepParticle * nucleon_cluster = event->HitNucleon();
  assert(nucleon_cluster);
  TLorentzVector p4cluster(*nucleon_cluster->GetP4());

  // get neutrino & its 4-momentum
  GHepParticle * neutrino = event->Probe();
  assert(neutrino);
  TLorentzVector p4v(*neutrino->P4());

  // get final state primary lepton & its 4-momentum
  GHepParticle * fsl = event->FinalStatePrimaryLepton();
  assert(fsl);
  TLorentzVector p4l(*fsl->P4());

  // calculate 4-momentum transfer
  TLorentzVector q = p4v - p4l;

  // calculate binding energy // for now just value from neut
  double ebind = 0.0250; // i "think" it's gev...

  // calculate recoil nucleon cluster 4-momentum
  TLorentzVector p4cluster_recoil = p4cluster + q - ebind*2.;

  // take off binding energy
  p4cluster_recoil.SetE(p4cluster_recoil.E());

  // work-out recoil nucleon cluster code
  LOG("MEC", pINFO) << "Interaction summary";
  LOG("MEC", pINFO) << *event->Summary();
  int recoil_nucleon_cluster_pdg = event->Summary()->RecoilNucleonPdg();

  // vtx position in nucleus coord system
  TLorentzVector v4(*neutrino->X4());

  std::cout << "*** cluster pf = " << p4cluster_recoil.E() << std::endl;

  // add to the event record
  event->AddParticle(
    recoil_nucleon_cluster_pdg, kIStDecayedState, 
    2, -1, -1, -1, p4cluster_recoil, v4); 


}
//_________________________________________________________________________
void MECGenerator::DecayNucleonCluster  (GHepRecord * event) const
{

// Perform a phase-space decay of the nucleon cluster and add its decay
// products in the event record
//


  // get di-nucleon cluster
  int nucleon_cluster_id = 5;
  GHepParticle * nucleon_cluster = event->Particle(nucleon_cluster_id);
  assert(nucleon_cluster);

  // get decay products
  PDGCodeList pdgv = this->NucleonClusterConstituents(nucleon_cluster->Pdg());
  LOG("MEC", pINFO) << "Decay product IDs: " << pdgv;

  // Get the decay product masses
  vector<int>::const_iterator pdg_iter;
  int i = 0;
  double * mass = new double[pdgv.size()];
  double   sum  = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
    int pdgc = *pdg_iter;
    double m = PDGLibrary::Instance()->Find(pdgc)->Mass();
    mass[i++] = m;
    sum += m;
  }

  TLorentzVector * p4d = nucleon_cluster->GetP4();
  TLorentzVector * v4d = nucleon_cluster->GetX4();

  // Set the decay
  bool permitted = fPhaseSpaceGenerator.SetDecay(*p4d, pdgv.size(), mass);
  if(!permitted) {

     // clean-up 
     delete [] mass;
     delete p4d;
     delete v4d; 
     // throw exception
     genie::exceptions::EVGThreadException exception;
     exception.SetReason("Decay not permitted kinematically");
     exception.SwitchOnFastForward();
     throw exception;
  }

  // Get the maximum weight
  double wmax = -1;
  for(int i=0; i<200; i++) {
     double w = fPhaseSpaceGenerator.Generate();   
     wmax = TMath::Max(wmax,w);
  }
  assert(wmax>0);
  wmax *= 2;

  RandomGen * rnd = RandomGen::Instance();
  bool accept_decay=false;
  unsigned int itry=0;
  while(!accept_decay) 
  {
     itry++;

     if(itry > controls::kMaxUnweightDecayIterations) {

       // clean up
       delete [] mass;
       delete p4d;
       delete v4d;
       // throw exception
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select decay after N attempts");
       exception.SwitchOnFastForward();
       throw exception;
     }
     double w  = fPhaseSpaceGenerator.Generate();   
     if(w > wmax) {
 
     }
     double gw = wmax * rnd->RndDec().Rndm();
     accept_decay = (gw<=w);

  } //!accept_decay

  // Insert the decay products in the event record
  TLorentzVector v4(*v4d); 
  GHepStatus_t ist = kIStHadronInTheNucleus;
  int idp = 0;
  for(pdg_iter = pdgv.begin(); pdg_iter != pdgv.end(); ++pdg_iter) {
     int pdgc = *pdg_iter;
     TLorentzVector * p4fin = fPhaseSpaceGenerator.GetDecay(idp);
     std::cout <<"*** part pf = " << p4fin->P() << std::endl;
     event->AddParticle(pdgc, ist, nucleon_cluster_id,-1,-1,-1, *p4fin, v4);
     idp++;
  }

  // Clean-up
  delete [] mass;
  delete p4d;
  delete v4d;
}
//_________________________________________________________________________
PDGCodeList MECGenerator::NucleonClusterConstituents(int pdgc) const
{
  bool allowdup = true;
  PDGCodeList pdgv(allowdup);

  if(pdgc == kPdgClusterNN) { 
     pdgv.push_back(kPdgNeutron);
     pdgv.push_back(kPdgNeutron);
  }
  else
  if(pdgc == kPdgClusterNP) { 
     pdgv.push_back(kPdgNeutron);
     pdgv.push_back(kPdgProton);
  }
  else
  if(pdgc == kPdgClusterPP) { 
     pdgv.push_back(kPdgProton);
     pdgv.push_back(kPdgProton);
  }
  else 
  {
     LOG("MEC", pERROR) 
        << "Unknown di-nucleon cluster PDG code (" << pdgc << ")";
  }
 
  return pdgv;
}
//___________________________________________________________________________
void MECGenerator::Configure(const Registry & config)   
{
  Algorithm::Configure(config);
  this->LoadConfig();
} 
//___________________________________________________________________________ 
void MECGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void MECGenerator::LoadConfig(void)
{
  fNuclModel = 0;
      
  RgKey nuclkey = "NuclearModel";
  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);
}
//___________________________________________________________________________
















