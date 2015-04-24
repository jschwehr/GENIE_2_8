//________________________________________________________
/*

generator for making NIEVES mec events
 - using hadron tensor look up tables

working version 
2014/09/15 ~ J.Schwehr
2015/01/14 all neutrino flavors, oxygen and carbon ~J.Schwehr
2015/02/06 fixed hadron system, added spline generation code ~R. Gran

** hardcoded maxxsec for debugging

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
#include "MECTensor/MECTensorGenerator.h"
#include "Numerical/RandomGen.h"
#include "Numerical/BLI2D.h"
#include "Nuclear/NuclearModelI.h"
#include "PDG/PDGCodes.h"
#include "PDG/PDGUtils.h"
#include "PDG/PDGLibrary.h"
#include "Utils/KineUtils.h"
#include "Utils/PrintUtils.h"
#include "MECTensor/MECLoadHadTensor.h"
#include <iostream>

using namespace genie;
using namespace genie::utils;
using namespace genie::constants;
using namespace genie::controls;


//___________________________________________________________________________
MECTensorGenerator::MECTensorGenerator() :
EventRecordVisitorI("genie::MECTensorGenerator")
{

}
//___________________________________________________________________________
MECTensorGenerator::MECTensorGenerator(string config) :
EventRecordVisitorI("genie::MECTensorGenerator", config)
{

}
//___________________________________________________________________________
MECTensorGenerator::~MECTensorGenerator()
{

}
//___________________________________________________________________________
void MECTensorGenerator::ProcessEventRecord(GHepRecord * event) const
{

  this -> SelectLeptonKinematics (event);
  this -> AddTargetRemnant       (event);
  this -> GenerateInitialHadrons (event);
  //this -> RecoilNucleonCluster   (event);
  this -> DecayNucleonCluster    (event);

}
//___________________________________________________________________________
void MECTensorGenerator::SelectLeptonKinematics (GHepRecord * event) const
{
  
  int delta = 1; // 1: all, 2: no delta, 3: only delta

  // -- Constants --------------------------------- //
  double Q3Max = 1.2;
  
    // -- Event Properties -----------------------------//
  Interaction * interaction = event->Summary();
  InitialState * init_state = interaction->InitStatePtr();
  double Enu = interaction->InitState().ProbeE(kRfHitNucRest);
  double LepMass = interaction->FSPrimLepton()->Mass();
  int NuPDG = interaction->InitState().ProbePdg();
  int TgtPDG = interaction->InitState().TgtPdg();
  // interacton vtx
  TLorentzVector v4(*event->Probe()->X4());
  TLorentzVector tempp4(0,0,0,0);
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
  MECLoadHadTensor * hadtensor = MECLoadHadTensor::Instance(TgtPDG, NuPDG);

  
  // Uncomment this hack to generate a spline, waiting to hook the call to the real spline generation scheme.
  // Rik for testing, does this give the total cross section ?
  //std::cout << "RIK total " << hadtensor->TotalXsecAtE(TgtPDG, NuPDG, 3.0) << std::endl;
  /*
  double energies[61] = { 0.01000, 0.012115, 0.014678, 0.017783, 0.021544, 0.026102, 
			  0.031623, 0.038312, 0.046416, 0.056234, 0.068129, 0.082540,
			  0.1000, 0.12115, 0.14678, 0.17783, 0.21544, 0.26102, 
			  0.31623, 0.38312, 0.46416, 0.56234, 0.68129, 0.82540,
			  1.000, 1.2115, 1.4678, 1.7783, 2.1544, 2.6102, 
			  3.1623, 3.8312, 4.6416, 5.6234, 6.8129, 8.2540,
			  10.00, 12.115, 14.678, 17.783, 21.544, 26.102, 
			  31.623, 38.312, 46.416, 56.234, 68.129, 82.540,
			  100.0, 121.15, 146.78, 177.83, 215.44, 261.02, 
			  316.23, 383.12, 464.16, 562.34, 681.29, 825.40,
			  1000.0 };
			  
  for(int iii=0; iii<61; iii++){
    std::cout << "<knot> <E> " << energies[iii] << " </E> <xsec> " << hadtensor->TotalXsecAtE(TgtPDG, NuPDG, energies[iii]) << " </xsec> </knot> " << std::endl;
    //std::cout << "<knot> <E> " << energies[iii] << " </E> <xsec> " << hadtensor->TotalXsecAtE(TgtPDG, NuPDG, energies[iii]) << " </xsec> </knot> " << std::endl;

  } 
  
  */  
  

  TMax = Enu - LepMass ;
  
  if(Enu < Q3Max){
    TMin = 0 ;
    CosthMin = -1 ; 
  }
  else{
    TMin = TMath::Sqrt( TMath::Power(LepMass,2) + TMath::Power((Enu-Q3Max),2) ) - LepMass ;
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
           << "Couldn't select a valid Tmu, CosTheta pair after " 
           << iter << " iterations";
        event->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select lepton kinematics");
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

      // get max xsec and xsec(e) (set hit nucleon)
      
      if (delta == 1){ // all
	double XSecMax = 5.2e-11; //hadtensor->MaxXSecAll(TgtPDG, NuPDG, Enu);
	std::cout << "~*~ T, Costh: " << T << ", " << Costh << std::endl;
	double XSec = hadtensor->XSecFullAll(TgtPDG, NuPDG, Enu, T, Costh);
	accept = XSec > XSecMax*rnd->RndKine().Rndm();
	std::cout << "~~~ Xsec, Max, Accept: " << XSec << ", " << XSecMax << ", " << accept << std::endl; 
	// choose initial nucleons
	//RIK  HERE is where I would get two more Xsec from the Delta component
	//and prepare to tag the initial or final nucleon state.
	// determine if delta!
	double XSecDelta = hadtensor->XSecDeltaAll(TgtPDG, NuPDG, Enu, T, Costh);
	bool isPDD = rnd->RndKine().Rndm() > XSecDelta/XSec; //
      
	if(accept){
	  double XSecPN= hadtensor->XSecFullpn(TgtPDG, NuPDG, Enu, T, Costh);
	  double myrand = rnd->RndKine().Rndm();
	  double pnFraction = XSecPN / XSec;
	  std::cout << "RIK test for pn " << XSecPN << " " << XSec << " " << pnFraction << " " << myrand << std::endl;
	  if ( myrand <= pnFraction){ // rnd->RndKine().Rndm() <= XSecPN / XSec ){
	    // add cluster to event record
	    event->AddParticle(kPdgClusterNP, kIStNucleonTarget, 1, -1, -1, -1, tempp4, v4);
	    init_state->TgtPtr()->SetHitNucPdg(kPdgClusterNP);
	    //event->HitNucleon()->SetPdgCode(kPdgClusterNP);
	  }
	  else {
	    if (NuPDG > 0) {
	      event->AddParticle(kPdgClusterNN, kIStNucleonTarget, 1, -1, -1, -1, tempp4, v4);
	      init_state->TgtPtr()->SetHitNucPdg(kPdgClusterNN);
	      //event->HitNucleon()->SetPdgCode(kPdgClusterNN);
	    }
	    else {
	      event->AddParticle(kPdgClusterPP, kIStNucleonTarget, 1, -1, -1, -1, tempp4, v4);
	      init_state->TgtPtr()->SetHitNucPdg(kPdgClusterPP); 
	      //event->HitNucleon()->SetPdgCode(kPdgClusterPP);
	    }
	  }
	  // set scattering type to pdd
	  if (isPDD){
	    // set scattering type... somehow...
	  }
	  
	} // end if accept
      }// end if delta ==1
      
      /*
      if (delta == 2){ // no delta
	double XSecMaxDelta = hadtensor->MaxXSecDelta(TgtPDG, NuPDG, Enu);
	double XSecMaxAll = hadtensor->MaxXSecAll(TgtPDG, NuPDG, Enu);
	double XSecMaxNoDelta = XSecMaxAll - XSecMaxDelta;

	std::cout << "~*~ T, Costh: " << T << ", " << Costh << std::endl;
	double XSec = hadtensor->XSecDeltaAll(TgtPDG, NuPDG, Enu, T, Costh);
	accept = XSec > XSecMax*rnd->RndKine().Rndm();
	std::cout << "~~~ Xsec, Max, Accept: " << XSec << ", " << XSecMax << ", " << accept << std::endl; 
	
	// choose initial nucleons
	// determine hit nucleon based on pn / all xsec
	double XSecPN= hadtensor->XSecDeltapn(TgtPDG, NuPDG, Enu, T, Costh);
	if ( rnd->RndKine().Rndm() <= XSecPN / XSec ){
	  event->HitNucleon()->SetPdgCode(kPdgClusterNP);
	}
	else {
	  if (NuPDG > 0) event->HitNucleon()->SetPdgCode(kPdgClusterNN);
	  else event->HitNucleon()->SetPdgCode(kPdgClusterPP);
	}
      }
      */

      if (delta == 3){ // only delta
	double XSecMax = 5.2e-11;//hadtensor->MaxXSecDelta(TgtPDG, NuPDG, Enu);
	std::cout << "~*~ T, Costh: " << T << ", " << Costh << std::endl;
	double XSec = hadtensor->XSecDeltaAll(TgtPDG, NuPDG, Enu, T, Costh);
	accept = XSec > XSecMax*rnd->RndKine().Rndm();
	std::cout << "~~~ Xsec, Max, Accept: " << XSec << ", " << XSecMax << ", " << accept << std::endl; 
	
	// choose initial nucleons
	// determine hit nucleon based on pn / all xsec
	double XSecPN= hadtensor->XSecDeltapn(TgtPDG, NuPDG, Enu, T, Costh);
	if ( rnd->RndKine().Rndm() <= XSecPN / XSec ){
	  event->HitNucleon()->SetPdgCode(kPdgClusterNP);
	}
	else {
	  if (NuPDG > 0) event->HitNucleon()->SetPdgCode(kPdgClusterNN);
	  else event->HitNucleon()->SetPdgCode(kPdgClusterPP);
	}
      }

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
  //TLorentzVector v4(*event->Probe()->X4()); -- moved to start of code

  int momidx = event->ProbePosition();

  // -- SANITY CHECK -- //
  if (Elep > Enu) std::cout << "~*~*~ Energy Problem! Elep > Enu " << std::endl;
 

  // -- Store Values ------------------------------------------//

  // -- Interaction: Q2
  Q0 = Enu - Elep;
  Q2 = TMath::Power(Q3,2) - TMath::Power(Q0,2);

  interaction->KinePtr()->SetQ2(Q2, true);

  std::cout << "RIK q0 q3 q2 " << Q0 << " " << Q3 << " " << Q2 << " " << std::endl;

  // -- Lepton
  event->AddParticle( pdgc, kIStStableFinalState, momidx, -1, -1, -1, p4l, v4);

  // 
  std::cout << "~~~ LEPTON DONE ~~~" << std::endl;

}
//____________________________________________________________________

//___________________________________________________________________________
// coppied verbatum //
void MECTensorGenerator::AddTargetRemnant(GHepRecord * event) const
{
// Add the remnant nucleus (= initial nucleus - nucleon cluster) in the
// event record.


  std::cout << "~~~ Adding Remnant ~~~" << std::endl;
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

  std::cout << "~~~ REMNANT DONE ~~~" << std::endl;

}
//_________________________________________________________________________
void MECTensorGenerator::GenerateInitialHadrons  (GHepRecord * event) const
{
  // Earlier version of the code separated the GenerateInitialHadrons from
  // generating the recoil hadrons.
  // But we need a kinematic limits accept/reject loop here, so its one method.

  std::cout << "~~~ PLAYING WITH HADRONS ~~~" << std::endl;
  // -- Inputs: Q4, Fermi Momentum ------------------------N--
  // get neutrino & its 4-momentum
  GHepParticle * neutrino = event->Probe();
  assert(neutrino);
  TLorentzVector p4nu(*neutrino->P4());

  // get final state primary lepton & its 4-momentum
  GHepParticle * fsl = event->FinalStatePrimaryLepton();
  assert(fsl);
  TLorentzVector p4l(*fsl->P4());

  // calculate 4-momentum transfer from these
  TLorentzVector Q4 = p4nu - p4l;

  // get the target two-nucleon cluster and nucleus.
  // the remnant nucleus is apparently set, except for its momentum.
  GHepParticle * target_nucleus = event->TargetNucleus();
  assert(target_nucleus);
  GHepParticle * initial_nucleon_cluster = event->HitNucleon();
  assert(initial_nucleon_cluster);
  GHepParticle * remnant_nucleus = event->RemnantNucleus();
  assert(remnant_nucleus);

  // -- get nucleons from fermi gas
  // instantiate an empty local target nucleus, so I can use existing methods
  // to get a momentum from the prevailing Fermi-motion distribution 
  Target tgt(target_nucleus->Pdg());
  //RIK NucleonClusterConstituents is an implementation within this class.
  //RIK the nucleon_cluster->Pdg() was set earlier based on the pn probability.
  //RIK  as yet, this method doesn't know who Fermi is
  PDGCodeList pdgv = this->NucleonClusterConstituents(initial_nucleon_cluster->Pdg());
  assert(pdgv.size()==2);


  // These things need to be saved through to the end of the accept loop.
  bool accept = false;
  TVector3 p31i;
  TVector3 p32i;
  unsigned int iter = 0;

  int initial_nucleon_cluster_pdg = initial_nucleon_cluster->Pdg();
  int final_nucleon_cluster_pdg = 0;
  // Get the outgoing cluster mass.
  // the method in the next line does not work for some reason.
  //final_nucleon_cluster_pdg = event->Summary()->RecoilNucleonPdg();
  //std::cout << "RIK final_nucleon " << final_nucleon_cluster_pdg << " and initial " << initial_nucleon_cluster->Pdg() << std::endl;
  //std::cout << "RIK final_nucleon " << final_nucleon_cluster_pdg << " and initial " << initial_nucleon_cluster_pdg << std::endl;
  //heisenbug here.  May be related to shoehorning this into just the 200 or 202 initial states, not 201.

  // This is lingering here because I once had a problem with RecoilNucleonPdg() method.
  // override that method for now, but really go and fix it.
  if(neutrino->Pdg() > 0){
    if(initial_nucleon_cluster->Pdg() == kPdgClusterNP)final_nucleon_cluster_pdg = kPdgClusterPP;
    else if(initial_nucleon_cluster->Pdg() == kPdgClusterNN)final_nucleon_cluster_pdg = kPdgClusterNP;
    else LOG("MEC", pERROR) << "Wrong pdg for a CC neutrino MEC interaction" << initial_nucleon_cluster->Pdg();
  } else if(neutrino->Pdg() < 0){
    if(initial_nucleon_cluster->Pdg() == kPdgClusterNP)final_nucleon_cluster_pdg = kPdgClusterNN;
    else if(initial_nucleon_cluster->Pdg() == kPdgClusterPP)final_nucleon_cluster_pdg = kPdgClusterNP;
    else LOG("MEC", pERROR) << "Wrong pdg for a CC anti-neutrino MEC interaction" << initial_nucleon_cluster->Pdg();
  }
  

  if(final_nucleon_cluster_pdg != kPdgClusterPP && final_nucleon_cluster_pdg != kPdgClusterNN
     && final_nucleon_cluster_pdg != kPdgClusterNP && final_nucleon_cluster_pdg != initial_nucleon_cluster_pdg){
    LOG("MEC", pERROR) << "Wrong pdg for a CC neutrino MEC interaction initial "
		       << initial_nucleon_cluster->Pdg() << " final " << final_nucleon_cluster_pdg;
    // Replace this crude fail with the right one for GENIE.
    assert(false);
  }

  TLorentzVector p4initial_cluster;
  TLorentzVector p4final_cluster;
  TLorentzVector p4remnant_nucleus;
    
  //===========================================================================
  // Choose two nucleons from the prevailing fermi-motion distribution.
  // Some produce kinematically unallowed combinations initial cluster and Q2
  // Find out, and if so choose them again with this accept/reject loop.
  // Some kinematics are especially tough 
  while(!accept){
    iter++;
    if(iter > kRjMaxIterations) {
      // error if try too many times
      LOG("MEC", pWARN)
           << "Couldn't select a valid W, Q^2 pair after " 
           << iter << " iterations";
        event->EventFlags()->SetBitNumber(kKineGenErr, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select initial hadron kinematics");
        exception.SwitchOnFastForward();
        throw exception;
    }

    // generate nucleons
    // tgt is a Target object for local use, just waiting to be filled.
    // this sets the pdg of each nucleon and its momentum from a Fermi-gas
    // Nieves et al. would use a local Fermi gas here, not this, but ok.
    tgt.SetHitNucPdg(pdgv[0]);
    fNuclModel->GenerateNucleon(tgt);
    p31i = fNuclModel->Momentum3();
    tgt.SetHitNucPdg(pdgv[1]);
    fNuclModel->GenerateNucleon(tgt);
    p32i = fNuclModel->Momentum3();

    // Now write down the initial cluster four-vector for this choice
    TVector3 p3i = p31i + p32i;
    double mass2 = PDGLibrary::Instance()->Find( initial_nucleon_cluster_pdg )->Mass();
    mass2 *= mass2;
    double energy = TMath::Sqrt(p3i.Mag2() + mass2);
    p4initial_cluster.SetPxPyPzE(p3i.Px(),p3i.Py(),p3i.Pz(),energy);
    std::cout << "RIK initial "
	      << p4initial_cluster.Px() << " " << p4initial_cluster.Py() << " " << p4initial_cluster.Pz() << " "
	      << p4initial_cluster.E() << " " << p4initial_cluster.M() << " " << initial_nucleon_cluster_pdg << std::endl;
    

    // calculate binding energy // for now just value from neut
    // cast it as a TLorentzVector with only an energy component, for the upcoming addition.
    // ebind is intrinsically negative here.
    // Is that the convention of the GENIE implementation of the local Fermi gas or spectral function?
    // This could be different for protons and neutrons, or for local Fermi gas, or per nucleus.
    double ebind1 = -0.0250;
    double ebind2 = -0.0250;
    TLorentzVector tLVebind(0., 0., 0., ebind1 + ebind2);
    
    //RIK why is this the right place to subtract ebin ?
    // its okay.  physically, I'm subtracting it from q.
    // the energy transfer to the nucleon is 50 MeV less.
    // the energy cost to put this cluster on-shell.
    
    // calculate recoil nucleon cluster 4-momentum
    // tLVebind is intrinsically negative, as it was before.
    p4final_cluster = p4initial_cluster + Q4 + tLVebind;
    //TLorentzVector p4cluster_recoil = p4cluster + Q4 + ebind*2.0;  // + tLVebind;
    
    // The test here is that the resulting four-vector corresponds to a high-enough invariant mass.

    //std::cout << "RIK Q4 " << Q4.Px() << " " << Q4.Py() << " " << Q4.Pz() << " " << Q4.E() << std::endl;

    //std::cout << "RIK accept/reject " << p4final_cluster.Px() << " "
    //	      << p4final_cluster.Py() << " "
    //  	      << p4final_cluster.Pz() << " "
    //  	      << p4final_cluster.E() << " "
    //  	      << p4final_cluster.M() << " "
    //      << PDGLibrary::Instance()->Find( final_nucleon_cluster_pdg )->Mass() << " " << final_nucleon_cluster_pdg << std::endl;
    
    if(p4final_cluster.M() < PDGLibrary::Instance()->Find( final_nucleon_cluster_pdg )->Mass()){
      accept = false;      
    } else {
      accept = true;
    }

  }

  // Now everything above should be written to ghep.
  // First the initial state nucleons.

  
  // Initial nucleon cluster
  // These two are redundant maybe? 
  initial_nucleon_cluster->SetMomentum(p4initial_cluster);

  //event->Summary()->InitStatePtr()->TgtPtr()->SetHitNucP4(p4initial_cluster);

  double Mi  = PDGLibrary::Instance()->Find(target_nucleus->Pdg() )-> Mass(); // initial nucleus mass
  remnant_nucleus->SetMomentum(-1.0*p4initial_cluster.Px(),
			       -1.0*p4initial_cluster.Py(),
			       -1.0*p4initial_cluster.Pz(),
			       Mi - p4initial_cluster.E());


  // Now the final nucleon cluster.
  
  //RIK I remember this can be used in a local Fermi gas implementation too.
  TLorentzVector v4(*neutrino->X4());
  
  // add to the event record
  event->AddParticle( final_nucleon_cluster_pdg, kIStDecayedState, 
		      2, -1, -1, -1, p4final_cluster, v4); 

  std::cout << "~~~ HADRONS DONE ~~~" << std::endl;

}

/*
//_________________________________________________________________________
void MECTensorGenerator::RecoilNucleonCluster    (GHepRecord * event) const
{
    // RIK removed this old method.
}
*/


//_________________________________________________________________________
void MECTensorGenerator::DecayNucleonCluster  (GHepRecord * event) const
{
  std::cout << "~~~ DECAY CLUSTER ~~~" << std::endl;

// Perform a phase-space decay of the nucleon cluster and add its decay
// products in the event record

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
    std::cout << "RIK fPhaseSpaceGenerator says this is not permitted." << std::endl;
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
       std::cout << "RIK fPhaseSpaceGenerator couldn't do this." << std::endl;
       // clean up
       delete [] mass;
       delete p4d;
       delete v4d;
       // throw exception
       LOG("MEC", pWARN)
	 << "Couldn't decay after " 
	 << itry << " iterations";
       genie::exceptions::EVGThreadException exception;
       exception.SetReason("Couldn't select decay after N attempts");
       exception.SwitchOnFastForward();
       throw exception;
     }
     double w  = fPhaseSpaceGenerator.Generate();   
     if(w > wmax) {
       //RIK Jackie deleted a warning here...
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
     event->AddParticle(pdgc, ist, nucleon_cluster_id,-1,-1,-1, *p4fin, v4);
     idp++;
  }

  // Clean-up
  delete [] mass;
  delete p4d;
  delete v4d;

  std::cout << "~~~ CLUSTER DECAY DONE ~~~" << std::endl;
}
//_________________________________________________________________________
PDGCodeList MECTensorGenerator::NucleonClusterConstituents(int pdgc) const
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
void MECTensorGenerator::Configure(const Registry & config)   
{
  Algorithm::Configure(config);
  this->LoadConfig();
} 
//___________________________________________________________________________ 
void MECTensorGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//___________________________________________________________________________
void MECTensorGenerator::LoadConfig(void)
{
  fNuclModel = 0;
      
  RgKey nuclkey = "NuclearModel";
  fNuclModel = dynamic_cast<const NuclearModelI *> (this->SubAlg(nuclkey));
  assert(fNuclModel);
}
//___________________________________________________________________________
















