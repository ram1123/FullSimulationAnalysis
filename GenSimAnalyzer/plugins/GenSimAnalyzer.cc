// -*- C++ -*-
//
// Package:    GEN-SIM/GenSimAnalyzer
// Class:      GenSimAnalyzer
// 
/**\class GenSimAnalyzer GenSimAnalyzer.cc GEN-SIM/GenSimAnalyzer/plugins/GenSimAnalyzer.cc

 Description: [	This file reads the edm file, output of GEN-SIM stage, and gives us root file 
 		with certain basic information about the particles]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  ramkrishna sharma
//         Created:  Thu, 26 Feb 2015 14:31:50 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// Header file for TriggerResults
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"


#include <vector>

#include "TTree.h"
#include "TFile.h"
#include <cmath>
#include "TLorentzVector.h"
#include <fstream>
#include <iostream>
#include "TVector.h"
 
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//
// class declaration
//

class GenSimAnalyzer : public edm::EDAnalyzer {
   public:
      explicit GenSimAnalyzer(const edm::ParameterSet&);
      ~GenSimAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      void SetBranches();
      void Clear();


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      virtual void AddBranch(std::vector<double>*, std::string name);
      virtual void AddBranch(std::vector<int>*, std::string name);
      virtual void AddBranch(int* vec, std::string name);

      // ----------member data ---------------------------
      
      edm::InputTag GenPart_;
      edm::InputTag GenJetAk4_;
      edm::InputTag GenMetTru_;
      edm::InputTag GenMetCal_;
      edm::InputTag trigTag_;

      // Few Boolean InputTags
      bool      Verbose_;
      bool      SortGen_;
      bool      wantLocalFile_; 
      bool      wantRFIOFile_;  
      std::string       loutputFile_;
      std::string       rfoutputFile_;

      edm::Service<TFileService> fs;
      TFile * outputFile_;
      TTree* tree;

      // Event Properties                                                                                                                             
      int run;
      int event;
      int orbit;
      int bx;
      int lumis;
      int isData;
 

      std::vector<double> GSLepVx_	;
      std::vector<double> GSLepVy_	;
      std::vector<double> GSLepVz_	;
      std::vector<double> GSLepCharge_	;
      std::vector<double> GSLepMass_	;
      std::vector<double> GSLepPt_	;
      std::vector<double> GSLepEta_	;
      std::vector<double> GSLepPhi_	;
      std::vector<double> GSLepTheta_	;

      std::vector<double> GSJetVx_	;
      std::vector<double> GSJetVy_	;
      std::vector<double> GSJetVz_	;
      std::vector<double> GSJetCharge_	;
      std::vector<double> GSJetMass_	;
      std::vector<double> GSJetPt_	;
      std::vector<double> GSJetEta_	;
      std::vector<double> GSJetPhi_	;
      std::vector<double> GSJetTheta_	;
       
      std::vector<double> GSMetTruVx_   ;
      std::vector<double> GSMetTruVy_   ;
      std::vector<double> GSMetTruVz_   ;
      std::vector<double> GSMetTruCharge_       ;
      std::vector<double> GSMetTruMass_ ;
      std::vector<double> GSMetTruPt_   ;
      std::vector<double> GSMetTruEta_  ;
      std::vector<double> GSMetTruPhi_  ;
      std::vector<double> GSMetTruTheta_        ;
       
      std::vector<double> GSMetCalVx_   ;
      std::vector<double> GSMetCalVy_   ;
      std::vector<double> GSMetCalVz_   ;
      std::vector<double> GSMetCalCharge_       ;
      std::vector<double> GSMetCalMass_ ;
      std::vector<double> GSMetCalPt_   ;
      std::vector<double> GSMetCalEta_  ;
      std::vector<double> GSMetCalPhi_  ;
      std::vector<double> GSMetCalTheta_        ;
       
      std::vector<double> GSMetCNPVx_   ;
      std::vector<double> GSMetCNPVy_   ;
      std::vector<double> GSMetCNPVz_   ;
      std::vector<double> GSMetCNPCharge_       ;
      std::vector<double> GSMetCNPMass_ ;
      std::vector<double> GSMetCNPPt_   ;
      std::vector<double> GSMetCNPEta_  ;
      std::vector<double> GSMetCNPPhi_  ;
      std::vector<double> GSMetCNPTheta_        ;

      std::vector<double> GSMetTrueVx_  ;
      std::vector<double> GSMetTrueVy_  ;
      std::vector<double> GSMetTrueVz_  ;
      std::vector<double> GSMetTrueCharge_      ;
      std::vector<double> GSMetTrueMass_        ;
      std::vector<double> GSMetTruePt_  ;
      std::vector<double> GSMetTrueEta_ ;
      std::vector<double> GSMetTruePhi_ ;
      std::vector<double> GSMetTrueTheta_       ;
     };

class PtGreater {
         public:
         template <typename T> bool operator () (const T& i, const T& j) {
         return (i->pt() > j->pt());
         }
};


//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenSimAnalyzer::GenSimAnalyzer(const edm::ParameterSet& iConfig):
//
//	List of Input Tags
//
//	Note :	For every InputTag defined or a variables we need to also define it in the 
//		class GenSimAnalyzer
//
//
GenPart_( iConfig.getParameter<edm::InputTag>( "GenPart" ) ),
GenJetAk4_(iConfig.getParameter<edm::InputTag>( "GenJetAk4")),
GenMetTru_(iConfig.getParameter<edm::InputTag>( "GenMetTru")),
GenMetCal_(iConfig.getParameter<edm::InputTag>( "GenMetCal")),
trigTag_(iConfig.getParameter<edm::InputTag>("trigTag")),
//genMetCNP_(iConfig.getParameter<edm::InputTag>( "genMetCNP")),
//
//	Few General Boolean or Variables
//
Verbose_(iConfig.getUntrackedParameter<bool>("Verbose",0)),
SortGen_(iConfig.getUntrackedParameter<bool>("SortGen",0)),
wantLocalFile_(iConfig.getUntrackedParameter<int>("wantLocalFile",1)),
wantRFIOFile_(iConfig.getUntrackedParameter<int>("wantRFIOFile",0)),
loutputFile_(iConfig.getUntrackedParameter<std::string>("loutputFile", "gsftrack.root")),
rfoutputFile_(iConfig.getUntrackedParameter<std::string>("rfoutputFile", "/uscms_data/d2/sushil/CMSSW/MonoPhoton/CMSSW_3_8_6/src/QCDFakeRate/Analyser/test/gsftrack.root"))
{
   //now do what ever initialization is needed
}


GenSimAnalyzer::~GenSimAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
   if(wantLocalFile_) 
   {
   delete outputFile_;
   }
   if(wantRFIOFile_) 
   {
   delete outputFile_;
   }
}


//
// member functions
//

// ------------ method called for each event  ------------
void
GenSimAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   //vectors to pat objects
   run    = iEvent.id().run();
   event  = iEvent.id().event();
   orbit  = iEvent.orbitNumber();
   bx     = iEvent.bunchCrossing();
   lumis  = iEvent.luminosityBlock();
   isData = iEvent.isRealData();

   // get the TriggerResults handle
   edm::Handle<edm::TriggerResults> trigResults;
   if (not iEvent.getByLabel(trigTag_, trigResults)) {
   std::cout << ">>> TRIGGER collection does not exist !!!\n";
   return;
   }
   const edm::TriggerNames & trigNames = iEvent.triggerNames(*trigResults);

   for (unsigned int i=0; i<trigResults->size(); i++)
   {
   std::string trigName = trigNames.triggerName(i);
  int trigResult = trigResults->accept(i); //bool not to use
  //if (Verbose_)
  if (trigResult)
  if (Verbose_)
  std::cout <<"Name of Trigger = " << trigName <<"  Trigger Result = " << trigResult <<"  Trigger Number = " << i << std::endl;
  }	// for (unsigned int i=0; i<trigResults->size(); i++)


// get the handle
Handle<reco::GenParticleCollection> genParticles;
iEvent.getByLabel("genParticles", genParticles);



#if SortGen_
	const reco::GenParticleCollection* genColl= &(*genParticles);
	std::vector<const reco::GenParticle *> sortedPtrs;
	sortedPtrs.reserve(genColl->size());
	for (const reco::GenParticle &g : *genColl) { sortedPtrs.push_back(&g); }
	     std::sort(sortedPtrs.begin(), sortedPtrs.end(),PtGreater());

	for (auto const & genPtr : sortedPtrs) {
	   auto const & geni = *genPtr;

	if (Verbose_)
	std::cout<<"pt = "<<geni.pt()<<std::endl;

	}
#else

   // loop over particles
   for(size_t i = 0; i < genParticles->size(); ++ i) {
   
   	// the reference p to the i-th particle:
   	const GenParticle & p = (*genParticles)[i];
   
   	// get pdgId:
   	int id = p.pdgId();
   
   #if 0
   	
   	// get status:
   	int st = p.status();
   	// get pointer to mother (reco::Candidate type!):
   	const Candidate * mom = p.mother();
   	// get pt, eta, phi, mass:
   	double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
   	GSMuonPt_.push_back(p.pt());
   	// get vertex components:
   	double vx = p.vx(), vy = p.vy(), vz = p.vz();
   	// get charge:
   	int charge = p.charge();
   	// get no. of daughters:
   	int n = p.numberOfDaughters();	
	if (Verbose_)
   	std::cout<<"id = "<<id<<"\tst = "<<st<<"\t pt , eta, phi, mass = "<<pt<<eta<<phi<<mass<<"\tvx,vy,vz = "<<vx<<vy<<vz<<"\tcharge = "<<charge<<"\t n = "<<n<<"\tMom = "<<&mom<<std::endl;
   #else
   	if ( abs(id)==11 || abs(id)==13)
   	{
   	GSLepPt_.push_back(p.pt());	
   	GSLepEta_.push_back(p.eta());	
   	GSLepPhi_.push_back(p.phi());	
   	GSLepVx_.push_back(p.vx());	
   	GSLepVy_.push_back(p.vy());	
   	GSLepVz_.push_back(p.vz());	
   	GSLepCharge_.push_back(p.charge());	
   	GSLepMass_.push_back(p.mass());	
	if (Verbose_)
   	std::cout<<"In GenParticle Loop"<<std::endl;	 
   	}
   #endif
   
   }	//for(size_t i = 0; i < genParticles->size(); ++ i) 
#endif

// get the handle
Handle<reco::GenJetCollection> genjets;
iEvent.getByLabel("ak5GenJets",genjets);

if (genjets.isValid()) {
	const reco::GenJetCollection* mygenjets = &(*genjets);
	std::vector<const reco::GenJet *> sortedPtrs;
	sortedPtrs.reserve(mygenjets->size());

	for (const reco::GenJet  &g : *mygenjets){
	sortedPtrs.push_back(&g);
	}

	std::sort(sortedPtrs.begin(), sortedPtrs.end(),PtGreater());


	if(mygenjets->size() >= 4)
	for (auto const & genPtr : sortedPtrs) {
	auto const & geni = *genPtr;

	GSJetPt_.push_back( geni.pt());	
	GSJetEta_.push_back(geni.eta());	
	GSJetPhi_.push_back(geni.phi());	
	GSJetVx_.push_back(geni.vx());	
	GSJetVy_.push_back(geni.vy());	
	GSJetVz_.push_back(geni.vz());	
	GSJetCharge_.push_back(geni.charge());	
	GSJetMass_.push_back(geni.mass());	
	//std::cout<<"jet px = "<<geni.pt()<<std::endl;
	}

}//if (genjets.isValid()) 


   
   Handle<reco::GenMETCollection> genMetTru;
   iEvent.getByLabel(GenMetTru_ ,genMetTru);
   if (genMetTru.isValid()) {
   	typedef reco::GenMETCollection::const_iterator gmiter;
   	for ( gmiter i=genMetTru->begin(); i!=genMetTru->end(); i++) {
   
   	GSMetTruPt_.push_back( i->pt());	
   	GSMetTruEta_.push_back(i->eta());	
   	GSMetTruPhi_.push_back(i->phi());	
   	GSMetTruVx_.push_back(i->vx());	
   	GSMetTruVy_.push_back(i->vy());	
   	GSMetTruVz_.push_back(i->vz());	
   	GSMetTruCharge_.push_back(i->charge());	
   	GSMetTruMass_.push_back(i->mass());	
   	
   	}
   }	//if (genMetTru.isValid()) {
   
   Handle<reco::GenMETCollection> genMetCal;
   iEvent.getByLabel(GenMetCal_ ,genMetCal);
   if (genMetCal.isValid()) {
   	typedef reco::GenMETCollection::const_iterator gmiter;
   	for ( gmiter i=genMetCal->begin(); i!=genMetCal->end(); i++) {
   
   	GSMetCalPt_.push_back( i->pt());	
   	GSMetCalEta_.push_back(i->eta());	
   	GSMetCalPhi_.push_back(i->phi());	
   	GSMetCalVx_.push_back(i->vx());	
   	GSMetCalVy_.push_back(i->vy());	
   	GSMetCalVz_.push_back(i->vz());	
   	GSMetCalCharge_.push_back(i->charge());	
   	GSMetCalMass_.push_back(i->mass());	
   	
   	}
   }	//if (genMetCal.isValid()) {
   
//   Handle<reco::GenMETCollection> genMetCNP;
//   iEvent.getByLabel(genMetCNP_ ,genMetCNP);
//   if (genMetCNP.isValid()) {
//   	typedef reco::GenMETCollection::const_iterator gmiter;
//   	for ( gmiter i=genMetCNP->begin(); i!=genMetCNP->end(); i++) {
//   
//   	GSMetCNPPt_.push_back( i->pt());	
//   	GSMetCNPEta_.push_back(i->eta());	
//   	GSMetCNPPhi_.push_back(i->phi());	
//   	GSMetCNPVx_.push_back(i->vx());	
//   	GSMetCNPVy_.push_back(i->vy());	
//   	GSMetCNPVz_.push_back(i->vz());	
//   	GSMetCNPCharge_.push_back(i->charge());	
//   	GSMetCNPMass_.push_back(i->mass());	
//   	
//   	}
//   }	//if (genMetCNP.isValid()) {

tree->Fill();

}


// ---------Add branches ---------------------
void GenSimAnalyzer::AddBranch(std::vector<double>* vec, std::string name){
	tree->Branch(name.c_str(),vec);
}
void GenSimAnalyzer::AddBranch(std::vector<int>* vec, std::string name){
	tree->Branch(name.c_str(),vec);
}
void GenSimAnalyzer::AddBranch(int* vec, std::string name){
	tree->Branch(name.c_str(),vec);
}

void GenSimAnalyzer::SetBranches(){
	//AddBranch(&,	"");
	AddBranch(&run,		"run");
	AddBranch(&event,	"event");
	AddBranch(&orbit,	"orbit");
	AddBranch(&bx,		"bx");
	AddBranch(&lumis,	"lumis");
	AddBranch(&isData,	"isData");

	AddBranch(&GSLepPt_,	"GSLepPt");
	AddBranch(&GSLepEta_,	"GSLepEta");
	AddBranch(&GSLepPhi_,	"GSLepPhi");
	AddBranch(&GSJetPt_,	"GSJetPt");
	AddBranch(&GSJetEta_,	"GSJetEta");
	AddBranch(&GSJetPhi_,	"GSJetPhi");

	AddBranch(&GSMetTruPt_,	"GSMetPt");
	AddBranch(&GSMetTruEta_,	"GSMetEta");
	AddBranch(&GSMetTruPhi_,	"GSMetPhi");
	AddBranch(&GSMetCalPt_,	"GSMetPt");
	AddBranch(&GSMetCalEta_,	"GSMetEta");
	AddBranch(&GSMetCalPhi_,	"GSMetPhi");
	AddBranch(&GSMetCNPPt_,	"GSMetPt");
	AddBranch(&GSMetCNPEta_,	"GSMetEta");
	AddBranch(&GSMetCNPPhi_,	"GSMetPhi");
}

void GenSimAnalyzer::Clear(){
	GSLepPt_.clear();
	GSLepEta_.clear();
	GSLepPhi_.clear();
	GSJetPt_.clear();
	GSJetEta_.clear();
	GSJetPhi_.clear();
	GSMetTruPt_.clear();
	GSMetTruEta_.clear();
	GSMetTruPhi_.clear();
	GSMetCalPt_.clear();
	GSMetCalEta_.clear();
	GSMetCalPhi_.clear();
	GSMetCNPPt_.clear();
	GSMetCNPEta_.clear();
	GSMetCNPPhi_.clear();
}
// ------------ method called once each job just before starting event loop  ------------
void 
GenSimAnalyzer::beginJob()
{
   std::cout<<"Inside beginJob()"<<std::endl;

    if(wantLocalFile_)
    {
	std::cout<<"inside rootFilename_"<<std::endl;

	outputFile_ = new TFile(loutputFile_.c_str(),"RECREATE");    
        outputFile_->SetCompressionLevel(2);
	tree = new TTree("RawSim","RAW-SIM Info");
    }
    if(wantRFIOFile_)
    {
	std::cout<<"inside rootFilename_"<<std::endl;

	outputFile_ = new TFile(rfoutputFile_.c_str(),"RECREATE");    
        outputFile_->SetCompressionLevel(2);
	tree = new TTree("RawSim","RAW-SIM Info");
    }

	SetBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenSimAnalyzer::endJob() 
{
    if(wantLocalFile_) 
    {
	outputFile_->Write();
	outputFile_->Close();
    }
    if(wantRFIOFile_) 
    {
	outputFile_->Write();
	outputFile_->Close();
    }
}

// ------------ method called when starting to processes a run  ------------
/*
void 
GenSimAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
GenSimAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
GenSimAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
GenSimAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenSimAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GenSimAnalyzer);
