// -*- C++ -*-
//
// Package:    AODSIM/AODSIMAnalyzer
// Class:      AODSIMAnalyzer
// 
/**\class AODSIMAnalyzer AODSIMAnalyzer.cc AODSIM/AODSIMAnalyzer/plugins/AODSIMAnalyzer.cc

 Description: [This is written for converting the edm format file to the root format]

 Implementation:
     May 24, 2015	:	Need to put the TriggerResults
*/
//
// Original Author:  Ram Krishna Sharma
//         Created:  Sat, 23 May 2015 19:14:00 GMT
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
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"

#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// Header file for TriggerResults
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

// Header file for Pileup info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

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

#include "FullSimAnalysis/AODSIMAnalyzer/plugins/AODSIMAnalyzer.h"
//
// class declaration
//

// Moved to AODSIMAnalyzer.h
//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
AODSIMAnalyzer::AODSIMAnalyzer(const edm::ParameterSet& iConfig) :
//
//	List of Input Tags
//
//	Note :	For every InputTag defined or a variables we need to also define it in the 
//		class AODSIMAnalyzer
//
//
GenPart_( iConfig.getParameter<edm::InputTag>( "GenPart" ) ),
PileUp_( iConfig.getParameter<edm::InputTag>( "PileUp" ) ),
GenJetAk4_(iConfig.getParameter<edm::InputTag>( "GenJetAk4")),
GenJetsAk4NoNu_(iConfig.getParameter<edm::InputTag>( "GenJetsAk4NoNu")),
GenMetTru_(iConfig.getParameter<edm::InputTag>( "GenMetTru")),
GenMetCal_(iConfig.getParameter<edm::InputTag>( "GenMetCal")),
trigTag_(iConfig.getParameter<edm::InputTag>("trigTag")),
triggerEventTag_(iConfig.getParameter<edm::InputTag>("triggerEventTag")),
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


AODSIMAnalyzer::~AODSIMAnalyzer()
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
AODSIMAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;

   //vectors to pat objects
   runNumber		= iEvent.id().run();
   eventNumber		= iEvent.id().event();
   orbitNumber		= iEvent.orbitNumber();
   bunchCrossing	= iEvent.bunchCrossing();
   lumiBlockNumber	= iEvent.luminosityBlock();
   isData		= iEvent.isRealData();

   if (iEvent.isRealData() ){
     runPileupinfo_ = false;
   }
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


  if(runPileupinfo_)
  {
  // Few References: https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_4_1/doc/html/d9/d53/classPileupSummaryInfo.html
  //
  Handle<std::vector<PileupSummaryInfo>>  PupInfo;
  iEvent.getByLabel(PileUp_, PupInfo);

  std::vector<PileupSummaryInfo>::const_iterator PVI;

  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
  if (Verbose_)
  std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
  //pileup_bunchXing = PVI->getBunchCrossing();
  //numberOfPUVertices      = PVI->getPU_NumInteractions();
  pileup_bunchXing_.push_back(PVI->getBunchCrossing());
  numberOfPUVertices_.push_back(PVI->getPU_NumInteractions());
  numberOfPUVerticesMixingTruth_.push_back(PVI->getTrueNumInteractions());
  }
  }	//	if(runPileupinfo_) 

  // Get the TriggerEvent
//  Handle<trigger::TriggerEvent> triggerEventHandle;
//  iEvent.getByLabel(triggerEventTag_,triggerEventHandle);
//
//  if (!triggerEventHandle.isValid()) {
//  	std::cout << "Error in getting TriggerEvent product from Event!" << std::endl;
//	return;
//  }
   // A HLT path consists of many different modules - producers and filters     
   // The event can be rejected at any filter stage along the path     
   //                                                                                                                                                 
   // As well as just the basic pass/fail info,       
   // the TriggerResults object stores the index of the module in the path        
   // which made the final decision
   // In the case where the event was accepted, this index is therefore just    
   // the index of the last module along the path  









   // get the GenParticleCollection handle
   Handle<reco::GenParticleCollection> genParticles;
   iEvent.getByLabel(GenPart_, genParticles);

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
   	genMuonPt_.push_back(p.pt());
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
   iEvent.getByLabel(GenJetAk4_ ,genjets);

   Handle<reco::PFJetCollection> pfjets;
   iEvent.getByLabel( "ak4PFJets", pfjets);
// PFJetCollection pfjc=*(pfjets.product());
   
   if (genjets.isValid()) {
   	const reco::GenJetCollection* mygenjets = &(*genjets);
   	std::vector<const reco::GenJet *> sortedPtrs;
   	sortedPtrs.reserve(mygenjets->size());
  	
	nJets = mygenjets->size();

   	for (const reco::GenJet  &g : *mygenjets){
   	sortedPtrs.push_back(&g);
   	}
   
   	std::sort(sortedPtrs.begin(), sortedPtrs.end(),PtGreater());
   
   
   	if(mygenjets->size() >= 4)
   	for (auto const & genPtr : sortedPtrs) {
   	auto const & geni = *genPtr;
   
   	GSJetPt_	.push_back(geni.pt()	);	
   	GSJetEta_	.push_back(geni.eta()	);	
   	GSJetPhi_	.push_back(geni.phi()	);	
   	GSJetVx_	.push_back(geni.vx()	);	
   	GSJetVy_	.push_back(geni.vy()	);	
   	GSJetVz_	.push_back(geni.vz()	);	
   	GSJetCharge_	.push_back(geni.charge());	
   	GSJetMass_	.push_back(geni.mass()	);	

	//Store some Calo Jet Information:
//	if( const CaloJet* caloJet = dynamic_cast<const CaloJet*>(&(*(geni.get()) )) )
//	GSJetEMFraction_.push_back(caloJet->emEnergyFraction());
//	else GSJetEMFraction_.push_back(-1.0);



	if (Verbose_)
   	std::cout<<"jet px = "<<geni.pt()<<std::endl;
	//std::cout<<"Ecal energy = "<<geni.emEnergy()<<std::endl;
	//std::cout<<"HCal energy = "<<geni.hadEnergy()<<std::endl;
	//std::cout<<"Invisible energy = "<<geni.invisibleEnergy()<<std::endl;
	//std::cout<<"auxiliary energy = "<<geni.auxiliaryEnergy()<<std::endl;
	std::cout<<"jet Area = "<<geni.jetArea()<<std::endl;
   	}
   
   }//if (genjets.isValid()) 


   // get the handle
   Handle<reco::GenJetCollection> genjetsNoNu;
   iEvent.getByLabel(GenJetsAk4NoNu_ ,genjetsNoNu);
   
   if (genjetsNoNu.isValid()) {
   	const reco::GenJetCollection* mygenjetsNoNu = &(*genjetsNoNu);
   	std::vector<const reco::GenJet *> sortedPtrs;
   	sortedPtrs.reserve(mygenjetsNoNu->size());

	nJetsNoNu = mygenjetsNoNu->size();
   
   	for (const reco::GenJet  &g : *mygenjetsNoNu){
   	sortedPtrs.push_back(&g);
   	}
   
   	std::sort(sortedPtrs.begin(), sortedPtrs.end(),PtGreater());
   
   
   	if(mygenjetsNoNu->size() >= 4)
   	for (auto const & genPtr : sortedPtrs) {
   	auto const & geni = *genPtr;
   
   	GSJetNoNuPt_.push_back( geni.pt());	
   	GSJetNoNuEta_.push_back(geni.eta());	
   	GSJetNoNuPhi_.push_back(geni.phi());	
   	GSJetNoNuVx_.push_back(geni.vx());	
   	GSJetNoNuVy_.push_back(geni.vy());	
   	GSJetNoNuVz_.push_back(geni.vz());	
   	GSJetNoNuCharge_.push_back(geni.charge());	
   	GSJetNoNuMass_.push_back(geni.mass());	
	if (Verbose_)
   	std::cout<<"jet px = "<<geni.pt()<<std::endl;
   	}
   
   }//if (genjetsNoNu.isValid()) 
   
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
void AODSIMAnalyzer::AddBranch(std::vector<double>* vec, std::string name){
	tree->Branch(name.c_str(),vec);
}
void AODSIMAnalyzer::AddBranch(std::vector<int>* vec, std::string name){
	tree->Branch(name.c_str(),vec);
}
void AODSIMAnalyzer::AddBranch(int* vec, std::string name){
	tree->Branch(name.c_str(),vec);
}

void AODSIMAnalyzer::SetBranches(){
	//AddBranch(&,	"");
	AddBranch(&runNumber,		"runNumber");
	AddBranch(&eventNumber,	"eventNumber");
	AddBranch(&orbitNumber,	"orbitNumber");
	AddBranch(&bunchCrossing,		"bunchCrossing");
	AddBranch(&lumiBlockNumber,	"lumiBlockNumber");
	AddBranch(&isData,	"isData");

	// Pile Up info
	AddBranch(&numberOfPUVertices_,"numberOfPUVertices");
	AddBranch(&numberOfPUVerticesMixingTruth_,"numberOfPUVerticesMixingTruth");
	AddBranch(&pileup_bunchXing_,"pileup_bunchXing");
	//
	AddBranch(&GSLepPt_,	"GSLepPt");
	AddBranch(&GSLepEta_,	"GSLepEta");
	AddBranch(&GSLepPhi_,	"GSLepPhi");
	
	AddBranch(&nJets,	"nJets");
	AddBranch(&GSJetPt_,	"GSJetPt");
	AddBranch(&GSJetEta_,	"GSJetEta");
	AddBranch(&GSJetPhi_,	"GSJetPhi");

	AddBranch(&nJetsNoNu,	"nJetsNoNu");
	AddBranch(&GSJetNoNuPt_,	"GSJetNoNuPt");
	AddBranch(&GSJetNoNuEta_,	"GSJetNoNuEta");
	AddBranch(&GSJetNoNuPhi_,	"GSJetNoNuPhi");
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

void AODSIMAnalyzer::Clear(){		// Clear only those variable which is decleared as vector
	pileup_bunchXing_.clear();
	numberOfPUVertices_.clear();
	numberOfPUVerticesMixingTruth_.clear();
	GSLepPt_.clear();
	GSLepEta_.clear();
	GSLepPhi_.clear();
	GSJetPt_.clear();
	GSJetEta_.clear();
	GSJetPhi_.clear();
	GSJetNoNuPt_.clear();
	GSJetNoNuEta_.clear();
	GSJetNoNuPhi_.clear();
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
AODSIMAnalyzer::beginJob()
{
    std::cout<<"Inside beginJob()"<<std::endl;

    if(wantLocalFile_)
    {
	std::cout<<"inside rootFilename_"<<std::endl;

	outputFile_ = new TFile(loutputFile_.c_str(),"RECREATE");    
        outputFile_->SetCompressionLevel(2);
	tree = new TTree("RawSim","AODSIM Info");
    }
    if(wantRFIOFile_)
    {
	std::cout<<"inside rootFilename_"<<std::endl;

	outputFile_ = new TFile(rfoutputFile_.c_str(),"RECREATE");    
        outputFile_->SetCompressionLevel(2);
	tree = new TTree("RawSim","AODSIM Info");
    }


	SetBranches();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
AODSIMAnalyzer::endJob() 
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
AODSIMAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
AODSIMAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
AODSIMAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
AODSIMAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
AODSIMAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(AODSIMAnalyzer);
