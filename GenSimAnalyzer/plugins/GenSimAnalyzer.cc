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

#include <vector>

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <fstream>
#include <iostream>


#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"


//#include "FWCore/ParameterSet/interface/InputTag.h"


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

      // ----------member data ---------------------------
      //edm::InputTag genParticles;
      edm::Service<TFileService> fs;
      TFile * file;
      TTree* tree;

      std::vector<double> genLepVx_	;
      std::vector<double> genLepVy_	;
      std::vector<double> genLepVz_	;
      std::vector<double> genLepCharge_	;
      std::vector<double> genLepMass_	;
      std::vector<double> genLepPt_	;
      std::vector<double> genLepEta_	;
      std::vector<double> genLepPhi_	;
      std::vector<double> genLepTheta_	;

      std::vector<double> genJetVx_	;
      std::vector<double> genJetVy_	;
      std::vector<double> genJetVz_	;
      std::vector<double> genJetCharge_	;
      std::vector<double> genJetMass_	;
      std::vector<double> genJetPt_	;
      std::vector<double> genJetEta_	;
      std::vector<double> genJetPhi_	;
      std::vector<double> genJetTheta_	;

      std::vector<double> genMetVx_	;
      std::vector<double> genMetVy_	;
      std::vector<double> genMetVz_	;
      std::vector<double> genMetCharge_	;
      std::vector<double> genMetMass_	;
      std::vector<double> genMetPt_	;
      std::vector<double> genMetEta_	;
      std::vector<double> genMetPhi_	;
      std::vector<double> genMetTheta_	;



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
GenSimAnalyzer::GenSimAnalyzer(const edm::ParameterSet& iConfig)
//:genParticles(iConfig.getUntrackedParameter<edm::InputTag>("GeneratorParticles"))
{
   //now do what ever initialization is needed
   #if 0
   bool wantLocalFile_=1;
   #else
   #endif


}


GenSimAnalyzer::~GenSimAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

   #if 0
   if(wantLocalFile_)
   {
	lf->Close();
	delete lf;
   }
   #else
   #endif

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
// get the handle
Handle<reco::GenParticleCollection> genParticles;
iEvent.getByLabel("genParticles", genParticles);


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

//	std::cout<<"id = "<<id<<"\tst = "<<st<<"\t pt , eta, phi, mass = "<<pt<<eta<<phi<<mass<<"\tvx,vy,vz = "<<vx<<vy<<vz<<"\tcharge = "<<charge<<"\t n = "<<n<<"\tMom = "<<&mom<<std::endl;

#else
	if ( abs(id)==11 || abs(id)==13)
	{
	genLepPt_.push_back(p.pt());	
	genLepEta_.push_back(p.eta());	
	genLepPhi_.push_back(p.phi());	
	genLepVx_.push_back(p.vx());	
	genLepVy_.push_back(p.vy());	
	genLepVz_.push_back(p.vz());	
	genLepCharge_.push_back(p.charge());	
	genLepMass_.push_back(p.mass());	
	//std::cout<<"test";	 
	}
#endif

}	//for(size_t i = 0; i < genParticles->size(); ++ i) 


#if 1
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

	genJetPt_.push_back( geni.pt());	
	genJetEta_.push_back(geni.eta());	
	genJetPhi_.push_back(geni.phi());	
	genJetVx_.push_back(geni.vx());	
	genJetVy_.push_back(geni.vy());	
	genJetVz_.push_back(geni.vz());	
	genJetCharge_.push_back(geni.charge());	
	genJetMass_.push_back(geni.mass());	
	//std::cout<<"jet px = "<<geni.pt()<<std::endl;
	}

}//if (genjets.isValid()) 

Handle<reco::GenMETCollection> genmets;
iEvent.getByLabel("genMetTrue",genmets);
if (genmets.isValid()) {
	typedef reco::GenMETCollection::const_iterator gmiter;
	for ( gmiter i=genmets->begin(); i!=genmets->end(); i++) {

	genMetPt_.push_back( i->pt());	
	genMetEta_.push_back(i->eta());	
	genMetPhi_.push_back(i->phi());	
	genMetVx_.push_back(i->vx());	
	genMetVy_.push_back(i->vy());	
	genMetVz_.push_back(i->vz());	
	genMetCharge_.push_back(i->charge());	
	genMetMass_.push_back(i->mass());	
	
	}
}

#else


	std::cout<<"Ramkrishna"<<std::endl;
	
#endif

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif

	tree->Fill();

}


// ---------Add branches ---------------------
void GenSimAnalyzer::AddBranch(std::vector<double>* vec, std::string name){
	tree->Branch(name.c_str(),vec);
}

void GenSimAnalyzer::SetBranches(){
	AddBranch(&genLepPt_,	"genLepPt");
	AddBranch(&genLepEta_,	"genLepEta");
	AddBranch(&genLepPhi_,	"genLepPhi");
	AddBranch(&genJetPt_,	"genJetPt");
	AddBranch(&genJetEta_,	"genJetEta");
	AddBranch(&genJetPhi_,	"genJetPhi");
	AddBranch(&genMetPt_,	"genMetPt");
	AddBranch(&genMetEta_,	"genMetEta");
	AddBranch(&genMetPhi_,	"genMetPhi");
}

void GenSimAnalyzer::Clear(){
	genLepPt_.clear();
	genLepEta_.clear();
	genLepPhi_.clear();
	genJetPt_.clear();
	genJetEta_.clear();
	genJetPhi_.clear();
	genMetPt_.clear();
	genMetEta_.clear();
	genMetPhi_.clear();
}
// ------------ method called once each job just before starting event loop  ------------
void 
GenSimAnalyzer::beginJob()
{

	file = new TFile("outfile.root","recreate");    
	tree = new TTree("tree","tree");


	SetBranches();
   #if 0
	tree  = fs->make<TTree>("tree", "tree");
	std::cout<<"inside begin job"<<std::endl;
	if(wantLocalFile_)
	{
	std::cout<<"inside rootFilename_"<<std::endl;
	lf     = new TFile(loutputFile_.c_str(), "RECREATE");
	tree   = new TTree("myEvent","a tree with histograms");
	}
   #else
   #endif
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenSimAnalyzer::endJob() 
{
	file->Write();
	file->Close();
	#if 0
	if(wantLocalFile_)
	lf->WriteTObject(tree);

	delete tree;
	std::cout<<"Written the root file"<<std::endl;
   #else
   #endif
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
