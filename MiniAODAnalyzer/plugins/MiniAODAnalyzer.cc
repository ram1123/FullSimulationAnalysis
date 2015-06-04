// -*- C++ -*-
//
// Package:    FullSimAnalysis/MiniAODAnalyzer
// Class:      MiniAODAnalyzer
// 
/**\class MiniAODAnalyzer MiniAODAnalyzer.cc FullSimAnalysis/MiniAODAnalyzer/plugins/MiniAODAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Ram Krishna Sharma
//         Created:  Thu, 04 Jun 2015 18:30:14 GMT
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

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//
// class declaration
//

class MiniAODAnalyzer : public edm::EDAnalyzer {
   public:
      explicit MiniAODAnalyzer(const edm::ParameterSet&);
      ~MiniAODAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;
      edm::EDGetTokenT<pat::ElectronCollection> electronToken_;
      edm::EDGetTokenT<pat::TauCollection> tauToken_;
      edm::EDGetTokenT<pat::PhotonCollection> photonToken_;
      edm::EDGetTokenT<pat::JetCollection> jetToken_;
      edm::EDGetTokenT<pat::JetCollection> fatjetToken_;
      edm::EDGetTokenT<pat::METCollection> metToken_;
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
MiniAODAnalyzer::MiniAODAnalyzer(const edm::ParameterSet& iConfig):
    vtxToken_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
    muonToken_(consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
    electronToken_(consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons"))),
    tauToken_(consumes<pat::TauCollection>(iConfig.getParameter<edm::InputTag>("taus"))),
    photonToken_(consumes<pat::PhotonCollection>(iConfig.getParameter<edm::InputTag>("photons"))),
    jetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jets"))),
    fatjetToken_(consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("fatjets"))),
    metToken_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("mets")))
{
   //now do what ever initialization is needed

}


MiniAODAnalyzer::~MiniAODAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MiniAODAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
    edm::Handle<reco::VertexCollection> vertices;
    iEvent.getByToken(vtxToken_, vertices);
    if (vertices->empty()) return; // skip the event if no PV found
    const reco::Vertex &PV = vertices->front();

    edm::Handle<pat::MuonCollection> muons;
    iEvent.getByToken(muonToken_, muons);
    for (const pat::Muon &mu : *muons) {
        if (mu.pt() < 5 || !mu.isLooseMuon()) continue;
        printf("muon with pt %4.1f, dz(PV) %+5.3f, POG loose id %d, tight id %d\n",
                mu.pt(), mu.muonBestTrack()->dz(PV.position()), mu.isLooseMuon(), mu.isTightMuon(PV));
    }

    edm::Handle<pat::ElectronCollection> electrons;
    iEvent.getByToken(electronToken_, electrons);
    for (const pat::Electron &el : *electrons) {
        if (el.pt() < 5) continue;
        printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), pass conv veto %d\n",
        //printf("elec with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes), lost hits %d, pass conv veto %d\n",
                    el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.passConversionVeto());
                    //el.pt(), el.superCluster()->eta(), el.sigmaIetaIeta(), el.full5x5_sigmaIetaIeta(), el.gsfTrack()->trackerExpectedHitsInner().numberOfLostHits(), el.passConversionVeto());
    }

    edm::Handle<pat::PhotonCollection> photons;
    iEvent.getByToken(photonToken_, photons);
    for (const pat::Photon &pho : *photons) {
        if (pho.pt() < 20 or pho.chargedHadronIso()/pho.pt() > 0.3) continue;
        printf("phot with pt %4.1f, supercluster eta %+5.3f, sigmaIetaIeta %.3f (%.3f with full5x5 shower shapes)\n",
                    pho.pt(), pho.superCluster()->eta(), pho.sigmaIetaIeta(), pho.full5x5_sigmaIetaIeta());
    }


    edm::Handle<pat::TauCollection> taus;
    iEvent.getByToken(tauToken_, taus);
    for (const pat::Tau &tau : *taus) {
        if (tau.pt() < 20) continue;
        printf("tau  with pt %4.1f, dxy signif %.1f, ID(byMediumCombinedIsolationDeltaBetaCorr3Hits) %.1f, lead candidate pt %.1f, pdgId %d \n",
                    tau.pt(), tau.dxy_Sig(), tau.tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"), tau.leadCand()->pt(), tau.leadCand()->pdgId());
    }


    edm::Handle<pat::JetCollection> jets;
    iEvent.getByToken(jetToken_, jets);
    int ijet = 0;
    for (const pat::Jet &j : *jets) {
        if (j.pt() < 20) continue;
        printf("jet  with pt %5.1f (raw pt %5.1f), eta %+4.2f, btag CSV %.3f, CISV %.3f, pileup mva disc %+.2f\n",
            j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), std::max(0.f,j.bDiscriminator("combinedSecondaryVertexBJetTags")), std::max(0.f,j.bDiscriminator("combinedInclusiveSecondaryVertexBJetTags")), j.userFloat("pileupJetId:fullDiscriminant"));
        if ((++ijet) == 1) { // for the first jet, let's print the leading constituents
            std::vector<reco::CandidatePtr> daus(j.daughterPtrVector());
            std::sort(daus.begin(), daus.end(), [](const reco::CandidatePtr &p1, const reco::CandidatePtr &p2) { return p1->pt() > p2->pt(); }); // the joys of C++11
            for (unsigned int i2 = 0, n = daus.size(); i2 < n && i2 <= 3; ++i2) {
                const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*daus[i2]);
                printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", i2,cand.pt(),cand.dz(PV.position()),cand.pdgId());
            }
        }
    }


    edm::Handle<pat::JetCollection> fatjets;
    iEvent.getByToken(fatjetToken_, fatjets);
    for (const pat::Jet &j : *fatjets) {
        printf("AK8j with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed, %5.1f softdrop, %5.1f pruned, %5.1f trimmed, %5.1f filtered. CMS TopTagger %.1f\n",
            j.pt(), j.pt()*j.jecFactor("Uncorrected"), j.eta(), j.mass(), j.userFloat("ak8PFJetsCHSSoftDropMass"), j.userFloat("ak8PFJetsCHSPrunedMass"), j.userFloat("ak8PFJetsCHSTrimmedMass"), j.userFloat("ak8PFJetsCHSFilteredMass"), j.userFloat("cmsTopTagPFJetsCHSMassAK8"));




        // To get the constituents of the AK8 jets, you have to loop over all of the
        // daughters recursively. To save space, the first two constituents are actually
        // the Soft Drop SUBJETS, which will then point to their daughters.
        // The remaining constituents are those constituents removed by soft drop but
        // still in the AK8 jet.
	std::vector<reco::Candidate const *> constituents;
        for ( unsigned ida = 0; ida < j.numberOfDaughters(); ++ida ) {
	  reco::Candidate const * cand = j.daughter(ida);
	  if ( cand->numberOfDaughters() == 0 )
	    constituents.push_back( cand ) ;
	  else {
	    for ( unsigned jda = 0; jda < cand->numberOfDaughters(); ++jda ) {
	      reco::Candidate const * cand2 = cand->daughter(jda);
	      constituents.push_back( cand2 );
	    }
	  }
	}
	std::sort( constituents.begin(), constituents.end(), [] (reco::Candidate const * ida, reco::Candidate const * jda){return ida->pt() > jda->pt();} );

	for ( unsigned int ida = 0; ida < constituents.size(); ++ida ) {
	  const pat::PackedCandidate &cand = dynamic_cast<const pat::PackedCandidate &>(*constituents[ida]);
	  printf("         constituent %3d: pt %6.2f, dz(pv) %+.3f, pdgId %+3d\n", ida,cand.pt(),cand.dz(PV.position()),cand.pdgId());
	}


	auto wSubjets = j.subjets("SoftDrop");
	for ( auto const & iw : wSubjets ) {
	  printf("   w subjet with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed\n",
		 iw->pt(), iw->pt()*iw->jecFactor("Uncorrected"), iw->eta(), iw->mass() );
	
	}
	auto tSubjets = j.subjets("CMSTopTag");
	for ( auto const & it : tSubjets ) {
	  printf("   t subjet with pt %5.1f (raw pt %5.1f), eta %+4.2f, mass %5.1f ungroomed\n",
		 it->pt(), it->pt()*it->jecFactor("Uncorrected"), it->eta(), it->mass() );
	
	}
    }
 
    edm::Handle<pat::METCollection> mets;
    iEvent.getByToken(metToken_, mets);
    const pat::MET &met = mets->front();
    printf("MET: pt %5.1f, phi %+4.2f, sumEt (%.1f). genMET %.1f. MET with JES up/down: %.1f/%.1f\n",
        met.pt(), met.phi(), met.sumEt(),
        met.genMET()->pt(),
        met.shiftedPt(pat::MET::JetEnUp), met.shiftedPt(pat::MET::JetEnDown));

    printf("\n");

}


// ------------ method called once each job just before starting event loop  ------------
void 
MiniAODAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MiniAODAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
MiniAODAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
MiniAODAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
MiniAODAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
MiniAODAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniAODAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniAODAnalyzer);
