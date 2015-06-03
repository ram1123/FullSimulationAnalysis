/*
 * =====================================================================================
 *
 *       Filename:  RawSimAnalyzer.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  06/02/15 14:55:20
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Ramkrishna Sharma (Ram), ramkrishna.sharma71@gmail.com
 *   Organization:  University Of Delhi, Delhi, India
 *
 * =====================================================================================
 */


#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include "TTree.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

using namespace edm;
using namespace std;
using namespace reco;
using namespace ROOT::Math::VectorUtil;


class RawSimAnalyzer : public edm::EDAnalyzer {
   public:
      explicit RawSimAnalyzer(const edm::ParameterSet&);
      ~RawSimAnalyzer();

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
      edm::InputTag PileUp_;
      edm::InputTag GenJetAk4_;
      edm::InputTag GenJetsAk4NoNu_;
      edm::InputTag GenMetTru_;
      edm::InputTag GenMetCal_;
      edm::InputTag trigTag_;
      edm::InputTag triggerEventTag_;
//      edm::InputTag genMetCNP_;

      // Few Boolean InputTags
      bool	Verbose_;
      bool	SortGen_;
      bool	runPileupinfo_;
      bool	wantLocalFile_;	
      bool	wantRFIOFile_;	
      std::string	loutputFile_;
      std::string	rfoutputFile_;


      edm::Service<TFileService> fs;
      TFile * outputFile_;
      TTree* tree;


      // Event Properties
      int runNumber;
      int eventNumber;
      int orbitNumber;
      int bunchCrossing;
      int lumiBlockNumber;
      int isData;

      // Basic Trigger Information
      std::vector<std::string> TriggerNames_,*HLTNamesSet;
      std::vector<Int_t> HLTTableBitMatch_;
      std::string *HLTTableName;
      std::vector<UInt_t> *HLTPrescaleFactors;
      std::vector<Bool_t> *HLTriggerResults;
 
      // Pile up Info
      std::vector<double> pileup_bunchXing_	;
      std::vector<double> numberOfPUVertices_	;	// the number of pileup interactions that have been added to the event from this bunch crossing
      std::vector<double> numberOfPUVerticesMixingTruth_	;	// the *true* mean number of pileup interactions for this event from which each bunch crossing has been sampled; same for all bunch crossings in an event


      // Basics GenParticles information
      std::vector<double> GSLepVx_	;
      std::vector<double> GSLepVy_	;
      std::vector<double> GSLepVz_	;
      std::vector<double> GSLepCharge_	;
      std::vector<double> GSLepMass_	;
      std::vector<double> GSLepPt_	;
      std::vector<double> GSLepEta_	;
      std::vector<double> GSLepPhi_	;
      std::vector<double> GSLepTheta_	;
      std::vector<double> GSLepStatus_	;
      std::vector<double> GSLepMother_	;
      std::vector<double> GSLepEle_	;
      std::vector<double> GSLepMuon_	;


      // Basic Jet informations
      int		  nJets		;

      //math::TLorentzVector jetP4_	;
      std::vector<double> GSJetPt_	;
      std::vector<double> GSJetEta_	;
      std::vector<double> GSJetPhi_	;
      std::vector<double> GSJetEMFraction_	;
      std::vector<double> GSJetChargedEmEnergyFraction_	;
      std::vector<double> GSJetNeutralEmEnergyFraction_	;
      std::vector<double> GSJetChargedHadronEnergyFraction_	;
      std::vector<double> GSJetNeutralHadronEnergyFraction_	;
      std::vector<double> GSJetChargedMultiplicity_	;
      std::vector<double> GSJetMass_	;
      std::vector<int>	  GSJetnConstituents_	;
      std::vector<int>	  GSJetnTracks_	;
      math::XYZVector 	  GSJetVertex_	;
      std::vector<double> GSJetVertexChi2_	;
      std::vector<double> GSJetVertexChi2Ndof_	;
      std::vector<double> GSJetVertexNormalizedChi2_	;
      std::vector<double> GSJetVx_	;
      std::vector<double> GSJetVy_	;
      std::vector<double> GSJetVz_	;
      std::vector<double> GSJetCharge_	;
      std::vector<double> GSJetTheta_	;


      int		  nJetsNoNu	;
      std::vector<double> GSJetNoNuVx_	;
      std::vector<double> GSJetNoNuVy_	;
      std::vector<double> GSJetNoNuVz_	;
      std::vector<double> GSJetNoNuCharge_	;
      std::vector<double> GSJetNoNuMass_	;
      std::vector<double> GSJetNoNuPt_	;
      std::vector<double> GSJetNoNuEta_	;
      std::vector<double> GSJetNoNuPhi_	;
      std::vector<double> GSJetNoNuTheta_	;

      std::vector<double> GSMetTruVx_	;
      std::vector<double> GSMetTruVy_	;
      std::vector<double> GSMetTruVz_	;
      std::vector<double> GSMetTruCharge_	;
      std::vector<double> GSMetTruMass_	;
      std::vector<double> GSMetTruPt_	;
      std::vector<double> GSMetTruEta_	;
      std::vector<double> GSMetTruPhi_	;
      std::vector<double> GSMetTruTheta_	;

      std::vector<double> GSMetCalVx_	;
      std::vector<double> GSMetCalVy_	;
      std::vector<double> GSMetCalVz_	;
      std::vector<double> GSMetCalCharge_	;
      std::vector<double> GSMetCalMass_	;
      std::vector<double> GSMetCalPt_	;
      std::vector<double> GSMetCalEta_	;
      std::vector<double> GSMetCalPhi_	;
      std::vector<double> GSMetCalTheta_	;

      std::vector<double> GSMetCNPVx_	;
      std::vector<double> GSMetCNPVy_	;
      std::vector<double> GSMetCNPVz_	;
      std::vector<double> GSMetCNPCharge_	;
      std::vector<double> GSMetCNPMass_	;
      std::vector<double> GSMetCNPPt_	;
      std::vector<double> GSMetCNPEta_	;
      std::vector<double> GSMetCNPPhi_	;
      std::vector<double> GSMetCNPTheta_	;


      std::vector<double> GSMetTrueVx_	;
      std::vector<double> GSMetTrueVy_	;
      std::vector<double> GSMetTrueVz_	;
      std::vector<double> GSMetTrueCharge_	;
      std::vector<double> GSMetTrueMass_	;
      std::vector<double> GSMetTruePt_	;
      std::vector<double> GSMetTrueEta_	;
      std::vector<double> GSMetTruePhi_	;
      std::vector<double> GSMetTrueTheta_	;


};

class PtGreater {
         public:
         template <typename T> bool operator () (const T& i, const T& j) {
         return (i->pt() > j->pt());
         }
};

