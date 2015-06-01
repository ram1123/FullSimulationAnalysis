import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/r/rasharma/work/WW_Scattering/AnalyzerFiles/GEN-SIM_Analysis/CMSSW_7_1_16/src/GEN-SIM/GenSimAnalyzer/qqToqqWWTolnu4q-RunIIWinter15GS_LL_1.root'
    )
)

process.demo = cms.EDAnalyzer('GenSimAnalyzer',
	
	#################       InputTags       #################################
	GenPart         =       cms.InputTag("genParticles"),
	GenJetAk4       =       cms.InputTag("ak4GenJets"),
	GenMetTru       =       cms.InputTag("genMetTrue"),
	GenMetCal       =       cms.InputTag("genMetCalo"),
	trigTag         =       cms.InputTag("TriggerResults","","SIM"),

	#################       Some Booleans   #################################

	Verbose         =       cms.untracked.bool(True),
	SortGen         =       cms.untracked.bool(True),
	wantLocalFile   =       cms.untracked.int32(0),
	wantRFIOFile    =       cms.untracked.int32(1),
	loutputFile     =       cms.untracked.string("RawSim.root"),
	rfoutputFile    =       cms.untracked.string("/tmp/rasharma/GenSim.root")
)


process.p = cms.Path(process.demo)
