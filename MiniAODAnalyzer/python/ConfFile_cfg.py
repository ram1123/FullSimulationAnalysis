import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/r/rasharma/work/WW_Scattering/CMSFullSimulation/CMSSW_7_4_1_patch1/src/B2G-RunIISpring15DR74-00001.root'
    )
)

hltLabel = 'HLT'

process.demo = cms.EDAnalyzer('MiniAODAnalyzer',

	#################	InputTags	#################################

	GenPart		=	cms.InputTag("genParticles"),
	PileUp		=	cms.InputTag("addPileupInfo"),
	GenJetAk4	=	cms.InputTag("ak4GenJets"),
	GenJetsAk4NoNu	=	cms.InputTag("ak4GenJetsNoNu"),
	GenMetTru	=	cms.InputTag("genMetTrue"),
	GenMetCal	=	cms.InputTag("genMetCalo"),
	trigTag		=	cms.InputTag("TriggerResults","","HLT"),	# In place of HLT we can put also SIM, which will give triggers at simulation level
	#trigTag	=	cms.InputTag("TriggerResults::HLT"),
	triggerEventTag	=	cms.InputTag("hltTriggerSummaryAOD","",hltLabel),
	#GenMetCNP	=	cms.InputTag("genMetCaloAndNonPrompt"),

	#################	Some Booleans	#################################

	Verbose		=	cms.untracked.bool(False),
	SortGen		=	cms.untracked.bool(False),
	#runPileupinfo	=	cms.untracked.bool(True),
	wantLocalFile	=	cms.untracked.int32(0),
	wantRFIOFile	=	cms.untracked.int32(1),
	loutputFile	=	cms.untracked.string("RawSim.root"),
	rfoutputFile	=	cms.untracked.string("/tmp/rasharma/RawSim.root")
)

#process.Tracer = cms.Service("Tracer")

process.p = cms.Path(process.demo)
