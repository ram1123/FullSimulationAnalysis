import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        #'file:myfile.root'
	'file:/afs/cern.ch/user/r/rasharma/work/WW_Scattering/CMSFullSimulation/CMSSW_7_4_1_patch1/src/B2G-RunIISpring15DR74-00001.root'
	#'/store/cmst3/user/gpetrucc/miniAOD/v1/TT_Tune4C_13TeV-pythia8-tauola_PU_S14_PAT.root'
	#'root://xrootd.unl.edu//store/mc/Phys14DR/WJetsToLNu_13TeV-madgraph-pythia8-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/02215B44-2D70-E411-90A3-0025905A60B8.root'
    )
)

process.demo = cms.EDAnalyzer('MiniAODAnalyzer',

	#################	InputTags	#################################
	vertices	= 	cms.InputTag("offlineSlimmedPrimaryVertices"),
	muons		= 	cms.InputTag("slimmedMuons"),
	electrons	= 	cms.InputTag("slimmedElectrons"),
	taus		= 	cms.InputTag("slimmedTaus"),
	photons		= 	cms.InputTag("slimmedPhotons"),
	jets		= 	cms.InputTag("slimmedJets"),
	fatjets		= 	cms.InputTag("slimmedJetsAK8"),
	mets		= 	cms.InputTag("slimmedMETs"),

	#################	Some Booleans	#################################
	Verbose		=	cms.untracked.bool(False),
	wantLocalFile	=	cms.untracked.int32(0),
	wantRFIOFile	=	cms.untracked.int32(0),
	loutputFile	=	cms.untracked.string("RawSim.root"),
	rfoutputFile	=	cms.untracked.string("/tmp/rasharma/RawSim.root")
)


process.p = cms.Path(process.demo)
