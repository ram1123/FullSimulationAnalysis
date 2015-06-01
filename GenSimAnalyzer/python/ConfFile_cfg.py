import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
        'file:/eos/uscms/store/user/rasharma/WWScattering/Full_Simulation_Data/Gen_Sim/qqToqqWlWl/GEN-SIM_qqToqqWlWl_NoPtEtaCut/150222_110947/0000/qqToqqWlWl_RunII_13TeV_NoPtEtaCut_97.root'
    )
)

process.demo = cms.EDAnalyzer('GenSimAnalyzer'
)


process.p = cms.Path(process.demo)
