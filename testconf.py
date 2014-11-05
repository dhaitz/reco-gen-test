import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
    'file:/storage/8/dhaitz/testfiles/rundepMC.root'
    #'file:/storage/8/dhaitz/testfiles/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6--Summer12_DR53X-NoPileUp_START53_V7N-v1--AODSIM.root'
    #'file:/storage/8/dhaitz/testfiles/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6--Summer12_DR53X-PU_S10_START53_V7A-v1--AODSIM.root'
    #'file:/storage/8/dhaitz/testfiles/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6--Summer12_DR53X-PU_S10_START53_V7A-v1--AODSIM.root'
    #'file:/storage/8/dhaitz/testfiles/TTJets_FullLeptMGDecays_8TeV-madgraph__Summer12_DR53X-PU_RD1_START53_V7N-v1__AODSIM.root'
    #'root://xrootd.unl.edu//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_RD1_START53_V7N-v1/20000/024197EF-59CF-E211-9DD3-001E67397CB0.root'
))

process.demo = cms.EDAnalyzer('DemoAnalyzer')

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('recogen.root'))

process.p = cms.Path(process.demo)
