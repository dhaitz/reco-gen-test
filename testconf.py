import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# defaults and global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START53_V27::All'


process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
    #'file:/storage/8/dhaitz/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball__Summer12_DR53X-PU_RD1_START53_V7N-v1__AODSIM.root'
    #'file:/storage/8/dhaitz/testfiles/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6--Summer12_DR53X-NoPileUp_START53_V7N-v1--AODSIM.root'
    #'file:/storage/8/dhaitz/testfiles/QCD_Pt_30_80_EMEnriched_TuneZ2star_8TeV_pythia6--Summer12_DR53X-PU_S10_START53_V7A-v1--AODSIM.root'
    #'file:/storage/8/dhaitz/testfiles/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6--Summer12_DR53X-PU_S10_START53_V7A-v1--AODSIM.root'
    #'file:/storage/8/dhaitz/testfiles/TTJets_FullLeptMGDecays_8TeV-madgraph__Summer12_DR53X-PU_RD1_START53_V7N-v1__AODSIM.root'
    'file:/storage/8/dhaitz/testfiles/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6__Summer12_DR53X-PU_RD1_START53_V7N-v1__AODSIM.root'
))



# corrected jets
process.load('JetMETCorrections.Configuration.JetCorrectionProducers_cff')
process.load('JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff')

process.ak5PFL1L2L3 = cms.ESProducer(
    'JetCorrectionESChain',
    correctors = cms.vstring('ak5PFL1Fastjet','ak5PFL2Relative','ak5PFL3Absolute')
)
process.ak5PFJetsL1 = cms.EDProducer(
    'PFJetCorrectionProducer',
    src        = cms.InputTag('ak5PFJets'),
    correctors = cms.vstring('ak5PFL1Fastjet')
)

process.jets = cms.Sequence(process.ak5PFJetsL1L2L3 * process.ak5PFJetsL1)



process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.NoNuGenJets = cms.Path(
    process.genParticlesForJetsNoNu *
    process.genParticlesForJets *
    process.ak5GenJetsNoNu)



# analyzer and output
process.l1l2l3 = cms.EDAnalyzer('DemoAnalyzer',
    jetCollection = cms.InputTag("ak5PFJetsL1L2L3")
)
process.l1 = cms.EDAnalyzer('DemoAnalyzer',
    jetCollection = cms.InputTag("ak5PFJetsL1")
)
process.raw = cms.EDAnalyzer('DemoAnalyzer',
    jetCollection = cms.InputTag("ak5PFJets")
)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string('recogen.root')
)

process.p = cms.Path(
    process.jets *
    process.l1 *
    process.l1l2l3 *
    process.raw
)

process.schedule = cms.Schedule(
process.NoNuGenJets,
    process.p
)
