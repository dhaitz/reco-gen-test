import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

# configure input file
opt = VarParsing ('analysis')
opt.register ('sample', 'dy', VarParsing.multiplicity.singleton,
            VarParsing.varType.string, "qcd or dy")
opt.parseArguments()
leptonIso = (opt.sample == 'dy')

files_dict = {
    'dy':[
        "root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_RD1_START53_V7N-v1/20000/A25A61FC-C9CE-E211-A4B9-003048D45FEE.root",
        "root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_RD1_START53_V7N-v1/20001/749933D9-8AD0-E211-AF0E-003048D46060.root",
        "root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_RD1_START53_V7N-v1/20001/4612855C-48D0-E211-A792-001E673983F4.root",
        "root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_RD1_START53_V7N-v1/20001/92C7EEBF-60D0-E211-87CB-001E673983F4.root",
        "root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_RD1_START53_V7N-v1/20001/D8D9B0BE-6FD0-E211-A11F-001E67396D56.root",
        "root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_RD1_START53_V7N-v1/20002/AA66CC97-9DD1-E211-B8DA-00259020081C.root",
        #"file:/storage/8/dhaitz/testfiles/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM.root"
    ],
    'qcd':[
        #"root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_RD1_START53_V7N-v1/20000/2A762626-836E-E311-B825-00266CFFA2B8.root",
        #"root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_RD1_START53_V7N-v1/20000/E0018766-756C-E311-9522-0025901D4D64.root"
        #"root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6/AODSIM/PU_RD1_START53_V7N-v1/20000/ECF7B78F-736C-E311-8AA4-00266CF25320.root",
        "file:/storage/8/dhaitz/testfiles/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM.root"
    ],
    #'qcd':'file:/storage/8/dhaitz/testfiles/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM.root',
    #'qcdnopu':'file:/storage/8/dhaitz/testfiles/QCD_Pt-15to3000_TuneZ2star_Flat_8TeV_pythia6--Summer12_DR53X-NoPileUp_START53_V7N-v1--AODSIM.root',
    #'dymm': 'file:/storage/8/dhaitz/DYToMuMu_M-20_CT10_TuneZ2star_v2_8TeV-powheg-pythia6__Summer12_DR53X-PU_RD1_START53_V7N-v1__AODSIM.root',
    #'tt': 'file:/storage/8/dhaitz/testfiles/TTJets_FullLeptMGDecays_8TeV-madgraph__Summer12_DR53X-PU_RD1_START53_V7N-v1__AODSIM.root',
    #'powheg': 'file:/storage/8/dhaitz/testfiles/DYToEE_M-20_CT10_TuneZ2star_8TeV-powheg-pythia6--Summer12_DR53X-PU_S10_START53_V7A-v1--AODSIM.root',
    #'ww': 'file:/storage/8/dhaitz/WW_TuneZ2star_8TeV_pythia6_tauola__Summer12_DR53X-PU_S10_START53_V7A-v1__AODSIM.root',
    #'qcd3050': 'file:/storage/8/dhaitz/testfiles/QCD_Pt-30to50_TuneZ2star_8TeV_pythia6__Summer12_DR53X-PU_S10_START53_V7A-v2__AODSIM.root',
}
#small:
#"root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_RD1_START53_V7N-v1/20002/3ECBF087-CCD1-E211-BAA6-001E67397AE4.root",
#"root://cms-xrd-global.cern.ch//store/mc/Summer12_DR53X/DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball/AODSIM/PU_RD1_START53_V7N-v1/20001/221CDD64-A8CF-E211-A4C4-001E673968A6.root",

process = cms.Process("Demo")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# Reduce amount of messages -----------------------------------------------
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.default = cms.untracked.PSet(
    ERROR = cms.untracked.PSet(limit = cms.untracked.int32(5)))
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

# defaults and global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'START53_V27::All'


process.source = cms.Source("PoolSource", fileNames = cms.untracked.vstring(
  *files_dict[opt.sample]
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

# gen jets
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')
process.NoNuGenJets = cms.Path(
    process.genParticlesForJetsNoNu *
    process.genParticlesForJets *
    process.ak5GenJetsNoNu)


# analyzer and output
process.raw = cms.EDAnalyzer('DemoAnalyzer',
    jetCollection = cms.InputTag("ak5PFJets"),
    leptonIso = cms.untracked.bool(leptonIso),
    folder = cms.InputTag("raw")
)
process.l1 = process.raw.clone(
    jetCollection = cms.InputTag("ak5PFJetsL1"),
    folder = cms.InputTag("l1")
)
process.l1l2l3 = process.l1.clone(
    jetCollection = cms.InputTag("ak5PFJetsL1L2L3"),
    folder = cms.InputTag("l1l2l3")
)

# output file
process.TFileService = cms.Service("TFileService",
        fileName = cms.string('recogen_%s.root' % opt.sample)
)

# path, schedule
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
