import FWCore.ParameterSet.Config as cms

from MEM.MEMProducer.MEMParams_cfi import *

MEMProducer = cms.EDProducer('MEMProducer',
    Verbosity = cms.untracked.int32(0),
    
    electrons = cms.InputTag('slimmedElectrons'),
    muons = cms.InputTag('slimmedMuons'),
    taus = cms.InputTag('slimmedTaus'),
    mets = cms.InputTag('slimmedMETs'),
    jets = cms.InputTag('slimmedJets'),
    
    mets2 = cms.InputTag('slimmedMETs'),

    srcSig = cms.InputTag("METSignificance", "METSignificance"),
    srcCov = cms.InputTag("METSignificance", "METCovariance"),
    
    configName = cms.string('MEMProducer.py'),

    parameters = cms.PSet(MEMParams)
)
