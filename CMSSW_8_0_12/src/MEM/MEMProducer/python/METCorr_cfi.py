import FWCore.ParameterSet.Config as cms

from MEM.MEMProducer.MEMParams_cfi import *

METSignificance2 = cms.EDProducer ("ExtractMETSignificance",
    mets2=cms.InputTag("slimmedMETs"),

    srcSig = cms.InputTag("METSignificance", "METSignificance"),
    srcCov = cms.InputTag("METSignificance", "METCovariance"),

    configName = cms.string('ExtractMETSignificance.py'),

    parameters = cms.PSet(MEMParams)
)
