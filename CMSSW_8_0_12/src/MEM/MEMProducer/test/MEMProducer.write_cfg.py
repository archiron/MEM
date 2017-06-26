import FWCore.ParameterSet.Config as cms
import os, sys

if len(sys.argv) > 1:
#    print "master step 2 - arg. 0 :", sys.argv[0]
#    print "master step 2 - arg. 1 :", sys.argv[1]
#    print "master step 2 - arg. 2 :", sys.argv[2]
    max_number = int(sys.argv[2])
if max_number == 0:
    max_number = -1
print "max number of evts to treat : ", max_number

process = cms.Process("OWNPARTICLES")

# MET TESTS TOOLS
try: IsMC
except NameError:
    IsMC=True
#try: USE_NOHFMET
#except NameError:
#    USE_NOHFMET=False
#PFMetName = "slimmedMETs"
#if USE_NOHFMET: PFMetName = "slimmedMETsNoHF"
#try: APPLYMETCORR
#except NameError:
#    APPLYMETCORR=True
#METfiltersProcess = "PAT" if IsMC else "RECO"
# END MET TESTS TOOLS

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load("Configuration.Geometry.GeometryIdeal_cff") # new
process.load("Configuration.StandardSequences.MagneticField_cff") # new

process.MessageLogger = cms.Service("MessageLogger",
        destinations = cms.untracked.vstring(                           
                                               'critical_ID_0',
                                             'detailedInfo_ID_0',
                                               'warni_ID_0'
        ),
       critical_ID_0       = cms.untracked.PSet(
                       threshold = cms.untracked.string('ERROR') 
        ),
       detailedInfo_ID_0   = cms.untracked.PSet(
                       threshold  = cms.untracked.string('INFO') 
       ),
       warni_ID_0           = cms.untracked.PSet(
                       threshold  = cms.untracked.string('WARNING') 
        )
)                                                                       

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(max_number) )
max_skipped = 4 # 4 for ttH/00CD6E0F-003B-E611-9817-002590DE6E8A.root
max_skipped = 0
print "max number of evts to skip : ", max_skipped

process.source = cms.Source("PoolSource", 
skipEvents = cms.untracked.uint32(max_skipped),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring([
#        'file:/afs/cern.ch/cms/Tutorials/TWIKI_DATA/CMSDataAnaSch_RelValZMM536.root',
#        '/store/relval/CMSSW_7_6_2/RelValZEE_13/MINIAODSIM/76X_mcRun2_asymptotic_v12-v1/00000/C6A84240-999C-E511-8B1F-0CC47A4C8F0C.root',
#        '/store/relval/CMSSW_7_6_2/RelValZEE_13/MINIAODSIM/76X_mcRun2_asymptotic_v12-v1/00000/EC3E1E3F-999C-E511-AAA7-002618943843.root',
#        '/store/relval/CMSSW_7_6_2/RelValZEE_13/GEN-SIM-RECO/76X_mcRun2_asymptotic_v12-v1/00000/023B475F-989C-E511-99B5-0CC47A78A4B8.root',
#        '/store/relval/CMSSW_7_6_2/RelValZEE_13/GEN-SIM-RECO/76X_mcRun2_asymptotic_v12-v1/00000/0A693408-9A9C-E511-A1D7-0CC47A4C8E16.root',
#        '/store/relval/CMSSW_7_6_2/RelValZEE_13/GEN-SIM-RECO/76X_mcRun2_asymptotic_v12-v1/00000/12515D08-9A9C-E511-9CF4-0CC47A4C8F18.root'
# files below are for GEN-SIM-RECO
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/12515D08-9A9C-E511-9CF4-0CC47A4C8F18.root',
# files below are for MINIAODSIM
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/MiniAOD_800pre4/A62986EB-CCA5-E511-8121-0CC47A745294.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/MiniAOD_800pre4/F6EBA5EB-CCA5-E511-897F-0026189438BC.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/MiniAOD_760pre2/CC91629B-C736-E511-B8DB-003048FFD744.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/MiniAOD_760pre2/FECC1FAD-C736-E511-B533-00261894389A.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/RunIIFall15MiniAODv2/1AC87B8F-8BB8-E511-92F5-FA163E84A67A.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/RunIIFall15MiniAODv2/2EEEFC79-8BB8-E511-BD3B-FA163E4E4AF3.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/RunIIFall15MiniAODv2/3CD84692-8BB8-E511-B6D2-FA163EC03C57.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/RunIIFall15MiniAODv2/4AAC498F-8BB8-E511-A9E0-FA163E84A67A.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/RunIIFall15MiniAODv2/4C5D8578-8BB8-E511-A7CF-FA163EEB25F4.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/RunIISpring16MiniAODv2/VBF/00E98807-AE54-E611-9D79-3417EBE47C5E.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/RunIISpring16MiniAODv2/VBF/1AC87B8F-8BB8-E511-92F5-FA163E84A67A.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/RunIISpring16MiniAODv2/VBF_M130/6E856ED9-F6BD-E511-8361-AC853D9DACD3.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/RunIISpring16MiniAODv2/ttH/00CD6E0F-003B-E611-9817-002590DE6E8A.root',
#        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/ThomasRootFiles/D41D878F-69BE-E611-8B73-0CC47A4C8E98.root',
        'file:/grid_mnt/vol__vol_U__u/llr/info/chiron_u/ANALYZE/ThomasRootFiles/pickevents.root',

    ])
)

process.load("MEM.MEMProducer.MEMParams_cfi")
process.load("MEM.MEMProducer.METCorr_cfi")

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

# MET TESTS TOOLS
process.METSequence = cms.Sequence()

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(
	process,
	isData = (not IsMC),
)
# MET new version
#process.METSignificance2 = cms.EDProducer ("ExtractMETSignificance",
#    mets2=cms.InputTag("slimmedMETs","","OWNPARTICLES"),

#    srcSig = cms.InputTag("METSignificance", "METSignificance"),
#    srcCov = cms.InputTag("METSignificance", "METCovariance"),
#)
process.METSequence += process.fullPatMetSequence
process.METSequence += process.METSignificance2 

# MET old version
#process.load("RecoMET.METProducers.METSignificance_cfi") 
#process.load("RecoMET.METProducers.METSignificanceParams_cfi") 

#process.METSequence += process.METSignificance 
# END MET TESTS TOOLS

process.load("MEM.MEMProducer.MEMProducer_cfi")

process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('myOutputFile.ID_0.root'),
    fastCloning = cms.untracked.bool(False),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_MEMProducer_*_*', 
        'keep patElectrons_*_*_*',
        'keep patMuons_*_*_*',
        'keep patMETs_*_*_*',
        'keep patTaus_*_*_*',
        'keep patJets_*_*_*',
##        'keep patElectrons_MEMProducer_*_*',
#        'keep patMuons_MEMProducer_*_*',
#        'keep patMETs_MEMProducer_*_*',
      )

)
  
process.p = cms.Path(
#    process.METSignificance2 +
    process.METSequence #+ 
#    process.MEMProducer
    ) 

process.e = cms.EndPath(process.out)
