import FWCore.ParameterSet.Config as cms
import sys
import os
import math

# MEM Base
MEM_BASE=os.getcwd()
#MEM_BASE=os.environ['MEM_BASE']
CMSSWBASE = os.environ['CMSSW_BASE'] # donne le repertoire de travail
CMSSWBASECMSSWRELEASEBASE = os.environ['CMSSW_RELEASE_BASE'] # donne la release et l'architecture
CMSSWBASECMSSWVERSION = os.environ['CMSSW_VERSION'] # donne la release (CMSSW_7_1_0 par exemple)
print "*** MEM_BASE                  : %s" % MEM_BASE
print "*** CMSSWBASE                 : %s" % CMSSWBASE
print "*** CMSSWBASECMSSWRELEASEBASE : %s" % CMSSWBASECMSSWRELEASEBASE
print "*** CMSSWBASECMSSWVERSION     : %s" % CMSSWBASECMSSWVERSION
print "*** ttH VERSION ***"

# Math constants
PI          = math.pi

# LHC beam energy (Gev)
sqrtS = 13000

# PDF calibration (Gev)
# Unused ??? to delete (replaced with mHiggs or mZ )
Q = 125

# Particles (Gev)
# Local variable for python interface
# ???
mTau    = 1.776
mHiggs  = 125.09
mZ      = 91.187

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Input files
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Input Type : MonteCarlo/Run1
# InputType = "MonteCarlo"
# InputType = "Run1"
# InputFileList = [ "../../Data/ntuple_VBF_8TeV_JEC_METCorrtype1_11.root"]
# InputFileList = ["../../Data/ntuple_DYtestX_8TeV_JEC_1.root"]
# InputFileList = ["/data_CMS/cms/htautau/MEM/data/test/ntuple_DYtestX_8TeV_JEC_1.root"]
# Input Type : Run2
InputType = "PyRun2"

#ttH
InputFileList = ["/data_CMS/cms/strebler/ttH_prod_80X_09_2016/ntuples_splitted/nominal/ttH/sync_HTauTauTree_split_80X.root"]

#ttbar SLfromT
#InputFileList = ["/data_CMS/cms/strebler/ttH_prod_80X_06_2016/ntuples_splitted/TTJets/HTauTauTree_ttbar_SLfromT_split.root"]

#ttbar SLfromTbar
#InputFileList = ["/data_CMS/cms/strebler/ttH_prod_80X_06_2016/ntuples_splitted/TTJets/HTauTauTree_ttbar_SLfromTbar_split.root"]

#ttbar DL
#InputFileList = [" /data_CMS/cms/strebler/ttH_prod_80X_06_2016/ntuples_splitted/TTJets/HTauTauTree_ttbar_DL_split.root"]

treeName = "HTauTauTree_2lSS"
#treeName = "HTauTauTree_2lSS_lepMVA_CR"

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Events
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ttH Htautau 2lSS
#maxNbrOfEventsToRead = 117719
#maxNbrOfEventsToRead = 10000
maxNbrOfEventsToRead = 50

#ttbar DL 2lSS lepMVA CR
#maxNbrOfEventsToRead = 1257

#ttbar SLfromT 2lSS lepMVA CR
#maxNbrOfEventsToRead = 5303

#ttbar SLfromTbar 2lSS lepMVA CR
#maxNbrOfEventsToRead = 6882

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Random Number Generator
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
flagSameRNG = False

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Computation / Run
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# ???
flagTFLepTau = True 
flagTFHadTau = True
flagTFFake = True
flagTFMET    = True
flagTFJet1   = True
flagTFJet2   = True
flagTFBJet_leptop = True
flagTFBJet_hadtop = True
flagTF_fakelep = True
flagTF_fakeleptau = True
flagTFTop    = True
flagJac      = True
flagWME      = True

#  Verbose Modes
NoVerbose        = 0
ResultLevel      = 1 
IntegrationLevel = 2 
IntegrandLevel   = 3 
# Seclect verbose mode
verbose      = IntegrandLevel

#
#  ME version 1=full ME, 2=v2, 3=v3
MEversion = 1

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Output files
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
FileOfIntegrals = "ttH_Hnonbb_split_2lSS_MEM_160704"
#FileOfFunctionValues = ""
#FileOfFunctionValues = "fct-values-all-process-dy"
# FileOfJetValues      = "jet-values"
#

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# MC Integration
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
nbrOfDimttH    = 5
nbrOfDimttZ    = 5
nbrOfDimttW    = 7
nbrOfDimttjets = 4
nbrOfDimttbar_SL = 3
nbrOfDimttbar_DL = 5
nbrOfDimttZ_Zll    = 3

nbrOfDimttH_miss    = 7
nbrOfDimttZ_miss    = 7
nbrOfDimttjets_miss = 6
nbrOfDimttbar_SL_miss = 5
nbrOfDimttZ_Zll_miss    = 5

runttZ_integration = False
runttW_integration = False
runttjets_integration = False
runttbar_SL_integration = False
runttbar_DL_integration = False
runttbar_DL_fakelep_integration = False
runttZ_Zll_integration = True

run_missing_jet_integration = True
force_missing_jet_integration = False
force_missing_jet_integration_ifnoperm = True

#nbrOfPoints_ttZ = 10
nbrOfPoints_ttZ = 10
#nbrOfPoints_ttH = 10000
nbrOfPoints_ttH = 10
nbrOfPoints_ttW = 20
nbrOfPoints_ttjets = 10
#nbrOfPoints_ttbar_SL = 5
nbrOfPoints_ttbar_SL = 5
#nbrOfPoints_ttbar_DL = 10000
nbrOfPoints_ttbar_DL = 10
nbrOfPoints_ttZ_Zll = 5

#nbrOfPoints_ttH_miss = 20000
nbrOfPoints_ttZ_miss = 20
nbrOfPoints_ttH_miss = 20
#nbrOfPoints_ttZ_miss = 20
nbrOfPoints_ttjets_miss = 10
#nbrOfPoints_ttbar_SL_miss = 10000
nbrOfPoints_ttbar_SL_miss = 10
nbrOfPoints_ttZ_Zll_miss = 10

nbrOfPermut_per_jet = 4

# Integration boundaries 
phiNu_tlep  = [ 0, 2*PI ]
cosThetaNu_tlep = [-1.,1.]
phiNu_ttau  = [ 0, 2*PI ]
cosThetaNu_ttau = [-1.,1.]
phiNu_W_ttW  = [ 0, 2*PI ]
cosThetaNu_W_ttW = [-1.,1.]
 
CI_TFJet = 0.95 #Use 95% CI from the jet TF
use_pT_TFJet = True
use_top_compatibility_check = True

include_hadrecoil = True
force_nonzero_integral = False

#Missing jet parameters
eta_acceptance = 2.4
jet_radius = 0.4
dR_veto_jet_lep = 0.4
rel_iso_lep = 0.4
pT_cut = 25.

phi_missing_jet = [ 0, 2*PI ]
cosTheta_missing_jet = [-1.,1.]

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# MadGraph 5
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#MGParamCardFile = MEM_BASE + "/MGProcesses/param_card.dat"
MGParamCardFile = "/home/llr/info/chiron_u/CMSSW_8_0_12/src/MEM/MEMProducer/test/param_card.dat"

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# LHAPDF
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# LHAPDFFile = MEM_BASE + "/MEMData/cteq65.LHgrid"
LHAPDFFile = "cteq66"

# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# LUTs
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#LUTOneProngPi0File = MEM_BASE + "/MGMEM/h_1P1PZeroSpectrumZ_LUT.lut"
#LUTThreeProngsFile = MEM_BASE + "/MGMEM/h_3PSpectrumZ_LUT.lut"


MEMParams = cms.PSet(

Q = cms.double(Q),
sqrtS = cms.double(sqrtS),
mTau    = cms.double(mTau),
mHiggs  = cms.double(mHiggs),
mZ      = cms.double(mZ),
InputType = cms.string(InputType),
InputFileList = cms.vstring(InputFileList),
maxNbrOfEventsToRead = cms.int32(maxNbrOfEventsToRead) ,
flagSameRNG = cms.bool(flagSameRNG) ,
flagTFLepTau = cms.bool(flagTFLepTau) , 
flagTFHadTau = cms.bool(flagTFHadTau) ,
flagTFFake   = cms.bool(flagTFFake) ,
flagTFMET    = cms.bool(flagTFMET) ,
flagTFJet1   = cms.bool(flagTFJet1) ,
flagTFJet2   = cms.bool(flagTFJet2) ,
flagTFBJet_leptop   = cms.bool(flagTFBJet_leptop) ,
flagTFBJet_hadtop   = cms.bool(flagTFBJet_hadtop) ,
flagTF_fakelep      = cms.bool(flagTF_fakelep) ,
flagTF_fakeleptau   = cms.bool(flagTF_fakeleptau) ,
flagTFTop    = cms.bool(flagTFTop) ,
flagJac      = cms.bool(flagJac) ,
flagWME      = cms.bool(flagWME) ,
verbose = cms.int32(verbose),
MEversion = cms.int32(MEversion),
FileOfIntegrals = cms.string(FileOfIntegrals),
#FileOfFunctionValues = cms.string(FileOfFunctionValues),

nbrOfDimttH = cms.int32(nbrOfDimttH),
nbrOfDimttZ = cms.int32(nbrOfDimttZ),
nbrOfDimttW = cms.int32(nbrOfDimttW),
nbrOfDimttjets = cms.int32(nbrOfDimttjets),
nbrOfDimttbar_SL = cms.int32(nbrOfDimttbar_SL),
nbrOfDimttbar_DL = cms.int32(nbrOfDimttbar_DL),
nbrOfDimttZ_Zll = cms.int32(nbrOfDimttZ_Zll),
nbrOfDimttH_miss = cms.int32(nbrOfDimttH_miss),
nbrOfDimttZ_miss = cms.int32(nbrOfDimttZ_miss),
nbrOfDimttjets_miss = cms.int32(nbrOfDimttjets_miss),
nbrOfDimttbar_SL_miss = cms.int32(nbrOfDimttbar_SL_miss),
nbrOfDimttZ_Zll_miss = cms.int32(nbrOfDimttZ_Zll_miss),

runttZ_integration = cms.bool(runttZ_integration),
runttW_integration = cms.bool(runttW_integration),
runttjets_integration = cms.bool(runttjets_integration),
runttbar_SL_integration = cms.bool(runttbar_SL_integration),
runttbar_DL_integration = cms.bool(runttbar_DL_integration),
runttbar_DL_fakelep_integration = cms.bool(runttbar_DL_fakelep_integration),
runttZ_Zll_integration = cms.bool(runttZ_Zll_integration),

run_missing_jet_integration = cms.bool(run_missing_jet_integration),
force_missing_jet_integration = cms.bool(force_missing_jet_integration),
force_missing_jet_integration_ifnoperm = cms.bool(force_missing_jet_integration_ifnoperm),

nbrOfPoints_ttZ = cms.int32(nbrOfPoints_ttZ),
nbrOfPoints_ttH = cms.int32(nbrOfPoints_ttH),
nbrOfPoints_ttW = cms.int32(nbrOfPoints_ttW),
nbrOfPoints_ttjets = cms.int32(nbrOfPoints_ttjets),
nbrOfPoints_ttbar_SL = cms.int32(nbrOfPoints_ttbar_SL),
nbrOfPoints_ttbar_DL = cms.int32(nbrOfPoints_ttbar_DL),
nbrOfPoints_ttZ_Zll = cms.int32(nbrOfPoints_ttZ_Zll),

nbrOfPoints_ttZ_miss = cms.int32(nbrOfPoints_ttZ_miss),
nbrOfPoints_ttH_miss = cms.int32(nbrOfPoints_ttH_miss),
nbrOfPoints_ttjets_miss = cms.int32(nbrOfPoints_ttjets_miss),
nbrOfPoints_ttbar_SL_miss = cms.int32(nbrOfPoints_ttbar_SL_miss),
nbrOfPoints_ttZ_Zll_miss = cms.int32(nbrOfPoints_ttZ_Zll_miss),

nbrOfPermut_per_jet = cms.int32(nbrOfPermut_per_jet),

phiNu_tlep = cms.vdouble(phiNu_tlep),
cosThetaNu_tlep = cms.vdouble(cosThetaNu_tlep),
phiNu_ttau = cms.vdouble(phiNu_ttau),
cosThetaNu_ttau = cms.vdouble(cosThetaNu_ttau),
phiNu_W_ttW = cms.vdouble(phiNu_W_ttW),
cosThetaNu_W_ttW = cms.vdouble(cosThetaNu_W_ttW),
 
CI_TFJet = cms.double(CI_TFJet),
use_pT_TFJet = cms.bool(use_pT_TFJet),
use_top_compatibility_check = cms.bool(use_top_compatibility_check),

include_hadrecoil = cms.bool(include_hadrecoil),
force_nonzero_integral = cms.bool(force_nonzero_integral),

eta_acceptance = cms.double(eta_acceptance),
jet_radius = cms.double(jet_radius),
dR_veto_jet_lep = cms.double(dR_veto_jet_lep),
rel_iso_lep = cms.double(rel_iso_lep),
pT_cut = cms.double(pT_cut),

phi_missing_jet = cms.vdouble(phi_missing_jet),
cosTheta_missing_jet = cms.vdouble(cosTheta_missing_jet),

treeName = cms.string(treeName),

MGParamCardFile = cms.string(MGParamCardFile),
LHAPDFFile = cms.string(LHAPDFFile),

#--# LUTOneProngPi0File = MEM_BASE + "/MGMEM/h_1P1PZeroSpectrumZ_LUT.lut"
#--# LUTThreeProngsFile = MEM_BASE + "/MGMEM/h_3PSpectrumZ_LUT.lut"
    
 )

