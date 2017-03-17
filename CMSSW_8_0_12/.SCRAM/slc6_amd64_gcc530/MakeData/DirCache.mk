ifeq ($(strip $(MEM/MEMAlgo)),)
ALL_COMMONRULES += src_MEM_MEMAlgo_src
src_MEM_MEMAlgo_src_parent := MEM/MEMAlgo
src_MEM_MEMAlgo_src_INIT_FUNC := $$(eval $$(call CommonProductRules,src_MEM_MEMAlgo_src,src/MEM/MEMAlgo/src,LIBRARY))
MEMMEMAlgo := self/MEM/MEMAlgo
MEM/MEMAlgo := MEMMEMAlgo
MEMMEMAlgo_files := $(patsubst src/MEM/MEMAlgo/src/%,%,$(wildcard $(foreach dir,src/MEM/MEMAlgo/src ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
MEMMEMAlgo_BuildFile    := $(WORKINGDIR)/cache/bf/src/MEM/MEMAlgo/BuildFile
MEMMEMAlgo_LOC_FLAGS_CXXFLAGS   :=  -g -I /opt/exp_soft/vo.gridcl.fr/software/OpenMPI/slc6_amd64_gcc491/openmpi/1.6.5/include -I /opt/exp_soft/vo.gridcl.fr/software/OpenCL/Intel-SDK/include -pthread 
MEMMEMAlgo_LOC_FLAGS_LDFLAGS   :=  -g -v -I /opt/exp_soft/vo.gridcl.fr/software/OpenMPI/slc6_amd64_gcc491/openmpi/1.6.5/include -L/opt/exp_soft/vo.gridcl.fr/software/OpenMPI/slc6_amd64_gcc491/openmpi/1.6.5/lib -pthread 
MEMMEMAlgo_LOC_USE := self  root python lhapdf FWCore/Framework FWCore/PluginManager FWCore/ParameterSet DataFormats/Math MEMDataFormats/MEMDataFormats
MEMMEMAlgo_EX_LIB   := MEMMEMAlgo
MEMMEMAlgo_EX_USE   := $(foreach d,$(MEMMEMAlgo_LOC_USE),$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
MEMMEMAlgo_PACKAGE := self/src/MEM/MEMAlgo/src
ALL_PRODS += MEMMEMAlgo
MEMMEMAlgo_CLASS := LIBRARY
MEM/MEMAlgo_forbigobj+=MEMMEMAlgo
MEMMEMAlgo_INIT_FUNC        += $$(eval $$(call Library,MEMMEMAlgo,src/MEM/MEMAlgo/src,src_MEM_MEMAlgo_src,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_LIB),$(SCRAMSTORENAME_LOGS)))
endif
ifeq ($(strip $(1)),)
1 := self/src/MEM/MEMProducer/plugins
PLUGINS:=yes
1_files := $(patsubst src/MEM/MEMProducer/plugins/%,%,$(foreach file,*.cc,$(eval xfile:=$(wildcard src/MEM/MEMProducer/plugins/$(file)))$(if $(xfile),$(xfile),$(warning No such file exists: src/MEM/MEMProducer/plugins/$(file). Please fix src/MEM/MEMProducer/plugins/BuildFile.))))
1_BuildFile    := $(WORKINGDIR)/cache/bf/src/MEM/MEMProducer/plugins/BuildFile
1_LOC_FLAGS_CXXFLAGS   :=  -g -I /opt/exp_soft/vo.gridcl.fr/software/OpenMPI/slc6_amd64_gcc491/openmpi/1.6.5/include -I /opt/exp_soft/vo.gridcl.fr/software/OpenCL/Intel-SDK/include -pthread 
1_LOC_FLAGS_LDFLAGS   :=  -g -v -I /opt/exp_soft/vo.gridcl.fr/software/OpenMPI/slc6_amd64_gcc491/openmpi/1.6.5/include -L/opt/exp_soft/vo.gridcl.fr/software/OpenMPI/slc6_amd64_gcc491/openmpi/1.6.5/lib -pthread 
1_LOC_USE := self  FWCore/Framework FWCore/PluginManager FWCore/ParameterSet DataFormats/TrackReco DataFormats/Math DataFormats/HepMCCandidate DataFormats/PatCandidates MEM/MEMAlgo MEMDataFormats/MEMDataFormats
1_PRE_INIT_FUNC += $$(eval $$(call edmPlugin,1,1,$(SCRAMSTORENAME_LIB),src/MEM/MEMProducer/plugins))
1_PACKAGE := self/src/MEM/MEMProducer/plugins
ALL_PRODS += 1
MEM/MEMProducer_forbigobj+=1
1_INIT_FUNC        += $$(eval $$(call Library,1,src/MEM/MEMProducer/plugins,src_MEM_MEMProducer_plugins,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_LIB),$(SCRAMSTORENAME_LOGS)))
1_CLASS := LIBRARY
else
$(eval $(call MultipleWarningMsg,1,src/MEM/MEMProducer/plugins))
endif
ALL_COMMONRULES += src_MEM_MEMProducer_plugins
src_MEM_MEMProducer_plugins_parent := MEM/MEMProducer
src_MEM_MEMProducer_plugins_INIT_FUNC += $$(eval $$(call CommonProductRules,src_MEM_MEMProducer_plugins,src/MEM/MEMProducer/plugins,PLUGINS))
ifeq ($(strip $(MEMDataFormats/MEMDataFormats)),)
ALL_COMMONRULES += src_MEMDataFormats_MEMDataFormats_src
src_MEMDataFormats_MEMDataFormats_src_parent := MEMDataFormats/MEMDataFormats
src_MEMDataFormats_MEMDataFormats_src_INIT_FUNC := $$(eval $$(call CommonProductRules,src_MEMDataFormats_MEMDataFormats_src,src/MEMDataFormats/MEMDataFormats/src,LIBRARY))
MEMDataFormatsMEMDataFormats := self/MEMDataFormats/MEMDataFormats
MEMDataFormats/MEMDataFormats := MEMDataFormatsMEMDataFormats
MEMDataFormatsMEMDataFormats_files := $(patsubst src/MEMDataFormats/MEMDataFormats/src/%,%,$(wildcard $(foreach dir,src/MEMDataFormats/MEMDataFormats/src ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
MEMDataFormatsMEMDataFormats_BuildFile    := $(WORKINGDIR)/cache/bf/src/MEMDataFormats/MEMDataFormats/BuildFile
MEMDataFormatsMEMDataFormats_LOC_USE := self  DataFormats/Candidate DataFormats/Common DataFormats/EgammaReco DataFormats/CaloRecHit DataFormats/CaloTowers DataFormats/Math DataFormats/RecoCandidate DataFormats/TrackReco DataFormats/TrackerRecHit2D DataFormats/TrackingRecHit DataFormats/VertexReco DataFormats/GeometryCommonDetAlgo FWCore/MessageLogger FWCore/Framework FWCore/PluginManager FWCore/ParameterSet DataFormats/HepMCCandidate DataFormats/PatCandidates
MEMDataFormatsMEMDataFormats_LCGDICTS  := x 
MEMDataFormatsMEMDataFormats_PRE_INIT_FUNC += $$(eval $$(call LCGDict,MEMDataFormatsMEMDataFormats,src/MEMDataFormats/MEMDataFormats/src/classes.h,src/MEMDataFormats/MEMDataFormats/src/classes_def.xml,$(SCRAMSTORENAME_LIB),$(GENREFLEX_ARGS) --fail_on_warnings,))
MEMDataFormatsMEMDataFormats_EX_LIB   := MEMDataFormatsMEMDataFormats
MEMDataFormatsMEMDataFormats_EX_USE   := $(foreach d,$(MEMDataFormatsMEMDataFormats_LOC_USE),$(if $($(d)_EX_FLAGS_NO_RECURSIVE_EXPORT),,$d))
MEMDataFormatsMEMDataFormats_PACKAGE := self/src/MEMDataFormats/MEMDataFormats/src
ALL_PRODS += MEMDataFormatsMEMDataFormats
MEMDataFormatsMEMDataFormats_CLASS := LIBRARY
MEMDataFormats/MEMDataFormats_forbigobj+=MEMDataFormatsMEMDataFormats
MEMDataFormatsMEMDataFormats_INIT_FUNC        += $$(eval $$(call Library,MEMDataFormatsMEMDataFormats,src/MEMDataFormats/MEMDataFormats/src,src_MEMDataFormats_MEMDataFormats_src,$(SCRAMSTORENAME_BIN),,$(SCRAMSTORENAME_LIB),$(SCRAMSTORENAME_LOGS)))
endif
