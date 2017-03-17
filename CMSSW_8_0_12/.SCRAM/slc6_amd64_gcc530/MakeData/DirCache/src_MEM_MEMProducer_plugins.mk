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
