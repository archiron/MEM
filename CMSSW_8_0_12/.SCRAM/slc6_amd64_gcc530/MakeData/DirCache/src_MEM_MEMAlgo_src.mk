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
