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
