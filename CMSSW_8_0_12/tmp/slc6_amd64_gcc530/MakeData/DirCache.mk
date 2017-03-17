ALL_SUBSYSTEMS+=MEM
subdirs_src_MEM = src_MEM_MEMAlgo src_MEM_MEMProducer
ALL_PACKAGES += MEM/MEMAlgo
subdirs_src_MEM_MEMAlgo := src_MEM_MEMAlgo_test src_MEM_MEMAlgo_src
ALL_COMMONRULES += src_MEM_MEMAlgo_test
src_MEM_MEMAlgo_test_parent := MEM/MEMAlgo
src_MEM_MEMAlgo_test_INIT_FUNC += $$(eval $$(call CommonProductRules,src_MEM_MEMAlgo_test,src/MEM/MEMAlgo/test,TEST))
ALL_PACKAGES += MEM/MEMProducer
subdirs_src_MEM_MEMProducer := src_MEM_MEMProducer_plugins src_MEM_MEMProducer_python src_MEM_MEMProducer_test src_MEM_MEMProducer_doc
ifeq ($(strip $(PyMEMMEMProducer)),)
PyMEMMEMProducer := self/src/MEM/MEMProducer/python
src_MEM_MEMProducer_python_parent := 
ALL_PYTHON_DIRS += $(patsubst src/%,%,src/MEM/MEMProducer/python)
PyMEMMEMProducer_files := $(patsubst src/MEM/MEMProducer/python/%,%,$(wildcard $(foreach dir,src/MEM/MEMProducer/python ,$(foreach ext,$(SRC_FILES_SUFFIXES),$(dir)/*.$(ext)))))
PyMEMMEMProducer_LOC_USE := self  
PyMEMMEMProducer_PACKAGE := self/src/MEM/MEMProducer/python
ALL_PRODS += PyMEMMEMProducer
PyMEMMEMProducer_INIT_FUNC        += $$(eval $$(call PythonProduct,PyMEMMEMProducer,src/MEM/MEMProducer/python,src_MEM_MEMProducer_python,1,1,$(SCRAMSTORENAME_PYTHON),$(SCRAMSTORENAME_LIB),,))
else
$(eval $(call MultipleWarningMsg,PyMEMMEMProducer,src/MEM/MEMProducer/python))
endif
ALL_COMMONRULES += src_MEM_MEMProducer_python
src_MEM_MEMProducer_python_INIT_FUNC += $$(eval $$(call CommonProductRules,src_MEM_MEMProducer_python,src/MEM/MEMProducer/python,PYTHON))
ALL_COMMONRULES += src_MEM_MEMProducer_test
src_MEM_MEMProducer_test_parent := MEM/MEMProducer
src_MEM_MEMProducer_test_INIT_FUNC += $$(eval $$(call CommonProductRules,src_MEM_MEMProducer_test,src/MEM/MEMProducer/test,TEST))
ALL_SUBSYSTEMS+=MEMDataFormats
subdirs_src_MEMDataFormats = src_MEMDataFormats_MEMDataFormats
ALL_PACKAGES += MEMDataFormats/MEMDataFormats
subdirs_src_MEMDataFormats_MEMDataFormats := src_MEMDataFormats_MEMDataFormats_doc src_MEMDataFormats_MEMDataFormats_src
