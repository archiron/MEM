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
