ALL_COMMONRULES += src_MEM_MEMProducer_test
src_MEM_MEMProducer_test_parent := MEM/MEMProducer
src_MEM_MEMProducer_test_INIT_FUNC += $$(eval $$(call CommonProductRules,src_MEM_MEMProducer_test,src/MEM/MEMProducer/test,TEST))
