ALL_COMMONRULES += src_MEM_MEMAlgo_test
src_MEM_MEMAlgo_test_parent := MEM/MEMAlgo
src_MEM_MEMAlgo_test_INIT_FUNC += $$(eval $$(call CommonProductRules,src_MEM_MEMAlgo_test,src/MEM/MEMAlgo/test,TEST))
