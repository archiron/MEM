# include "MEM/MEMAlgo/interface/KernelEmu.h"

static size_t global_id;

void KernelEmu_initKernel() { global_id=0; };

size_t get_global_id( unsigned int dim) {
  /*
  size_t id = global_id; 
  // global_id++; 
  return id;
  */
  return 0;
}

size_t get_global_size( unsigned int dim) {
  // ??? size_t id = global_id; global_id++; return id;
}

size_t get_local_id( unsigned int dim) {
  // ???size_t id = global_id; global_id++; return id;
  return 0;
}

size_t get_local_size( unsigned int dim) {
  // ???size_t id = global_id; global_id++; return id;
  return 1;
}

size_t get_group_id( unsigned int dim) {
  return 0;
}

size_t get_num_groups( unsigned int dim) {
  return 1;
}

void barrier( unsigned int type) {
}

int   atomic_xchg ( int *p, int val) {}
int    atomic_inc    ( int *p) {}
int    atomic_dec    ( int *p) {}
