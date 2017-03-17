/* 
 * File:   KernelEmu.h
 * Author: grasseau
 *
 * Created on 7 juillet 2014, 17:56
 */

/* 
 * Need to use kernels outside OpenCL (for tests)
 */

#ifndef MEM_MEMAlgo_KERNELEMU_H
#define	MEM_MEMAlgo_KERNELEMU_H

# include <stddef.h>

#ifdef	__cplusplus
extern "C" {
#endif
# define CLK_LOCAL_MEM_FENCE 1
# define CLK_GLOBAL_MEM_FENCE 2

void KernelEmu_initKernel();

size_t get_global_id  ( unsigned int dim);
size_t get_global_size( unsigned int dim); 
size_t get_local_id  ( unsigned int dim);
size_t get_local_size ( unsigned int dim);
size_t get_group_id( unsigned int dim);
size_t get_num_groups( unsigned int dim);
void   barrier        ( unsigned int type);
int    atomic_xchg    ( int *p, int val);
int    atomic_inc    ( int *p);
int    atomic_dec    ( int *p);
#ifdef	__cplusplus
}
#endif

#endif	/* KERNELEMU_H */

