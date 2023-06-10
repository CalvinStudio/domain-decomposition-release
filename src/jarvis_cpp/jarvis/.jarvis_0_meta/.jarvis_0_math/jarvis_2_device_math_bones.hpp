#ifndef _JARVIS_DEVICE_MATH_BONES
#define _JARVIS_DEVICE_MATH_BONES
//**********************************Developer******************************************
// 2020.04.10 BY CAIWEI CALVIN CAI
//*************************************************************************************
#include "../../../.external_header/.public_cuda_header.h"
#include "jarvis_1_host_math_meat.hpp"
//************************************ARRAY********************************************
//
#define jarvis_default_cuda_stream (cudaStream_t)0
//
#define jarvis_handle_cuda_error(err) (handle_cuda_error(err, __FILE__, __LINE__))
//
#define jarvis_cuda_kernel_size(_cuda_threads_num) (_cuda_threads_num + jarvis_const::block_size - 1) / jarvis_const::block_size, jarvis_const::block_size
//
#define GPU_NUM_PER_NODE 4
//
#define Swap(type, a, b) \
	{                    \
		type temp;       \
		temp = a;        \
		a = b;           \
		b = temp;        \
	}
//
#define cuda_tic(mark)                            \
	cudaEvent_t start_time##mark, end_time##mark; \
	cudaEventCreate(&start_time##mark);           \
	cudaEventCreate(&end_time##mark);             \
	cudaEventRecord(start_time##mark, 0)
//
#define cuda_toc(mark, str)                                              \
	cudaEventRecord(end_time##mark, 0);                                  \
	cudaEventSynchronize(start_time##mark);                              \
	cudaEventSynchronize(end_time##mark);                                \
	float time##mark;                                                    \
	cudaEventElapsedTime(&time##mark, start_time##mark, end_time##mark); \
	std::cout << str << ": " << std::fixed << std::setprecision(4) << time##mark << " ms" << std::endl
//
namespace jarvis
{
	void handle_cuda_error(cudaError_t err, const char *file, int line);
	void jarvis_get_gpu_memory(int i_device);
	void print_gpu_info();
	__device__ bool cuJacobi(float *matrix, int dim, float *eigenvectors, float *eigenvalues, float precision, int max);
	__device__ bool CudaNextPermutation(int *p, int n);
}
#endif