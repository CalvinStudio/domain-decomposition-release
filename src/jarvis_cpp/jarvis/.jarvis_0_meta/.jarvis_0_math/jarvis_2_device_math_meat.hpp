#ifndef _JARVIS_DEVICE_MATH_MEAT
#define _JARVIS_DEVICE_MATH_MEAT
//**********************************Developer******************************************
// 2020.04.10 BY CAIWEI CALVIN CAI
//*************************************************************************************
#include "jarvis_2_device_math_bones.hpp"
namespace jarvis
{
	inline void handle_cuda_error(cudaError_t err, const char *file, int line)
	{
		if (err != cudaSuccess)
		{
			printf("\033[41;37m[error]: %s in %s at line %d\033[0m\n", cudaGetErrorString(err), file, line);
			std::abort();
		}
	}
	//
	inline void jarvis_get_gpu_memory(int i_device)
	{
		int deviceCount = 0;
		cudaError_t error_id = cudaGetDeviceCount(&deviceCount);
		if (deviceCount == 0)
		{
			std::cout << "\033[41;37mThe current PC does not support CUDA devices!\033[0m" << std::endl;
		}

		size_t gpu_total_size;
		size_t gpu_free_size;

		cudaError_t cuda_status = cudaMemGetInfo(&gpu_free_size, &gpu_total_size);

		if (cudaSuccess != cuda_status)
		{
			std::cout << "Error: cudaMemGetInfo fails : " << cudaGetErrorString(cuda_status) << std::endl;
			std::abort();
		}

		double total_memory = double(gpu_total_size) / (1024.0 * 1024.0);
		double free_memory = double(gpu_free_size) / (1024.0 * 1024.0);
		double used_memory = total_memory - free_memory;

		std::string free_memory_str = "[free\tcuda memory]:\t";
		if (free_memory < 10 * 1024.0)
		{
			std::cout << "[device\tnumber]:\t" << i_device << "\n"
					  << "[total\tcuda memory]:\t" << total_memory << " MB\n"
					  << "[used\tcuda memory]:\t" << used_memory << " MB\n"
					  << "\033[41;37m[free\tcuda memory]:\t" << free_memory << " MB\033[0m\n"
					  << std::endl;
		}
		else
		{
			std::cout << "[device\tnumber]:\t" << i_device << "\n"
					  << "[total\tcuda memory]:\t" << total_memory << " MB\n"
					  << "[used\tcuda memory]:\t" << used_memory << " MB\n"
					  << "[free\tcuda memory]:\t" << free_memory << " MB\n"
					  << std::endl;
		}

		if (free_memory < 10 * 1024.0)
		{
			printf("\033[41;37m[error]: free memory is not enough!\033[0m\n");
			std::abort();
		}
	}
	//
	inline void print_gpu_info()
	{
		cudaDeviceProp deviceProp;
		int deviceCount;
		cudaGetDeviceCount(&deviceCount);
		for (int i = 0; i < deviceCount; ++i)
		{
			cudaGetDeviceProperties(&deviceProp, i);
			std::cout << "Device " << i + 1 << ":" << std::endl;
			std::cout << "Graphics card model:" << deviceProp.name << std::endl;
			std::cout << "Total device global memory in MB:"
					  << deviceProp.totalGlobalMem / 1024 / 1024 << std::endl;
			std::cout << "The maximum available shared memory (in KB) in a thread block on the device:"
					  << deviceProp.sharedMemPerBlock / 1024 << std::endl;
			std::cout << "Number of 32-bit registers available for a thread block on the device:"
					  << deviceProp.regsPerBlock << std::endl;
			std::cout << "The maximum number of threads that a thread block on a device can contain:"
					  << deviceProp.maxThreadsPerBlock << std::endl;
			std::cout << "Version number of the device's compute capability:"
					  << deviceProp.major << "." << deviceProp.minor << std::endl;
			std::cout << "Number of multiprocessors on the device:" << deviceProp.multiProcessorCount << std::endl;
			std::cout << "canMapHostMemory:" << deviceProp.canMapHostMemory << std::endl;
			int canAccessPeer;

			cudaDeviceCanAccessPeer(&canAccessPeer, 0, 4);
			std::cout << "support peer to peer==%d" << canAccessPeer << std::endl;
		}
	}

	inline __device__ bool cuJacobi(float *matrix, int dim, float *eigenvectors, float *eigenvalues, float precision, int max)
	{
		for (int i = 0; i < dim; ++i)
		{
			eigenvectors[i * dim + i] = 1.0f;
			for (int j = 0; j < dim; j++)
			{
				if (i != j)
					eigenvectors[i * dim + j] = 0.0f;
			}
		}

		int nCount = 0; // current iteration
		while (1)
		{
			// find the largest element on the off-diagonal line of the matrix
			double dbMax = matrix[1];
			int nRow = 0;
			int nCol = 1;
			for (int i = 0; i < dim; ++i)
			{ // row
				for (int j = 0; j < dim; j++)
				{ // column
					double d = fabs(matrix[i * dim + j]);
					if ((i != j) && (d > dbMax))
					{
						dbMax = d;
						nRow = i;
						nCol = j;
					}
				}
			}

			if (dbMax < precision) // precision check
				break;
			if (nCount > max) // iterations check
				break;
			nCount++;

			double dbApp = matrix[nRow * dim + nRow];
			double dbApq = matrix[nRow * dim + nCol];
			double dbAqq = matrix[nCol * dim + nCol];
			// compute rotate angle
			double dbAngle = 0.5 * atan2(-2 * dbApq, dbAqq - dbApp);
			double dbSinTheta = sin(dbAngle);
			double dbCosTheta = cos(dbAngle);
			double dbSin2Theta = sin(2 * dbAngle);
			double dbCos2Theta = cos(2 * dbAngle);
			matrix[nRow * dim + nRow] = dbApp * dbCosTheta * dbCosTheta +
										dbAqq * dbSinTheta * dbSinTheta + 2 * dbApq * dbCosTheta * dbSinTheta;
			matrix[nCol * dim + nCol] = dbApp * dbSinTheta * dbSinTheta +
										dbAqq * dbCosTheta * dbCosTheta - 2 * dbApq * dbCosTheta * dbSinTheta;
			matrix[nRow * dim + nCol] = 0.5 * (dbAqq - dbApp) * dbSin2Theta + dbApq * dbCos2Theta;
			matrix[nCol * dim + nRow] = matrix[nRow * dim + nCol];

			for (int i = 0; i < dim; ++i)
			{
				if ((i != nCol) && (i != nRow))
				{
					int u = i * dim + nRow; // p
					int w = i * dim + nCol; // q
					dbMax = matrix[u];
					matrix[u] = matrix[w] * dbSinTheta + dbMax * dbCosTheta;
					matrix[w] = matrix[w] * dbCosTheta - dbMax * dbSinTheta;
				}
			}

			for (int j = 0; j < dim; j++)
			{
				if ((j != nCol) && (j != nRow))
				{
					int u = nRow * dim + j; // p
					int w = nCol * dim + j; // q
					dbMax = matrix[u];
					matrix[u] = matrix[w] * dbSinTheta + dbMax * dbCosTheta;
					matrix[w] = matrix[w] * dbCosTheta - dbMax * dbSinTheta;
				}
			}

			// compute eigenvector
			for (int i = 0; i < dim; ++i)
			{
				int u = i * dim + nRow; // p
				int w = i * dim + nCol; // q
				dbMax = eigenvectors[u];
				eigenvectors[u] = eigenvectors[w] * dbSinTheta + dbMax * dbCosTheta;
				eigenvectors[w] = eigenvectors[w] * dbCosTheta - dbMax * dbSinTheta;
			}
		}

		for (int i = 0; i < dim; ++i)
		{
			eigenvalues[i] = matrix[i * dim + i];
		}
		return true;
	}

	inline __device__ bool CudaNextPermutation(int *p, int n)
	{
		int last = n - 1;
		int i, j, k;
		i = last;
		while (i > 0 && p[i] < p[i - 1])
			i--;
		if (i == 0)
			return false;
		k = i;
		for (j = last; j >= i; j--)
			if (p[j] > p[i - 1] && p[j] < p[k])
				k = j;
		Swap(int, p[k], p[i - 1]);
		for (j = last, k = i; j > k; j--, k++)
			Swap(int, p[j], p[k]);
		return true;
	}
}
#endif