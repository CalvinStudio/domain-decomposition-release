#pragma once
#ifndef _FDTD3D_BLEND_ELASTIC_MODEL_MPI_MAIN
#define _FDTD3D_BLEND_ELASTIC_MODEL_MPI_MAIN
#include "fdtd3d_blend_elastic_model_mpi_2_device_meat.hpp"
using namespace jarvis;
int main(int argc, char *argv[])
{
    if (!argv[1])
    {
        std::cout << "Please input argfile!" << std::endl;
        std::abort();
    }
    //
    jarvis_mpi_init(argc, argv, mpi_rank, mpi_size);
    printf("Thread %d of %d\n", mpi_rank, mpi_size);
    //
    string arg_path = string(argv[1]) + "/model/argfile.json";
    jarvis_handle_cuda_error(cudaSetDevice(mpi_rank % GPU_NUM_PER_NODE)); //* Select the GPU device to start.
    //*The following codes operate on this GPU. Be sure to define CUDA objects later
    //
    cudaStream_t cal_stream;
    cudaStream_t data_stream;
    //
    cudaStreamCreate(&cal_stream);
    cudaStreamCreate(&data_stream);
    fdtd3d_blend_elastic_model_mpi_module fdtd_model_mpi;
    fdtd_model_mpi.link(cal_stream, data_stream);
    fdtd_model_mpi.mpi_read_arg_file(arg_path, mpi_rank, mpi_size);
    fdtd_model_mpi.mpi_sub_fdtd3d();
    jarvis_mpi_final;
    return 0;
}
#endif
