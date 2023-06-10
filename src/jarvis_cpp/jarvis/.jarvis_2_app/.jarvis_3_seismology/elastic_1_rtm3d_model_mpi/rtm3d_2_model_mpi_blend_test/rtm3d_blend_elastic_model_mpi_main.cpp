#ifndef _RTM3D_BLEND_ELASTIC_MODEL_MPI_MAIN_CPP
#define _RTM3D_BLEND_ELASTIC_MODEL_MPI_MAIN_CPP
#include "../rtm3d_1_model_mpi_blend/.rtm3d_blend_elastic_model_mpi_header_out.h"
using namespace jarvis;
int main(int argc, char *argv[])
{
    jarvis_mpi_init(argc, argv, mpi_rank, mpi_size);
    if (mpi_rank == 0)
    {
        if (!argv[1])
        {
            std::cout << "Please input argfile!" << std::endl;
            std::abort();
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    //
    string proj_path = string(argv[1]);
    //
    jarvis_handle_cuda_error(cudaSetDevice(mpi_rank % GPU_NUM_PER_NODE));
    //
    for (int i = 0; i < mpi_size; i++)
    {
        if (mpi_rank == i)
        {
            jarvis_get_gpu_memory(mpi_rank % GPU_NUM_PER_NODE);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    jarvis_global_mpicu_Stream_t jarvis_mpi_cuda_stream;
    ArgList arg_list;
    rtm3d_blend_model_mpi_module rtm_mpi;
    arg_list.mpi_read_arg_file(proj_path, mpi_rank, mpi_size);
    jarvis_mpi_cuda_stream.initialize(arg_list.x_divide_num, arg_list.y_divide_num, arg_list.z_divide_num,
                                       mpi_rank, mpi_size, arg_list.slots_per_node);
    rtm_mpi.link(&arg_list, &jarvis_mpi_cuda_stream);
    device_mpi_sub_rtm3d(&rtm_mpi);
    jarvis_mpi_final;
    return 0;
}
#endif