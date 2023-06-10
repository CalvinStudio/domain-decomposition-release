#ifndef _RTM3D_BLEND_ELASTIC_MODEL_MPI_DEVICE
#define _RTM3D_BLEND_ELASTIC_MODEL_MPI_DEVICE
#include "../rtm3d_1_model_mpi_blend/.rtm3d_blend_elastic_model_mpi_header_out.h"
namespace jarvis
{
    jarvis_global_mpicu_Stream_t *elastic_fdtd3d_module_Base::jarvis_mpi_cuda_stream_p = nullptr;
    ArgList *elastic_fdtd3d_module_Base::arg_list_p = nullptr;
    void device_mpi_sub_rtm3d(rtm3d_blend_model_mpi_module *rtm_mpi)
    {
        rtm_mpi->mpi_sub_rtm3d();
    }
}
#endif