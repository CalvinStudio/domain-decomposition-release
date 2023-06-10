#ifndef JARVIS_DEVICE_STREAM_BONES
#define JARVIS_DEVICE_STREAM_BONES
#include "jarvis_1_host_field_meat.hpp"
#include "../.jarvis_1_vector/.jarvis_2_device_vector_out.h"
namespace jarvis
{
    struct jarvis_cuda_Stream
    {
    private:
        cudaStream_t _cal_stream;
        cudaStream_t _d2h_stream;
        cudaStream_t _h2d_stream;

    public:
        jarvis_cuda_Stream();
        void sync_cal_stream();
        void sync_d2h_stream();
        void sync_h2d_stream();
        cudaStream_t &cal_stream();
        cudaStream_t &d2h_stream();
        cudaStream_t &h2d_stream();
    };
}
#endif