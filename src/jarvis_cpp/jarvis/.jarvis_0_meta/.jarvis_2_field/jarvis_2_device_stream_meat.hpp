#ifndef JARVIS_DEVICE_STREAM_MEAT
#define JARVIS_DEVICE_STREAM_MEAT
#include "jarvis_2_device_stream_bones.hpp"
namespace jarvis
{
    inline jarvis_cuda_Stream::jarvis_cuda_Stream()
    {
        cudaStreamCreate(&_cal_stream);
        cudaStreamCreate(&_d2h_stream);
        cudaStreamCreate(&_h2d_stream);
    }
    inline void jarvis_cuda_Stream::sync_cal_stream()
    {
        cudaStreamSynchronize(_cal_stream);
    }
    inline void jarvis_cuda_Stream::sync_d2h_stream()
    {
        cudaStreamSynchronize(_d2h_stream);
    }
    inline void jarvis_cuda_Stream::sync_h2d_stream()
    {
        cudaStreamSynchronize(_h2d_stream);
    }
    inline cudaStream_t &jarvis_cuda_Stream::cal_stream()
    {
        return _cal_stream;
    }
    inline cudaStream_t &jarvis_cuda_Stream::d2h_stream()
    {
        return _d2h_stream;
    }
    inline cudaStream_t &jarvis_cuda_Stream::h2d_stream()
    {
        return _h2d_stream;
    }
}
#endif