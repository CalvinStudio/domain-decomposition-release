#ifndef _RTM3D_BASE_MODULE_MEAT
#define _RTM3D_BASE_MODULE_MEAT
#include "rtm3d_iamge_kernel.hpp"
namespace jarvis
{
    inline void MigImage::cuda_cal_correlation_from_wave()
    {
        void *args_list[] = {&fdtd_restruct_p->sub_model_p->gridext,
                             &fdtd_restruct_p->sub_model_p->gridphy,
                             //
                             &fdtd_restruct_p->sub_wave_p->vx.ptr(),
                             &fdtd_restruct_p->sub_wave_p->vy.ptr(),
                             &fdtd_restruct_p->sub_wave_p->vz.ptr(),
                             &fdtd_restruct_p->sub_wave_p->vxp.ptr(),
                             &fdtd_restruct_p->sub_wave_p->vyp.ptr(),
                             &fdtd_restruct_p->sub_wave_p->vzp.ptr(),
                             //
                             &fdtd_reverse_p->sub_wave_p->vx.ptr(),
                             &fdtd_reverse_p->sub_wave_p->vy.ptr(),
                             &fdtd_reverse_p->sub_wave_p->vz.ptr(),
                             &fdtd_reverse_p->sub_wave_p->vxp.ptr(),
                             &fdtd_reverse_p->sub_wave_p->vyp.ptr(),
                             &fdtd_reverse_p->sub_wave_p->vzp.ptr(),
                             //
                             &Ip_fd.ptr(),
                             &image_pp.ptr(),
                             &image_ps.ptr(),
                             &image_sp.ptr(),
                             &image_ss.ptr()};
        cudaLaunchKernel((void *)cal_correlation_from_wave, jarvis_cuda_kernel_size(fdtd_restruct_p->sub_model_p->gridphy.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
    }
    //
    inline void MigImage::cuda_cal_image_from_tmp_image()
    {
        void *args_list[] = {&fdtd_restruct_p->sub_model_p->gridphy,
                             &Ip_fd_tmp.ptr(),
                             &Ipp_tmp.ptr(),
                             &Ips_tmp.ptr(),
                             &Isp_tmp.ptr(),
                             &Iss_tmp.ptr(),
                             &Ip_fd.ptr(),
                             &image_pp.ptr(),
                             &image_ps.ptr(),
                             &image_sp.ptr(),
                             &image_ss.ptr()};
        cudaLaunchKernel((void *)cal_correlation_from_tmp_image, jarvis_cuda_kernel_size(fdtd_restruct_p->sub_model_p->gridphy.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
    }
    //
    inline void MigImage::cuda_cal_correlatione_tmp_image()
    {
        void *args_list[] = {&fdtd_restruct_p->sub_model_p->gridext,
                             &fdtd_restruct_p->sub_model_p->gridphy,
                             &fdtd_restruct_p->sub_model_p->margin,
                             &fdtd_reverse_p->sub_wave_p->vx.ptr(),
                             &fdtd_reverse_p->sub_wave_p->vy.ptr(),
                             &fdtd_reverse_p->sub_wave_p->vz.ptr(),
                             &fdtd_reverse_p->sub_wave_p->vxp.ptr(),
                             &fdtd_reverse_p->sub_wave_p->vyp.ptr(),
                             &fdtd_reverse_p->sub_wave_p->vzp.ptr(),
                             &fdtd_restruct_p->sub_wave_p->vx.ptr(),
                             &fdtd_restruct_p->sub_wave_p->vy.ptr(),
                             &fdtd_restruct_p->sub_wave_p->vz.ptr(),
                             &fdtd_restruct_p->sub_wave_p->vxp.ptr(),
                             &fdtd_restruct_p->sub_wave_p->vyp.ptr(),
                             &fdtd_restruct_p->sub_wave_p->vzp.ptr(),
                             &Ip_fd_tmp.ptr(),
                             &Ipp_tmp.ptr(),
                             &Ips_tmp.ptr(),
                             &Isp_tmp.ptr(),
                             &Iss_tmp.ptr()};
        cudaLaunchKernel((void *)cal_correlation_tmp, jarvis_cuda_kernel_size(fdtd_restruct_p->sub_model_p->gridphy.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
    }
    //
    inline void MigImage::cuda_source_compensation()
    {
        void *args_list[] = {&fdtd_restruct_p->sub_model_p->gridphy,
                             &Ip_fd.ptr(),
                             &image_pp.ptr(),
                             &image_ps.ptr(),
                             &image_sp.ptr(),
                             &image_ss.ptr()};
        cudaLaunchKernel((void *)source_compensation, jarvis_cuda_kernel_size(fdtd_restruct_p->sub_model_p->gridphy.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
    }
    inline void MigImage::cuda_laplace_filter()
    {
        void *args_list[] = {&fdtd_restruct_p->sub_model_p->gridphy,
                             &image_pp.ptr(),
                             &image_ps.ptr(),
                             &image_sp.ptr(),
                             &image_ss.ptr(),
                             &Ipp_laplace.ptr(),
                             &Ips_laplace.ptr(),
                             &Isp_laplace.ptr(),
                             &Iss_laplace.ptr()};
        cudaLaunchKernel((void *)laplace_filter_cal_for_image, jarvis_cuda_kernel_size(fdtd_restruct_p->sub_model_p->gridphy.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
    }
    // //
    // void ADCIGs::cuda_cal_poynting_vec_angle()
    // {
    //     cal_poynting_vec_angle<<<jarvis_cuda_kernel_size(fdtd_restruct_p->sub_model_p->gridphy.n_elem), 0, jarvis_cuda_stream_p->cal_stream>>>(
    //         fdtd_restruct_p->sub_model_p->gridphy, fdtd_restruct_p->pyt.pyt_vec.ptr(), fdtd_reverse_p->pyt.pyt_vec.ptr(), pyt_angle.ptr());
    // }
    // //
    // void ADCIGs::cuda_cal_adcig_image()
    // {
    //     for (int i_angle = 0; i_angle < adcig_angle_num; i_angle++)
    //     {
    //         cal_adcig_image_one_angle<<<jarvis_cuda_kernel_size(fdtd_restruct_p->sub_model_p->gridphy.n_elem), 0, jarvis_cuda_stream_p->cal_stream>>>(
    //             fdtd_restruct_p->sub_model_p->gridext, fdtd_restruct_p->sub_model_p->gridphy,
    //             fdtd_restruct_p->sub_model_p->margin,
    //             pyt_angle.ptr(), i_angle * max_adcig_angle / adcig_angle_num,
    //             image_p->Ipp_tmp.ptr(), image_p->Ips_tmp.ptr(),
    //             image_p->Isp_tmp.ptr(), image_p->Iss_tmp.ptr(),
    //             adcig_image(i_angle).image_pp.ptr(), adcig_image(i_angle).image_ps.ptr(),
    //             adcig_image(i_angle).image_sp.ptr(), adcig_image(i_angle).image_ss.ptr());
    //     }
    // }
}
#endif