#ifndef _RTM3D_BASE_MODULE_HPP
#define _RTM3D_BASE_MODULE_HPP
#include ".rtm3d_image_header_in.h"
namespace jarvis
{
    struct MigImage
    {
        jarvis_global_mpicu_Stream_t *jarvis_mpi_cuda_stream_p = nullptr;
        fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_backward> *fdtd_restruct_p = 0;
        fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_reverse> *fdtd_reverse_p = 0;
        //
        Field<float,MemType::device> Ipp_tmp;
        Field<float,MemType::device> Ips_tmp;
        Field<float,MemType::device> Isp_tmp;
        Field<float,MemType::device> Iss_tmp;
        Field<float,MemType::device> Ip_fd_tmp;
        //
        Field<float,MemType::device> image_pp;
        Field<float,MemType::device> image_ps;
        Field<float,MemType::device> image_sp;
        Field<float,MemType::device> image_ss;
        Field<float,MemType::device> Ip_fd;

        Field<float,MemType::device> Ipp_laplace;
        Field<float,MemType::device> Ips_laplace;
        Field<float,MemType::device> Isp_laplace;
        Field<float,MemType::device> Iss_laplace;

        Field<float,MemType::device> Ipp_laplace_multi_shot;
        Field<float,MemType::device> Ips_laplace_multi_shot;
        Field<float,MemType::device> Isp_laplace_multi_shot;
        Field<float,MemType::device> Iss_laplace_multi_shot;
        //
        void link(jarvis_global_mpicu_Stream_t *_jarvis_mpi_cuda_stream_p,
                         fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_backward> *_fdtd_restruct_p,
                         fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_reverse> *_fdtd_reverse_p)
        {
            jarvis_mpi_cuda_stream_p = _jarvis_mpi_cuda_stream_p;
            fdtd_restruct_p = _fdtd_restruct_p;
            fdtd_reverse_p = _fdtd_reverse_p;
        }
        //
        void initialize()
        {
            Frame &_gridphy = fdtd_restruct_p->sub_model_p->gridphy;
            // Ipp_tmp.alloc(_gridphy);
            // Ips_tmp.alloc(_gridphy);
            // Isp_tmp.alloc(_gridphy);
            // Iss_tmp.alloc(_gridphy);
            // Ip_fd_tmp.alloc(_gridphy);
            //
            image_pp.alloc(_gridphy);
            image_ps.alloc(_gridphy);
            image_sp.alloc(_gridphy);
            image_ss.alloc(_gridphy);
            Ip_fd.alloc(_gridphy);
            //
            Ipp_laplace.alloc(_gridphy);
            Ips_laplace.alloc(_gridphy);
            Isp_laplace.alloc(_gridphy);
            Iss_laplace.alloc(_gridphy);
        }
        //
        //
        void alloc_image_only_host()
        {
            Frame &_gridphy = fdtd_restruct_p->sub_model_p->gridphy;
            //
            image_pp.alloc(_gridphy);
            image_ps.alloc(_gridphy);
            image_sp.alloc(_gridphy);
            image_ss.alloc(_gridphy);
            Ip_fd.alloc(_gridphy);
            //
            Ipp_laplace.alloc(_gridphy);
            Ips_laplace.alloc(_gridphy);
            Isp_laplace.alloc(_gridphy);
            Iss_laplace.alloc(_gridphy);
        }
        //
        void alloc_tmp_image_only_device()
        {
            Frame &_gridphy = fdtd_restruct_p->sub_model_p->gridphy;
            Ipp_tmp.alloc(_gridphy);
            Ips_tmp.alloc(_gridphy);
            Isp_tmp.alloc(_gridphy);
            Iss_tmp.alloc(_gridphy);
            Ip_fd_tmp.alloc(_gridphy);
        }
        void cuda_cal_correlation_from_wave();
        void cuda_cal_image_from_tmp_image();
        void cuda_cal_correlatione_tmp_image();
        void cuda_source_compensation();
        void cuda_laplace_filter();
    };
    //
    // struct ADCIGs
    // {
    //     //
    //     fdtd3d_blend_elastic_model_mpi_module *fdtd_restruct_p = 0;
    //     fdtd3d_blend_elastic_model_mpi_module *fdtd_reverse_p = 0;
    //     MigImage *image_p = 0;
    //     jarvis_cuda_Stream *jarvis_cuda_stream_p = 0;
    //     //
    //     int adcig_angle_num;
    //     float max_adcig_angle;
    //     Field<float,MemType::device> pyt_angle;
    //     MigImage adcig_image_tmp;
    //     Field<MigImage> adcig_image;
    //     //
    //     void link_mpi_cuda_stream(jarvis_cuda_Stream &_jarvis_cuda_stream)
    //     {
    //         jarvis_cuda_stream_p = &_jarvis_cuda_stream;
    //         adcig_image_tmp.link_mpi_cuda_stream(_jarvis_cuda_stream);
    //     }
    //     //
    //     void link(fdtd3d_blend_elastic_model_mpi_module *_fdtd_restruct_p,
    //                      fdtd3d_blend_elastic_model_mpi_module *_fdtd_reverse_p,
    //                      MigImage *_image_p)
    //     {
    //         fdtd_restruct_p = _fdtd_restruct_p;
    //         fdtd_reverse_p = _fdtd_reverse_p;
    //         image_p = _image_p;
    //         adcig_image_tmp.link(fdtd_restruct_p, fdtd_reverse_p);
    //     }
    //     //
    //     void init(float _max_adcig_angle, int _adcig_angle_num)
    //     {
    //         //
    //         pyt_angle.alloc( fdtd_restruct_p->sub_model_p->gridphy);
    //         max_adcig_angle = _max_adcig_angle;
    //         adcig_angle_num = _adcig_angle_num;
    //         adcig_image.alloc(_adcig_angle_num);
    //         adcig_image_tmp.alloc_tmp_image_only_device();
    //         for (int i = 0; i < _adcig_angle_num; i++)
    //         {
    //             adcig_image(i).link_mpi_cuda_stream(*jarvis_cuda_stream_p);
    //             adcig_image(i).link(fdtd_restruct_p, fdtd_reverse_p);
    //             adcig_image(i).alloc_image_only_host();
    //         }
    //     }
    //     //
    //     void cuda_cal_poynting_vec_angle();
    //     void cuda_cal_adcig_image();
    // };
}
#endif