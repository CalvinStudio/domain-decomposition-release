#pragma once
#ifndef _FDTD3D_PHY_MPI_MEAT_HOST_MEAT_HPP
#define _FDTD3D_PHY_MPI_MEAT_HOST_MEAT_HPP
#include "fdtd3d_phy_mpi_0.hpp"
namespace jarvis
{
    inline void mpicuPhy_Base::link(elasticGeoModel<domain::local> *_geomodel_p,
                                    SeisRecord *_seisrecord_p,
                                    elasticWave *_sub_ext_wavefield_p,
                                    SubGeometry *_geometry_p)
    {
        sub_model_p = _geomodel_p;
        seis_record_p = _seisrecord_p;
        ext_wavefield_p = _sub_ext_wavefield_p;
        sub_geometry_p = _geometry_p;
        jarvis_error_if_ptr_is_null(sub_model_p, "mpicuPhy::sub_model_p");
        jarvis_error_if_ptr_is_null(seis_record_p, "mpicuPhy::seis_record_p");
        jarvis_error_if_ptr_is_null(ext_wavefield_p, "mpicuPhy::ext_wavefield_p");
        jarvis_error_if_ptr_is_null(sub_geometry_p, "mpicuPhy::sub_geometry_p");
    }
    inline void mpicuPhy_Base::link_geomodel_para(Field<float, MemType::paged_device> *_lambda_p, Field<float, MemType::paged_device> *_mu_p, Field<float, MemType::paged_device> *_rho_p)
    {
        lambda_p = _lambda_p, mu_p = _mu_p, rho_p = _rho_p;
    }
    inline void mpicuPhy_Base::set_base_wave_mpi()
    {
        //*Assignment of static variables
        mpi_halo_vx.jarvis_mpi_cuda_stream_p = jarvis_mpi_cuda_stream_p;
        mpi_multi_halo_vel.jarvis_mpi_cuda_stream_p = jarvis_mpi_cuda_stream_p;

        int fd_half = geo_const::phy_fdorder_half;
        mpi_halo_vx.set_padding_halo_grid(mpi_multi_halo_vel, &ext_wavefield_p->vx, sub_model_p->gridphy, fd_half - 1, fd_half, fd_half - 1, fd_half, fd_half, fd_half - 1);
        mpi_halo_vy.set_padding_halo_grid(mpi_multi_halo_vel, &ext_wavefield_p->vy, sub_model_p->gridphy, fd_half - 1, fd_half, fd_half, fd_half - 1, fd_half - 1, fd_half);
        mpi_halo_vz.set_padding_halo_grid(mpi_multi_halo_vel, &ext_wavefield_p->vz, sub_model_p->gridphy, fd_half, fd_half - 1, fd_half - 1, fd_half, fd_half - 1, fd_half);
        //
        mpi_halo_sxx.set_padding_halo_grid(mpi_multi_halo_stress, &ext_wavefield_p->sxx, sub_model_p->gridphy, 0, 0, 0, 0, fd_half - 1, fd_half);
        mpi_halo_syy.set_padding_halo_grid(mpi_multi_halo_stress, &ext_wavefield_p->syy, sub_model_p->gridphy, 0, 0, fd_half - 1, fd_half, 0, 0);
        mpi_halo_szz.set_padding_halo_grid(mpi_multi_halo_stress, &ext_wavefield_p->szz, sub_model_p->gridphy, fd_half - 1, fd_half, 0, 0, 0, 0);
        //
        mpi_halo_sxy.set_padding_halo_grid(mpi_multi_halo_stress, &ext_wavefield_p->sxy, sub_model_p->gridphy, 0, 0, fd_half, fd_half - 1, fd_half, fd_half - 1);
        mpi_halo_sxz.set_padding_halo_grid(mpi_multi_halo_stress, &ext_wavefield_p->sxz, sub_model_p->gridphy, fd_half, fd_half - 1, 0, 0, fd_half, fd_half - 1);
        mpi_halo_syz.set_padding_halo_grid(mpi_multi_halo_stress, &ext_wavefield_p->syz, sub_model_p->gridphy, fd_half, fd_half - 1, fd_half, fd_half - 1, 0, 0);
    }
    inline void mpicuPhy_Base::set_base_wave_mpi_halo_wave()
    {
        mpi_halo_wave.top_vx = mpi_halo_vx.top_recv.device::ptr();
        mpi_halo_wave.top_vy = mpi_halo_vy.top_recv.device::ptr();
        mpi_halo_wave.top_vz = mpi_halo_vz.top_recv.device::ptr();
        mpi_halo_wave.top_sxz = mpi_halo_sxz.top_recv.device::ptr();
        mpi_halo_wave.top_syz = mpi_halo_syz.top_recv.device::ptr();
        mpi_halo_wave.top_szz = mpi_halo_szz.top_recv.device::ptr();
        //
        mpi_halo_wave.bottom_vx = mpi_halo_vx.bottom_recv.device::ptr();
        mpi_halo_wave.bottom_vy = mpi_halo_vy.bottom_recv.device::ptr();
        mpi_halo_wave.bottom_vz = mpi_halo_vz.bottom_recv.device::ptr();
        mpi_halo_wave.bottom_sxz = mpi_halo_sxz.bottom_recv.device::ptr();
        mpi_halo_wave.bottom_syz = mpi_halo_syz.bottom_recv.device::ptr();
        mpi_halo_wave.bottom_szz = mpi_halo_szz.bottom_recv.device::ptr();
        //
        mpi_halo_wave.front_vx = mpi_halo_vx.front_recv.device::ptr();
        mpi_halo_wave.front_vy = mpi_halo_vy.front_recv.device::ptr();
        mpi_halo_wave.front_vz = mpi_halo_vz.front_recv.device::ptr();
        mpi_halo_wave.front_sxy = mpi_halo_sxy.front_recv.device::ptr();
        mpi_halo_wave.front_syy = mpi_halo_syy.front_recv.device::ptr();
        mpi_halo_wave.front_syz = mpi_halo_syz.front_recv.device::ptr();
        //
        mpi_halo_wave.back_vx = mpi_halo_vx.back_recv.device::ptr();
        mpi_halo_wave.back_vy = mpi_halo_vy.back_recv.device::ptr();
        mpi_halo_wave.back_vz = mpi_halo_vz.back_recv.device::ptr();
        mpi_halo_wave.back_sxy = mpi_halo_sxy.back_recv.device::ptr();
        mpi_halo_wave.back_syy = mpi_halo_syy.back_recv.device::ptr();
        mpi_halo_wave.back_syz = mpi_halo_syz.back_recv.device::ptr();
        //
        mpi_halo_wave.right_vx = mpi_halo_vx.right_recv.device::ptr();
        mpi_halo_wave.right_vy = mpi_halo_vy.right_recv.device::ptr();
        mpi_halo_wave.right_vz = mpi_halo_vz.right_recv.device::ptr();
        mpi_halo_wave.right_sxx = mpi_halo_sxx.right_recv.device::ptr();
        mpi_halo_wave.right_sxy = mpi_halo_sxy.right_recv.device::ptr();
        mpi_halo_wave.right_sxz = mpi_halo_sxz.right_recv.device::ptr();
        //
        mpi_halo_wave.left_vx = mpi_halo_vx.left_recv.device::ptr();
        mpi_halo_wave.left_vy = mpi_halo_vy.left_recv.device::ptr();
        mpi_halo_wave.left_vz = mpi_halo_vz.left_recv.device::ptr();
        mpi_halo_wave.left_sxx = mpi_halo_sxx.left_recv.device::ptr();
        mpi_halo_wave.left_sxy = mpi_halo_sxy.left_recv.device::ptr();
        mpi_halo_wave.left_sxz = mpi_halo_sxz.left_recv.device::ptr();
    }
    inline void mpicuPhy_Base::set_grid()
    {
        grid_inside.set_ndl(sub_model_p->gridphy.n_rows - sub_pad.pad_left - sub_pad.pad_right,
                            sub_model_p->gridphy.n_cols - sub_pad.pad_front - sub_pad.pad_back,
                            sub_model_p->gridphy.n_slices - sub_pad.pad_top - sub_pad.pad_bottom,
                            sub_model_p->gridphy.d_rows, sub_model_p->gridphy.d_cols, sub_model_p->gridphy.d_slices,
                            sub_model_p->gridphy.l_rows + sub_pad.pad_right * sub_model_p->gridphy.d_rows,
                            sub_model_p->gridphy.l_cols + sub_pad.pad_front * sub_model_p->gridphy.d_cols,
                            sub_model_p->gridphy.l_slices + sub_pad.pad_top * sub_model_p->gridphy.d_slices);

        grid_top.set_ndl(sub_model_p->gridphy.n_rows, sub_model_p->gridphy.n_cols, sub_pad.pad_top,
                         sub_model_p->gridphy.d_rows, sub_model_p->gridphy.d_cols, sub_model_p->gridphy.d_slices,
                         sub_model_p->gridphy.l_rows,
                         sub_model_p->gridphy.l_cols,
                         sub_model_p->gridphy.l_slices);

        grid_bottom.set_ndl(sub_model_p->gridphy.n_rows, sub_model_p->gridphy.n_cols, sub_pad.pad_bottom,
                            sub_model_p->gridphy.d_rows, sub_model_p->gridphy.d_cols, sub_model_p->gridphy.d_slices,
                            sub_model_p->gridphy.l_rows,
                            sub_model_p->gridphy.l_cols,
                            sub_model_p->gridphy.r_slices - (sub_pad.pad_bottom - 1) * sub_model_p->gridphy.d_slices);

        grid_front.set_ndl(sub_model_p->gridphy.n_rows, sub_pad.pad_front, sub_model_p->gridphy.n_slices - sub_pad.pad_top - sub_pad.pad_bottom,
                           sub_model_p->gridphy.d_rows, sub_model_p->gridphy.d_cols, sub_model_p->gridphy.d_slices,
                           sub_model_p->gridphy.l_rows,
                           sub_model_p->gridphy.l_cols,
                           sub_model_p->gridphy.l_slices + sub_pad.pad_top * sub_model_p->gridphy.d_slices);

        grid_back.set_ndl(sub_model_p->gridphy.n_rows, sub_pad.pad_back, sub_model_p->gridphy.n_slices - sub_pad.pad_top - sub_pad.pad_bottom,
                          sub_model_p->gridphy.d_rows, sub_model_p->gridphy.d_cols, sub_model_p->gridphy.d_slices,
                          sub_model_p->gridphy.l_rows,
                          sub_model_p->gridphy.r_cols - (sub_pad.pad_back - 1) * sub_model_p->gridphy.d_cols,
                          sub_model_p->gridphy.l_slices + sub_pad.pad_top * sub_model_p->gridphy.d_slices);

        grid_right.set_ndl(sub_pad.pad_right, sub_model_p->gridphy.n_cols - sub_pad.pad_back - sub_pad.pad_front, sub_model_p->gridphy.n_slices - sub_pad.pad_top - sub_pad.pad_bottom,
                           sub_model_p->gridphy.d_rows, sub_model_p->gridphy.d_cols, sub_model_p->gridphy.d_slices,
                           sub_model_p->gridphy.l_rows,
                           sub_model_p->gridphy.l_cols + sub_pad.pad_front * sub_model_p->gridphy.d_cols,
                           sub_model_p->gridphy.l_slices + sub_pad.pad_top * sub_model_p->gridphy.d_slices);

        grid_left.set_ndl(sub_pad.pad_left, sub_model_p->gridphy.n_cols - sub_pad.pad_back - sub_pad.pad_front, sub_model_p->gridphy.n_slices - sub_pad.pad_top - sub_pad.pad_bottom,
                          sub_model_p->gridphy.d_rows, sub_model_p->gridphy.d_cols, sub_model_p->gridphy.d_slices,
                          sub_model_p->gridphy.r_rows - (sub_pad.pad_left - 1) * sub_model_p->gridphy.d_rows,
                          sub_model_p->gridphy.l_cols + sub_pad.pad_front * sub_model_p->gridphy.d_cols,
                          sub_model_p->gridphy.l_slices + sub_pad.pad_top * sub_model_p->gridphy.d_slices);
    }
    inline void mpicuPhy_Base::check_is_stable()
    {
        int point_num_of_T = 12;
        float fm_cal = sub_model_p->glb_vp_min / (sub_model_p->gridphy.d_rows * point_num_of_T);
        if (sub_model_p->glb_vs_min > 10)
        {
            fm_cal = sub_model_p->glb_vs_min / (sub_model_p->gridphy.d_rows * point_num_of_T);
        }
        float s = sqrt(1.0f / (sub_model_p->gridphy.d_rows * sub_model_p->gridphy.d_rows) +
                       1.0f / (sub_model_p->gridphy.d_cols * sub_model_p->gridphy.d_cols) +
                       1.0f / (sub_model_p->gridphy.d_slices * sub_model_p->gridphy.d_slices));
        Field<float, MemType::paged_device> dc;
        cu_cal_diff_coeff(dc, geo_const::phy_fdorder_half);
        float dt_cal = 1.0f / dc.abs_sum() / (s * max(sub_model_p->glb_vp_max, sub_model_p->glb_vs_max));
        dc.clear();

        if (jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank == 0)
        {
           // printf("\033[42;37m[recommend max dt]:\t\033[0m %f\n", dt_cal);
           // printf("\033[42;37m[used          dt]:\t\033[0m %f\n\n", seis_record_p->dt);
           // printf("\033[42;37m[recommend max fm]:\t\033[0m %f\n", fm_cal);
           // printf("\033[42;37m[used          fm]:\t\033[0m %f\n\n", seis_record_p->fm);
          //  printf("\033[42;37m[ntime]:\t\033[0m %d\n\n", seis_record_p->ntime);
          //  // 3
          //  if (sub_model_p->gridphy.d_rows > sub_model_p->glb_vp_min / (seis_record_p->fm * point_num_of_T) ||
            //    sub_model_p->gridphy.d_cols > sub_model_p->glb_vp_min / (seis_record_p->fm * point_num_of_T) ||
           //     sub_model_p->gridphy.d_slices > sub_model_p->glb_vp_min / (seis_record_p->fm * point_num_of_T))
           // {
            //    printf("\033[41;37mThe stability may not be enough, so it is recommended to modify the parameters!!!\033[0m\n");
            //    printf("\033[41;37mIt is recommended to reduce the grid spacing or main frequency\033[0m\n");
            //    std::abort();
            //}
          //  if (sub_model_p->glb_vs_min > 10)
           // {
             //   if (sub_model_p->gridphy.d_rows > sub_model_p->glb_vs_min / (seis_record_p->fm * point_num_of_T) ||
             //       sub_model_p->gridphy.d_cols > sub_model_p->glb_vs_min / (seis_record_p->fm * point_num_of_T) ||
             //       sub_model_p->gridphy.d_slices > sub_model_p->glb_vs_min / (seis_record_p->fm * point_num_of_T))
             //   {
            //        printf("\033[41;37mThe stability may not be enough, so it is recommended to modify the parameters!!!\033[0m\n");
            //        printf("\033[41;37mIt is recommended to reduce the grid spacing or main frequency\033[0m\n");
            //        std::abort();
           //     }
           // }
            //*Time stability
            //if (seis_record_p->dt > dt_cal)
           // {
            //    printf("\033[41;37mThe stability may not be enough, so it is recommended to modify the parameters!!!\033[0m\n");
            //    printf("\033[41;37mIt is recommended to reduce the sampling interval or increase the grid spacing\033[0m\n");
            //    std::abort();
            //}
        }

        for (int i = 0; i < jarvis_mpi_cuda_stream_p->mpi_frame.mpi_size; i++)
        {
            if (jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank == i)
            {
                sub_model_p->gridext.print_info("mpi=" + to_string(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank) + ":sub_gridext:");
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
    inline void mpicuPhy_Base::set_zero()
    {
        mpi_multi_halo_vel.set_zero();
        mpi_multi_halo_stress.set_zero();
    }
    inline void mpicuPhy_Base::clear()
    {
        mpi_multi_halo_vel.clear();
        mpi_multi_halo_stress.clear();
    }
    //
    //
    //
    //
    //
    inline void mpicuPhy<SimulateType::pure_forward>::initialize()
    {
        jarvis_error_if_ptr_is_null(jarvis_mpi_cuda_stream_p, "mpicuPhy::jarvis_mpi_cuda_stream_p");
        set_phy_frame();
        cu_cal_diff_coeff(this->phy_dc, geo_const::phy_fdorder_half);
        cu_cal_diff_coeff_list(this->phy_dc_list, geo_const::rec_width, geo_const::phy_fdorder_half);
        set_base_wave_mpi();
        mpi_multi_halo_vel.alloc_halo();
        mpi_multi_halo_stress.alloc_halo();
        mpi_multi_halo_vel.mpi_stream_init();
        mpi_multi_halo_stress.mpi_stream_init();
        set_base_wave_mpi_halo_wave();
        if (jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank == 0)
        {
            printf("current module:phy::number of all jarvis_mpi_stream: %d\n\n", jarvis_mpi_cuda_stream_p->mpi_graph.size());
        }
    }
    inline void mpicuPhy<SimulateType::pure_forward>::set_phy_frame()
    {
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_top)
            sub_pad.pad_top = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_bottom)
            sub_pad.pad_bottom = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_front)
            sub_pad.pad_front = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_back)
            sub_pad.pad_back = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_right)
            sub_pad.pad_right = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_left)
            sub_pad.pad_left = geo_const::phy_fdorder_half;
        set_grid();
    }
    //
    //
    //
    //
    //
    inline void mpicuPhy<SimulateType::rtm_forward>::set_rtm_wave_mpi()
    {
        int fd_half = geo_const::phy_fdorder_half;
        mpi_halo_sau.set_padding_halo_grid(mpi_multi_halo_stress, &ext_wavefield_p->sau, sub_model_p->gridphy, fd_half - 1, fd_half, fd_half - 1, fd_half, fd_half - 1, fd_half);
    }
    inline void mpicuPhy<SimulateType::rtm_forward>::set_rtm_wave_mpi_halo_wave()
    {
        mpi_halo_wave.top_sau = mpi_halo_sau.top_recv.device::ptr();
        mpi_halo_wave.bottom_sau = mpi_halo_sau.bottom_recv.device::ptr();
        mpi_halo_wave.front_sau = mpi_halo_sau.front_recv.device::ptr();
        mpi_halo_wave.back_sau = mpi_halo_sau.back_recv.device::ptr();
        mpi_halo_wave.right_sau = mpi_halo_sau.right_recv.device::ptr();
        mpi_halo_wave.left_sau = mpi_halo_sau.left_recv.device::ptr();
    }
    inline void mpicuPhy<SimulateType::rtm_forward>::initialize()
    {
        jarvis_error_if_ptr_is_null(jarvis_mpi_cuda_stream_p, "mpicuPhy::jarvis_mpi_cuda_stream_p");
        set_phy_frame();
        cu_cal_diff_coeff(this->phy_dc, geo_const::phy_fdorder_half);
        cu_cal_diff_coeff_list(this->phy_dc_list, geo_const::rec_width, geo_const::phy_fdorder_half);
        set_base_wave_mpi();
        set_rtm_wave_mpi();
        mpi_multi_halo_vel.alloc_halo();
        mpi_multi_halo_stress.alloc_halo();
        mpi_multi_halo_vel.mpi_stream_init();
        mpi_multi_halo_stress.mpi_stream_init();
        set_base_wave_mpi_halo_wave();
        set_rtm_wave_mpi_halo_wave();
        if (jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank == 0)
        {
            printf("current module:phy::number of all jarvis_mpi_stream: %d\n\n", jarvis_mpi_cuda_stream_p->mpi_graph.size());
        }
    }
    inline void mpicuPhy<SimulateType::rtm_forward>::set_phy_frame()
    {
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_top)
            sub_pad.pad_top = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_bottom)
            sub_pad.pad_bottom = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_front)
            sub_pad.pad_front = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_back)
            sub_pad.pad_back = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_right)
            sub_pad.pad_right = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_left)
            sub_pad.pad_left = geo_const::phy_fdorder_half;
        set_grid();
    }
    //
    //
    //
    //
    //
    inline void mpicuPhy<SimulateType::rtm_backward>::set_rtm_wave_mpi()
    {
        int fd_half = geo_const::phy_fdorder_half;
        mpi_halo_sau.set_padding_halo_grid(mpi_multi_halo_stress, &ext_wavefield_p->sau, sub_model_p->gridphy, fd_half - 1, fd_half, fd_half - 1, fd_half, fd_half - 1, fd_half);
    }
    inline void mpicuPhy<SimulateType::rtm_backward>::set_rtm_wave_mpi_halo_wave()
    {
        mpi_halo_wave.top_sau = mpi_halo_sau.top_recv.device::ptr();
        mpi_halo_wave.bottom_sau = mpi_halo_sau.bottom_recv.device::ptr();
        mpi_halo_wave.front_sau = mpi_halo_sau.front_recv.device::ptr();
        mpi_halo_wave.back_sau = mpi_halo_sau.back_recv.device::ptr();
        mpi_halo_wave.right_sau = mpi_halo_sau.right_recv.device::ptr();
        mpi_halo_wave.left_sau = mpi_halo_sau.left_recv.device::ptr();
    }
    inline void mpicuPhy<SimulateType::rtm_backward>::initialize()
    {
        jarvis_error_if_ptr_is_null(jarvis_mpi_cuda_stream_p, "mpicuPhy::jarvis_mpi_cuda_stream_p");
        set_phy_frame();
        cu_cal_diff_coeff(this->phy_dc, geo_const::phy_fdorder_half);
        cu_cal_diff_coeff_list(this->phy_dc_list, geo_const::rec_width, geo_const::phy_fdorder_half);
        set_base_wave_mpi();
        set_rtm_wave_mpi();
        mpi_multi_halo_vel.alloc_halo();
        mpi_multi_halo_stress.alloc_halo();
        mpi_multi_halo_vel.mpi_stream_init();
        mpi_multi_halo_stress.mpi_stream_init();
        set_base_wave_mpi_halo_wave();
        set_rtm_wave_mpi_halo_wave();
        if (jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank == 0)
        {
            printf("current module:phy::number of all jarvis_mpi_stream: %d\n\n", jarvis_mpi_cuda_stream_p->mpi_graph.size());
        }
    }
    inline void mpicuPhy<SimulateType::rtm_backward>::set_phy_frame()
    {
        sub_pad.pad_top = geo_const::phy_fdorder_half;
        sub_pad.pad_bottom = geo_const::phy_fdorder_half;
        sub_pad.pad_front = geo_const::phy_fdorder_half;
        sub_pad.pad_back = geo_const::phy_fdorder_half;
        sub_pad.pad_right = geo_const::phy_fdorder_half;
        sub_pad.pad_left = geo_const::phy_fdorder_half;
        set_grid();
    }
    //
    //
    //
    //
    //
    inline void mpicuPhy<SimulateType::rtm_reverse>::set_rtm_wave_mpi()
    {
        int fd_half = geo_const::phy_fdorder_half;
        mpi_halo_sau.set_padding_halo_grid(mpi_multi_halo_stress, &ext_wavefield_p->sau, sub_model_p->gridphy, fd_half - 1, fd_half, fd_half - 1, fd_half, fd_half - 1, fd_half);
    }
    inline void mpicuPhy<SimulateType::rtm_reverse>::set_rtm_wave_mpi_halo_wave()
    {
        mpi_halo_wave.top_sau = mpi_halo_sau.top_recv.device::ptr();
        mpi_halo_wave.bottom_sau = mpi_halo_sau.bottom_recv.device::ptr();
        mpi_halo_wave.front_sau = mpi_halo_sau.front_recv.device::ptr();
        mpi_halo_wave.back_sau = mpi_halo_sau.back_recv.device::ptr();
        mpi_halo_wave.right_sau = mpi_halo_sau.right_recv.device::ptr();
        mpi_halo_wave.left_sau = mpi_halo_sau.left_recv.device::ptr();
    }
    inline void mpicuPhy<SimulateType::rtm_reverse>::initialize()
    {
        jarvis_error_if_ptr_is_null(jarvis_mpi_cuda_stream_p, "mpicuPhy::jarvis_mpi_cuda_stream_p");
        set_phy_frame();
        cu_cal_diff_coeff(this->phy_dc, geo_const::phy_fdorder_half);
        cu_cal_diff_coeff_list(this->phy_dc_list, geo_const::rec_width, geo_const::phy_fdorder_half);
        set_base_wave_mpi();
        set_rtm_wave_mpi();
        mpi_multi_halo_vel.alloc_halo();
        mpi_multi_halo_stress.alloc_halo();
        mpi_multi_halo_vel.mpi_stream_init();
        mpi_multi_halo_stress.mpi_stream_init();
        set_base_wave_mpi_halo_wave();
        set_rtm_wave_mpi_halo_wave();
        if (jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank == 0)
        {
            printf("current module:phy::number of all jarvis_mpi_stream: %d\n\n", jarvis_mpi_cuda_stream_p->mpi_graph.size());
        }
    }
    inline void mpicuPhy<SimulateType::rtm_reverse>::set_phy_frame()
    {
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_top)
            sub_pad.pad_top = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_bottom)
            sub_pad.pad_bottom = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_front)
            sub_pad.pad_front = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_back)
            sub_pad.pad_back = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_right)
            sub_pad.pad_right = geo_const::phy_fdorder_half;
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_left)
            sub_pad.pad_left = geo_const::phy_fdorder_half;
        set_grid();
    }
}
#endif