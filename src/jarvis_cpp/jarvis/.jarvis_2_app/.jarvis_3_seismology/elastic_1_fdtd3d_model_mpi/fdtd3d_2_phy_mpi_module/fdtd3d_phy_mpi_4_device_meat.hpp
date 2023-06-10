#pragma once
#ifndef _FDTD3D_PHY_MPI_DEVICE_MEAT_HPP
#define _FDTD3D_PHY_MPI_DEVICE_MEAT_HPP
#include "fdtd3d_phy_mpi_3_stress_kernel.hpp"
namespace jarvis
{
    inline void mpicuPhy_Base::update_inner_vel(SimulateType _simulate_type, float _dt)
    {
        void *args_list[] = {&_simulate_type,
                             &sub_model_p->gridext,
                             &grid_inside,
                             &_dt,
                             &phy_dc.device::ptr(),
                             &rho_p->device::ptr(),
                             &ext_wavefield_p->org_wave,
                             &ext_wavefield_p->rtm_wave,
                             &ext_wavefield_p->fwi_wave};
        cudaLaunchKernel((void *)sub_update_vel_phy_domain_inside, jarvis_cuda_kernel_size(grid_inside.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
    }
    inline void mpicuPhy_Base::update_inner_stress(SimulateType _simulate_type, float _dt)
    {
        void *args_list[] = {&_simulate_type,
                             &sub_model_p->gridext,
                             &grid_inside,
                             &_dt,
                             &phy_dc.device::ptr(),
                             &lambda_p->device::ptr(),
                             &mu_p->device::ptr(),
                             &ext_wavefield_p->org_wave,
                             &ext_wavefield_p->rtm_wave};
        cudaLaunchKernel((void *)sub_update_stress_phy_domain_inside, jarvis_cuda_kernel_size(grid_inside.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
    }
    inline void mpicuPhy_Base::update_outer_vel(SimulateType _simulate_type, int _start_fdorder, float _dt)
    {
        void *args_list[] = {&_simulate_type,
                             &sub_model_p->gridext,
                             &sub_model_p->gridphy,
                             &grid_top,
                             &sub_pad,
                             &jarvis_mpi_cuda_stream_p->mpi_frame.near,
                             &_start_fdorder,
                             &_dt,
                             &phy_dc.device::ptr(),
                             &phy_dc_list.device::ptr(),
                             &rho_p->device::ptr(),
                             &ext_wavefield_p->org_wave,
                             &ext_wavefield_p->rtm_wave,
                             &ext_wavefield_p->fwi_wave,
                             &mpi_halo_wave};
        if (sub_pad.pad_top == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_top;
            cudaLaunchKernel((void *)sub_update_vel_phy_domain_edge_top, jarvis_cuda_kernel_size(grid_top.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
        if (sub_pad.pad_bottom == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_bottom;
            cudaLaunchKernel((void *)sub_update_vel_phy_domain_edge_bottom, jarvis_cuda_kernel_size(grid_bottom.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
        if (sub_pad.pad_front == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_front;
            cudaLaunchKernel((void *)sub_update_vel_phy_domain_edge_front, jarvis_cuda_kernel_size(grid_front.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
        if (sub_pad.pad_back == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_back;
            cudaLaunchKernel((void *)sub_update_vel_phy_domain_edge_back, jarvis_cuda_kernel_size(grid_back.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
        if (sub_pad.pad_right == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_right;
            cudaLaunchKernel((void *)sub_update_vel_phy_domain_edge_right, jarvis_cuda_kernel_size(grid_right.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
        if (sub_pad.pad_left == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_left;
            cudaLaunchKernel((void *)sub_update_vel_phy_domain_edge_left, jarvis_cuda_kernel_size(grid_left.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
    }
    inline void mpicuPhy_Base::update_outer_stress(SimulateType _simulate_type, int _start_fdorder, float _dt)
    {
        void *args_list[] = {&_simulate_type,
                             &sub_model_p->gridext,
                             &sub_model_p->gridphy,
                             &grid_top,
                             &sub_pad,
                             &jarvis_mpi_cuda_stream_p->mpi_frame.near,
                             &_start_fdorder,
                             &_dt,
                             &phy_dc.device::ptr(),
                             &phy_dc_list.device::ptr(),
                             &lambda_p->device::ptr(),
                             &mu_p->device::ptr(),
                             &ext_wavefield_p->org_wave,
                             &ext_wavefield_p->rtm_wave,
                             &mpi_halo_wave};
        if (sub_pad.pad_top == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_top;
            cudaLaunchKernel((void *)sub_update_stress_phy_domain_edge_top, jarvis_cuda_kernel_size(grid_top.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
        if (sub_pad.pad_bottom == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_bottom;
            cudaLaunchKernel((void *)sub_update_stress_phy_domain_edge_bottom, jarvis_cuda_kernel_size(grid_bottom.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
        if (sub_pad.pad_front == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_front;
            cudaLaunchKernel((void *)sub_update_stress_phy_domain_edge_front, jarvis_cuda_kernel_size(grid_front.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
        if (sub_pad.pad_back == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_back;
            cudaLaunchKernel((void *)sub_update_stress_phy_domain_edge_back, jarvis_cuda_kernel_size(grid_back.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
        if (sub_pad.pad_right == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_right;
            cudaLaunchKernel((void *)sub_update_stress_phy_domain_edge_right, jarvis_cuda_kernel_size(grid_right.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
        if (sub_pad.pad_left == geo_const::phy_fdorder_half)
        {
            args_list[3] = &grid_left;
            cudaLaunchKernel((void *)sub_update_stress_phy_domain_edge_left, jarvis_cuda_kernel_size(grid_left.n_elem), args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
    }
    //
    //
    //
    //
    //
    inline void mpicuPhy<SimulateType::pure_forward>::cuda_update_inside_vel()
    {
        update_inner_vel(SimulateType::pure_forward, seis_record_p->dt);
    }
    inline void mpicuPhy<SimulateType::pure_forward>::cuda_update_inside_stress()
    {
        update_inner_stress(SimulateType::pure_forward, seis_record_p->dt);
    }
    inline void mpicuPhy<SimulateType::pure_forward>::cuda_update_and_mpi_exchange_edge_vel()
    {
        update_outer_vel(SimulateType::pure_forward, geo_const::phy_fdorder_half, seis_record_p->dt);
        mpi_halo_vx.extract_into_halo();
        mpi_halo_vy.extract_into_halo();
        mpi_halo_vz.extract_into_halo();
        jarvis_mpi_cuda_stream_p->sync_cal_stream();
        jarvis_mpi_cuda_stream_p->mpi_start_batch_id(mpi_multi_halo_vel.operator_id);
    }
    inline void mpicuPhy<SimulateType::pure_forward>::cuda_update_and_mpi_exchange_edge_stress()
    {
        update_outer_stress(SimulateType::pure_forward, geo_const::phy_fdorder_half, seis_record_p->dt);
        mpi_halo_sxx.extract_into_halo();
        mpi_halo_syy.extract_into_halo();
        mpi_halo_szz.extract_into_halo();
        mpi_halo_sxy.extract_into_halo();
        mpi_halo_sxz.extract_into_halo();
        mpi_halo_syz.extract_into_halo();
        cudaStreamSynchronize(jarvis_mpi_cuda_stream_p->cal_stream());
        jarvis_mpi_cuda_stream_p->mpi_start_batch_id(mpi_multi_halo_stress.operator_id);
    }
    inline void mpicuPhy<SimulateType::pure_forward>::cuda_add_stress_source(int _it, int _i_shot)
    {
        void *args_list[] = {&_it,
                             &_i_shot,
                             &sub_model_p->gridext,
                             &sub_geometry_p->shot.vector<Point3D, MemType::device>::ptr(),
                             &seis_record_p->wave_let.vector<float, MemType::device>::ptr(),
                             &ext_wavefield_p->sxx.ptr(),
                             &ext_wavefield_p->syy.ptr(),
                             &ext_wavefield_p->szz.ptr()};
        cudaLaunchKernel((void *)sub_add_stress_source, 1, 1, args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
    }
    //
    //
    //
    //
    //
    inline void mpicuPhy<SimulateType::rtm_forward>::cuda_update_inside_vel()
    {
        update_inner_vel(SimulateType::rtm_forward, seis_record_p->dt);
    }
    inline void mpicuPhy<SimulateType::rtm_forward>::cuda_update_inside_stress()
    {
        update_inner_stress(SimulateType::rtm_forward, seis_record_p->dt);
    }
    inline void mpicuPhy<SimulateType::rtm_forward>::cuda_update_and_mpi_exchange_edge_vel()
    {
        update_outer_vel(SimulateType::rtm_forward, geo_const::phy_fdorder_half, seis_record_p->dt);
        mpi_halo_vx.extract_into_halo();
        mpi_halo_vy.extract_into_halo();
        mpi_halo_vz.extract_into_halo();
        jarvis_mpi_cuda_stream_p->sync_cal_stream();
        jarvis_mpi_cuda_stream_p->mpi_start_batch_id(mpi_multi_halo_vel.operator_id);
    }
    inline void mpicuPhy<SimulateType::rtm_forward>::cuda_update_and_mpi_exchange_edge_stress()
    {
        update_outer_stress(SimulateType::rtm_forward, geo_const::phy_fdorder_half, seis_record_p->dt);
        mpi_halo_sxx.extract_into_halo();
        mpi_halo_syy.extract_into_halo();
        mpi_halo_szz.extract_into_halo();
        mpi_halo_sxy.extract_into_halo();
        mpi_halo_sxz.extract_into_halo();
        mpi_halo_syz.extract_into_halo();
        mpi_halo_sau.extract_into_halo();
        cudaStreamSynchronize(jarvis_mpi_cuda_stream_p->cal_stream());
        jarvis_mpi_cuda_stream_p->mpi_start_batch_id(mpi_multi_halo_stress.operator_id);
    }
    inline void mpicuPhy<SimulateType::rtm_forward>::cuda_add_stress_source(int _it, int _i_shot)
    {
        void *args_list[] = {&_it,
                             &_i_shot,
                             &sub_model_p->gridext,
                             &sub_geometry_p->shot.vector<Point3D, MemType::device>::ptr(),
                             &seis_record_p->wave_let.vector<float, MemType::device>::ptr(),
                             &ext_wavefield_p->sxx.ptr(),
                             &ext_wavefield_p->syy.ptr(),
                             &ext_wavefield_p->szz.ptr()};
        cudaLaunchKernel((void *)sub_add_stress_source, 1, 1, args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
    }
    //
    //
    //
    //
    //
    inline void mpicuPhy<SimulateType::rtm_backward>::cuda_update_inside_vel()
    {
        update_inner_vel(SimulateType::rtm_backward, -seis_record_p->dt);
    }
    inline void mpicuPhy<SimulateType::rtm_backward>::cuda_update_inside_stress()
    {
        update_inner_stress(SimulateType::rtm_backward, -seis_record_p->dt);
    }
    inline void mpicuPhy<SimulateType::rtm_backward>::cuda_update_and_mpi_exchange_edge_vel()
    {
        update_outer_vel(SimulateType::rtm_backward, geo_const::rec_width, -seis_record_p->dt);
        mpi_halo_vx.extract_into_halo();
        mpi_halo_vy.extract_into_halo();
        mpi_halo_vz.extract_into_halo();
        jarvis_mpi_cuda_stream_p->sync_cal_stream();
        jarvis_mpi_cuda_stream_p->mpi_start_batch_id(mpi_multi_halo_vel.operator_id);
    }
    inline void mpicuPhy<SimulateType::rtm_backward>::cuda_update_and_mpi_exchange_edge_stress()
    {
        update_outer_stress(SimulateType::rtm_backward, geo_const::rec_width, -seis_record_p->dt);
        mpi_halo_sxx.extract_into_halo();
        mpi_halo_syy.extract_into_halo();
        mpi_halo_szz.extract_into_halo();
        mpi_halo_sxy.extract_into_halo();
        mpi_halo_sxz.extract_into_halo();
        mpi_halo_syz.extract_into_halo();
        mpi_halo_sau.extract_into_halo();
        cudaStreamSynchronize(jarvis_mpi_cuda_stream_p->cal_stream());
        jarvis_mpi_cuda_stream_p->mpi_start_batch_id(mpi_multi_halo_stress.operator_id);
    }
    inline void mpicuPhy<SimulateType::rtm_backward>::cuda_subtract_stress_source(int _it, int _i_shot)
    {
        void *args_list[] = {&_it,
                             &_i_shot,
                             &sub_model_p->gridext,
                             &sub_geometry_p->shot.vector<Point3D, MemType::device>::ptr(),
                             &seis_record_p->wave_let.vector<float, MemType::device>::ptr(),
                             &ext_wavefield_p->sxx.ptr(),
                             &ext_wavefield_p->syy.ptr(),
                             &ext_wavefield_p->szz.ptr()};
        cudaLaunchKernel((void *)sub_subtract_stress_source, 1, 1, args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
    }
    //
    //
    //
    //
    //
    inline void mpicuPhy<SimulateType::rtm_reverse>::cuda_update_inside_vel()
    {
        update_inner_vel(SimulateType::rtm_reverse, seis_record_p->dt);
    }
    inline void mpicuPhy<SimulateType::rtm_reverse>::cuda_update_inside_stress()
    {
        update_inner_stress(SimulateType::rtm_reverse, seis_record_p->dt);
    }
    inline void mpicuPhy<SimulateType::rtm_reverse>::cuda_update_and_mpi_exchange_edge_vel()
    {
        update_outer_vel(SimulateType::rtm_reverse, geo_const::phy_fdorder_half, seis_record_p->dt);
        mpi_halo_vx.extract_into_halo();
        mpi_halo_vy.extract_into_halo();
        mpi_halo_vz.extract_into_halo();
        jarvis_mpi_cuda_stream_p->sync_cal_stream();
        jarvis_mpi_cuda_stream_p->mpi_start_batch_id(mpi_multi_halo_vel.operator_id);
    }
    inline void mpicuPhy<SimulateType::rtm_reverse>::cuda_update_and_mpi_exchange_edge_stress()
    {
        update_outer_stress(SimulateType::rtm_reverse, geo_const::phy_fdorder_half, seis_record_p->dt);
        mpi_halo_sxx.extract_into_halo();
        mpi_halo_syy.extract_into_halo();
        mpi_halo_szz.extract_into_halo();
        mpi_halo_sxy.extract_into_halo();
        mpi_halo_sxz.extract_into_halo();
        mpi_halo_syz.extract_into_halo();
        mpi_halo_sau.extract_into_halo();
        cudaStreamSynchronize(jarvis_mpi_cuda_stream_p->cal_stream());
        jarvis_mpi_cuda_stream_p->mpi_start_batch_id(mpi_multi_halo_stress.operator_id);
    }
    inline void mpicuPhy<SimulateType::rtm_reverse>::cuda_rtm_reverse_add_vel_source(int it, SeisRecord *raw_seis)
    {
        if (sub_geometry_p->is_have_recv)
        {
            void *args_list[] = {&sub_model_p->gridext,
                                 &it,
                                 &seis_record_p->ntime,
                                 &sub_geometry_p->sub_recv.n_elem,
                                 &sub_geometry_p->sub_recv.vector<Point3D, MemType::device>::ptr(),
                                 &ext_wavefield_p->vx.ptr(),
                                 &ext_wavefield_p->vy.ptr(),
                                 &ext_wavefield_p->vz.ptr(),
                                 &raw_seis->vx_seis.ptr(),
                                 &raw_seis->vy_seis.ptr(),
                                 &raw_seis->vz_seis.ptr()};
            cudaLaunchKernel((void *)sub_rtm_reverse_add_vel_source, 1, sub_geometry_p->sub_recv.n_elem, args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
    }
}
#endif