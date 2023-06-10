#pragma once
#ifndef _FDTD3D_PML_MEAT_HPP
#define _FDTD3D_PML_MEAT_HPP
#include "fdtd3d_base_6_pml_bones.hpp"
namespace jarvis
{
    inline void Padding::print()
    {
        printf("pad_top: %d\n", pad_top);
        printf("pad_bottom: %d\n", pad_bottom);
        printf("pad_front: %d\n", pad_front);
        printf("pad_back: %d\n", pad_back);
        printf("pad_right: %d\n", pad_right);
        printf("pad_left: %d\n", pad_left);
    }
    //
    inline void cuPMLBase::link(elasticGeoModel<domain::local> *_geomodel_p, SeisRecord *_seisrecord_p,
                                elasticWave *_sub_ext_wavefield_p)
    {
        geomodel_p = _geomodel_p, seis_record_p = _seisrecord_p;
        ext_wavefield_p = _sub_ext_wavefield_p;
    }
    //
    inline void cuPMLBase::link_geomodel_para(Field<float, MemType::paged_device> *_lambda_p, Field<float, MemType::paged_device> *_mu_p, Field<float, MemType::paged_device> *_rho_p)
    {
        lambda_p = _lambda_p, mu_p = _mu_p, rho_p = _rho_p;
    }
    //
    inline void cuPMLBase::set_grid_dir()
    {
        jarvis_error_if_ptr_is_null(geomodel_p, "PML::geomodel_p");
        jarvis_error_if_ptr_is_null(seis_record_p, "PML::seis_record_p");
        jarvis_error_if_ptr_is_null(ext_wavefield_p, "PML::ext_wavefield_p");
        jarvis_error_if_ptr_is_null(jarvis_mpi_cuda_stream_p, "PML::jarvis_mpi_cuda_stream_p");

        elasticGeoModel<domain::local> &model = *geomodel_p;
        if (is_pml.is_top)
        {
            grid_top.set_ndl(model.gridext.n_rows, model.gridext.n_cols, model.margin.top_margin - geo_const::pml_fdorder_half,
                             model.gridext.d_rows, model.gridext.d_cols, model.gridext.d_slices,
                             model.gridext.l_rows,
                             model.gridext.l_cols,
                             model.gridext.l_slices + geo_const::pml_fdorder_half * model.gridext.d_slices);
        }
        if (is_pml.is_bottom)
        {
            grid_bottom.set_ndl(model.gridext.n_rows, model.gridext.n_cols, model.margin.bottom_margin - geo_const::pml_fdorder_half,
                                model.gridext.d_rows, model.gridext.d_cols, model.gridext.d_slices,
                                model.gridext.l_rows,
                                model.gridext.l_cols,
                                model.gridphy.r_slices + model.gridphy.d_slices);
        }
        if (is_pml.is_front)
        {
            grid_front.set_ndl(model.gridext.n_rows, model.margin.front_margin - geo_const::pml_fdorder_half, model.gridphy.n_slices,
                               model.gridext.d_rows, model.gridext.d_cols, model.gridext.d_slices,
                               model.gridext.l_rows,
                               model.gridext.l_cols + geo_const::pml_fdorder_half * model.gridext.d_cols,
                               model.gridphy.l_slices);
        }
        if (is_pml.is_back)
        {
            grid_back.set_ndl(model.gridext.n_rows, model.margin.back_margin - geo_const::pml_fdorder_half, model.gridphy.n_slices,
                              model.gridext.d_rows, model.gridext.d_cols, model.gridext.d_slices,
                              model.gridext.l_rows,
                              model.gridphy.r_cols + model.gridphy.d_cols,
                              model.gridphy.l_slices);
        }
        if (is_pml.is_left)
        {
            grid_left.set_ndl(model.margin.left_margin - geo_const::pml_fdorder_half, model.gridphy.n_cols, model.gridphy.n_slices,
                              model.gridext.d_rows, model.gridext.d_cols, model.gridext.d_slices,
                              model.gridphy.r_rows + model.gridphy.d_rows,
                              model.gridphy.l_cols,
                              model.gridphy.l_slices);
        }
        if (is_pml.is_right)
        {
            grid_right.set_ndl(model.margin.right_margin - geo_const::pml_fdorder_half, model.gridphy.n_cols, model.gridphy.n_slices,
                               model.gridext.d_rows, model.gridext.d_cols, model.gridext.d_slices,
                               model.gridext.l_rows + geo_const::pml_fdorder_half * model.gridext.d_rows,
                               model.gridphy.l_cols,
                               model.gridphy.l_slices);
        }
        if (jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank == 0)
        {
            grid_top.print_info("top");
            grid_bottom.print_info("bottom");
            grid_front.print_info("front");
            grid_back.print_info("back");
            grid_right.print_info("right");
            grid_left.print_info("left");
        }
    }
    //
    inline void cuPMLBase::cu_alloc_only_device_memory_wave_dir() // 2
    {
        if (is_pml.is_top)
            top.alloc_for_z_dir(grid_top);
        if (is_pml.is_bottom)
            bottom.alloc_for_z_dir(grid_bottom);
        if (is_pml.is_front)
            front.alloc_for_y_dir(grid_front);
        if (is_pml.is_back)
            back.alloc_for_y_dir(grid_back);
        if (is_pml.is_left)
            left.alloc_for_x_dir(grid_left);
        if (is_pml.is_right)
            right.alloc_for_x_dir(grid_right);
    }

    inline void cuPMLBase::print_is_pml(int mpi_rank)
    {
        std::cout << "mpi_rank_" << to_string(mpi_rank) << ":pml_up:" << is_pml.is_top << std::endl;
        std::cout << "mpi_rank_" << to_string(mpi_rank) << ":pml_down:" << is_pml.is_bottom << std::endl;
        std::cout << "mpi_rank_" << to_string(mpi_rank) << ":pml_left:" << is_pml.is_left << std::endl;
        std::cout << "mpi_rank_" << to_string(mpi_rank) << ":pml_right:" << is_pml.is_right << std::endl;
        std::cout << "mpi_rank_" << to_string(mpi_rank) << ":pml_front:" << is_pml.is_front << std::endl;
        std::cout << "mpi_rank_" << to_string(mpi_rank) << ":pml_back:" << is_pml.is_back << std::endl;
    }
    //
    inline void mpicuPMLBase::mpicuHalo::batch_set_pad_halo_grid(Position pos_flag, Frame &grid_pml, elasticWave *ext_wavefield_p,
                                                                 mpicu_halo_comm_Operator<float, MemBlock::multiple> &_mpi_multi_halo_vel,
                                                                 mpicu_halo_comm_Operator<float, MemBlock::multiple> &_mpi_multi_halo_stress)
    {
        mpi_multi_halo_vel_p = &_mpi_multi_halo_vel;
        mpi_multi_halo_stress_p = &_mpi_multi_halo_stress;
        int fd_half = geo_const::pml_fdorder_half;
        if (pos_flag == Position::top || pos_flag == Position::bottom)
        {
            halo_vx.set_padding_halo_grid(*mpi_multi_halo_vel_p, &ext_wavefield_p->vx, grid_pml, 0, 0, fd_half - 1, fd_half, fd_half, fd_half - 1);
            halo_vy.set_padding_halo_grid(*mpi_multi_halo_vel_p, &ext_wavefield_p->vy, grid_pml, 0, 0, fd_half, fd_half - 1, fd_half - 1, fd_half);
            halo_vz.set_padding_halo_grid(*mpi_multi_halo_vel_p, &ext_wavefield_p->vz, grid_pml, 0, 0, fd_half - 1, fd_half, fd_half - 1, fd_half);
            //
            halo_sxx.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->sxx, grid_pml, 0, 0, 0, 0, fd_half - 1, fd_half);
            halo_syy.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->syy, grid_pml, 0, 0, fd_half - 1, fd_half, 0, 0);
            halo_szz.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->szz, grid_pml, 0, 0, 0, 0, 0, 0);
            halo_sxy.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->sxy, grid_pml, 0, 0, fd_half, fd_half - 1, fd_half, fd_half - 1);
            halo_sxz.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->sxz, grid_pml, 0, 0, 0, 0, fd_half, fd_half - 1);
            halo_syz.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->syz, grid_pml, 0, 0, fd_half, fd_half - 1, 0, 0);
        }
        else if (pos_flag == Position::front || pos_flag == Position::back)
        {
            halo_vx.set_padding_halo_grid(*mpi_multi_halo_vel_p, &ext_wavefield_p->vx, grid_pml, fd_half - 1, fd_half, 0, 0, fd_half, fd_half - 1);
            halo_vy.set_padding_halo_grid(*mpi_multi_halo_vel_p, &ext_wavefield_p->vy, grid_pml, fd_half - 1, fd_half, 0, 0, fd_half - 1, fd_half);
            halo_vz.set_padding_halo_grid(*mpi_multi_halo_vel_p, &ext_wavefield_p->vz, grid_pml, fd_half, fd_half - 1, 0, 0, fd_half - 1, fd_half);
            //
            halo_sxx.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->sxx, grid_pml, 0, 0, 0, 0, fd_half - 1, fd_half);
            halo_syy.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->syy, grid_pml, 0, 0, 0, 0, 0, 0);
            halo_szz.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->szz, grid_pml, fd_half - 1, fd_half, 0, 0, 0, 0);
            halo_sxy.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->sxy, grid_pml, 0, 0, 0, 0, fd_half, fd_half - 1);
            halo_sxz.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->sxz, grid_pml, fd_half, fd_half - 1, 0, 0, fd_half, fd_half - 1);
            halo_syz.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->syz, grid_pml, fd_half, fd_half - 1, 0, 0, 0, 0);
        }
        else if (pos_flag == Position::right || pos_flag == Position::left)
        {
            halo_vx.set_padding_halo_grid(*mpi_multi_halo_vel_p, &ext_wavefield_p->vx, grid_pml, fd_half - 1, fd_half, fd_half - 1, fd_half, 0, 0);
            halo_vy.set_padding_halo_grid(*mpi_multi_halo_vel_p, &ext_wavefield_p->vy, grid_pml, fd_half - 1, fd_half, fd_half, fd_half - 1, 0, 0);
            halo_vz.set_padding_halo_grid(*mpi_multi_halo_vel_p, &ext_wavefield_p->vz, grid_pml, fd_half, fd_half - 1, fd_half - 1, fd_half, 0, 0);
            //
            halo_sxx.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->sxx, grid_pml, 0, 0, 0, 0, 0, 0);
            halo_syy.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->syy, grid_pml, 0, 0, fd_half - 1, fd_half, 0, 0);
            halo_szz.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->szz, grid_pml, fd_half - 1, fd_half, 0, 0, 0, 0);
            halo_sxy.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->sxy, grid_pml, 0, 0, fd_half, fd_half - 1, 0, 0);
            halo_sxz.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->sxz, grid_pml, fd_half, fd_half - 1, 0, 0, 0, 0);
            halo_syz.set_padding_halo_grid(*mpi_multi_halo_stress_p, &ext_wavefield_p->syz, grid_pml, fd_half, fd_half - 1, fd_half, fd_half - 1, 0, 0);
        }
    }
    //
    inline void mpicuPMLBase::initialize()
    {
        set_is_pml();
        set_grid_dir();
        cu_alloc_only_device_memory_wave_dir();
        cu_alloc_damp_coeff();
        calculate_damp_coeff();
        cu_copy_damp_coeff_h2d();
        set_sub_padding_and_grid();
        cu_cal_diff_coeff(pml_dc, geo_const::pml_fdorder_half);
        mpi_init();
    }
    //
    inline void mpicuPMLBase::mpi_init()
    {
        top_pml_halo.halo_vx.jarvis_mpi_cuda_stream_p = jarvis_mpi_cuda_stream_p;
        mpi_multi_halo_vel.jarvis_mpi_cuda_stream_p = jarvis_mpi_cuda_stream_p;

        top_pml_halo.batch_set_pad_halo_grid(Position::top, grid_top, ext_wavefield_p, mpi_multi_halo_vel, mpi_multi_halo_stress);
        bottom_pml_halo.batch_set_pad_halo_grid(Position::bottom, grid_bottom, ext_wavefield_p, mpi_multi_halo_vel, mpi_multi_halo_stress);
        front_pml_halo.batch_set_pad_halo_grid(Position::front, grid_front, ext_wavefield_p, mpi_multi_halo_vel, mpi_multi_halo_stress);
        back_pml_halo.batch_set_pad_halo_grid(Position::back, grid_back, ext_wavefield_p, mpi_multi_halo_vel, mpi_multi_halo_stress);
        right_pml_halo.batch_set_pad_halo_grid(Position::right, grid_right, ext_wavefield_p, mpi_multi_halo_vel, mpi_multi_halo_stress);
        left_pml_halo.batch_set_pad_halo_grid(Position::left, grid_left, ext_wavefield_p, mpi_multi_halo_vel, mpi_multi_halo_stress);
        //
        mpi_multi_halo_vel.alloc_halo();
        mpi_multi_halo_stress.alloc_halo();

        mpi_multi_halo_vel.mpi_stream_init();
        mpi_multi_halo_stress.mpi_stream_init();

        if (jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank == 0)
        {
            printf("current module:pml::number of all jarvis_mpi_stream: %d\n\n", jarvis_mpi_cuda_stream_p->mpi_graph.size());
        }
    }
    inline void mpicuPMLBase::set_sub_padding_and_grid()
    {
        if (geomodel_p->margin.front_margin >= geo_const::pml_fdorder_half)
            sub_pad.pad_front = geo_const::pml_fdorder_half;
        if (geomodel_p->margin.back_margin >= geo_const::pml_fdorder_half)
            sub_pad.pad_back = geo_const::pml_fdorder_half;
        if (geomodel_p->margin.right_margin >= geo_const::pml_fdorder_half)
            sub_pad.pad_right = geo_const::pml_fdorder_half;
        if (geomodel_p->margin.left_margin >= geo_const::pml_fdorder_half)
            sub_pad.pad_left = geo_const::pml_fdorder_half;
    }
    //
    inline void mpicuPMLBase::set_is_pml()
    {
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_top)
            this->is_pml.is_top = false;
        else
            this->is_pml.is_top = true;
        //
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_bottom)
            this->is_pml.is_bottom = false;
        else
            this->is_pml.is_bottom = true;
        //
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_front)
            this->is_pml.is_front = false;
        else
            this->is_pml.is_front = true;
        //
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_back)
            this->is_pml.is_back = false;
        else
            this->is_pml.is_back = true;
        //
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_right)
            this->is_pml.is_right = false;
        else
            this->is_pml.is_right = true;
        //
        if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_left)
            this->is_pml.is_left = false;
        else
            this->is_pml.is_left = true;
    }
    //
    inline void mpicuPMLBase::set_zero()
    {
        if (this->is_pml.is_top)
            this->top.set_zero();
        if (this->is_pml.is_bottom)
            this->bottom.set_zero();
        if (this->is_pml.is_front)
            this->front.set_zero();
        if (this->is_pml.is_back)
            this->back.set_zero();
        if (this->is_pml.is_right)
            this->right.set_zero();
        if (this->is_pml.is_left)
            this->left.set_zero();
        mpi_multi_halo_vel.set_zero();
        mpi_multi_halo_stress.set_zero();
    }
    //
    inline void mpicuPMLBase::mpicuHalo::batch_extract_halo_vel()
    {
        halo_vx.extract_into_halo();
        halo_vy.extract_into_halo();
        halo_vz.extract_into_halo();
    }
    //
    inline void mpicuPMLBase::mpicuHalo::batch_extract_halo_stress()
    {
        halo_sxx.extract_into_halo();
        halo_syy.extract_into_halo();
        halo_szz.extract_into_halo();
        halo_sxy.extract_into_halo();
        halo_sxz.extract_into_halo();
        halo_syz.extract_into_halo();
    }
}
#endif