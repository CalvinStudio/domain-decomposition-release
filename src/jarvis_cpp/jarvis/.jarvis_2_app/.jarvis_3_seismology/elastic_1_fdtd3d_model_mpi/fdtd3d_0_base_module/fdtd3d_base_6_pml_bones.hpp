#pragma once
#ifndef _FDTD3D_PML_BONES_HPP
#define _FDTD3D_PML_BONES_HPP
#include "fdtd3d_base_5_wave_recon_meat.hpp"
namespace jarvis
{
    struct Padding
    {
        uint8_t pad_top = 0, pad_bottom = 0;
        uint8_t pad_front = 0, pad_back = 0;
        uint8_t pad_right = 0, pad_left = 0;
        void print();
    };
    enum class Position
    {
        top = 1,
        bottom,
        left,
        right,
        front,
        back
    };
    class cuPMLBase : public elastic_fdtd3d_module_Base
    {
    public:
        elasticGeoModel<domain::local> *geomodel_p = nullptr;
        SeisRecord *seis_record_p = nullptr;
        elasticWave *ext_wavefield_p = nullptr;
        //
        Field<float, MemType::paged_device> *lambda_p = nullptr;
        Field<float, MemType::paged_device> *mu_p = nullptr;
        Field<float, MemType::paged_device> *rho_p = nullptr;
        //
        int pml_num;
        isPosition is_pml;
        //
        Frame grid_top;
        Frame grid_bottom;
        Frame grid_front;
        Frame grid_back;
        Frame grid_right;
        Frame grid_left;
        //
        MemoryWave top;
        MemoryWave bottom;
        MemoryWave front;
        MemoryWave back;
        MemoryWave right;
        MemoryWave left;
        //
        Field<float, MemType::paged_device> pml_dc;
        struct host_damp
        {
            Field<float, MemType::paged, MemBlock::puppet> top;
            Field<float, MemType::paged, MemBlock::puppet> bottom;
            Field<float, MemType::paged, MemBlock::puppet> front;
            Field<float, MemType::paged, MemBlock::puppet> back;
            Field<float, MemType::paged, MemBlock::puppet> right;
            Field<float, MemType::paged, MemBlock::puppet> left;
        };
        struct device_damp
        {
            Field<float, MemType::paged_device, MemBlock::puppet> top;
            Field<float, MemType::paged_device, MemBlock::puppet> bottom;
            Field<float, MemType::paged_device, MemBlock::puppet> front;
            Field<float, MemType::paged_device, MemBlock::puppet> back;
            Field<float, MemType::paged_device, MemBlock::puppet> right;
            Field<float, MemType::paged_device, MemBlock::puppet> left;
        };
        float R = 0.001;
        host_damp d_x, d_y, d_z, d_x_half, d_y_half, d_z_half;
        //
        void link(elasticGeoModel<domain::local> *_geomodel_p, SeisRecord *_seisrecord_p,
                         elasticWave *_sub_ext_wavefield_p);
        void link_geomodel_para(Field<float, MemType::paged_device> *_lambda_p, Field<float, MemType::paged_device> *_mu_p, Field<float, MemType::paged_device> *_rho_p);
        // FUNC
        void set_grid_dir();                         // 1st
        void cu_alloc_only_device_memory_wave_dir(); // 2nd
        void print_is_pml(int mpi_rank = 0);
        // damp
        virtual void cu_alloc_damp_coeff() = 0;    // 3rd
        virtual void calculate_damp_coeff() = 0;   // 4th
        virtual void cu_copy_damp_coeff_h2d() = 0; // 5th
    };
    //
    class mpicuPMLBase : public cuPMLBase
    {
    public:
        struct mpicuHalo
        {
            mpicu_halo_comm_Operator<float, MemBlock::multiple> *mpi_multi_halo_vel_p;
            mpicu_halo_comm_Operator<float, MemBlock::multiple> *mpi_multi_halo_stress_p;
            mpicu_halo_comm_Operator<float, MemBlock::puppet> halo_vx;
            mpicu_halo_comm_Operator<float, MemBlock::puppet> halo_vy;
            mpicu_halo_comm_Operator<float, MemBlock::puppet> halo_vz;
            mpicu_halo_comm_Operator<float, MemBlock::puppet> halo_sxx;
            mpicu_halo_comm_Operator<float, MemBlock::puppet> halo_syy;
            mpicu_halo_comm_Operator<float, MemBlock::puppet> halo_szz;
            mpicu_halo_comm_Operator<float, MemBlock::puppet> halo_sxy;
            mpicu_halo_comm_Operator<float, MemBlock::puppet> halo_sxz;
            mpicu_halo_comm_Operator<float, MemBlock::puppet> halo_syz;
            //
            void batch_set_pad_halo_grid(Position pos_flag, Frame &grid_pml, elasticWave *ext_wavefield_p,
                                         mpicu_halo_comm_Operator<float, MemBlock::multiple> &_mpi_multi_halo_vel,
                                         mpicu_halo_comm_Operator<float, MemBlock::multiple> &_mpi_multi_halo_stress);
            void batch_extract_halo_vel();
            void batch_extract_halo_stress();
        };
        //
        mpicuHalo top_pml_halo;
        mpicuHalo bottom_pml_halo;
        mpicuHalo front_pml_halo;
        mpicuHalo back_pml_halo;
        mpicuHalo right_pml_halo;
        mpicuHalo left_pml_halo;
        //
        Padding sub_pad;
        //
        mpicu_halo_comm_Operator<float, MemBlock::multiple> mpi_multi_halo_vel;
        mpicu_halo_comm_Operator<float, MemBlock::multiple> mpi_multi_halo_stress;
        //
        void initialize();
        void mpi_init();
        //
        void set_zero();
        //
    protected:
        void set_sub_padding_and_grid();
        void set_is_pml();
    };
}
#endif