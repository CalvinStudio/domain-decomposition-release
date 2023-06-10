#pragma once
#ifndef _FDTD3D_PHY_MPI_BONES_HPP
#define _FDTD3D_PHY_MPI_BONES_HPP
#include ".fdtd3d_phy_mpi_header_in.h"
namespace jarvis
{
    struct mpiHaloWave_p
    {
        float *top_vx, *top_vy, *top_vz, *top_sxz, *top_syz, *top_szz, *top_sau;
        float *bottom_vx, *bottom_vy, *bottom_vz, *bottom_sxz, *bottom_syz, *bottom_szz, *bottom_sau;
        float *front_vx, *front_vy, *front_vz, *front_sxy, *front_syy, *front_syz, *front_sau;
        float *back_vx, *back_vy, *back_vz, *back_sxy, *back_syy, *back_syz, *back_sau;
        float *right_vx, *right_vy, *right_vz, *right_sxx, *right_sxy, *right_sxz, *right_sau;
        float *left_vx, *left_vy, *left_vz, *left_sxx, *left_sxy, *left_sxz, *left_sau;
    };
    //
    class mpicuPhy_Base : public elastic_fdtd3d_module_Base
    {
    public:
        //**********MODULE INPUT**********
        elasticGeoModel<domain::local> *sub_model_p = nullptr;
        elasticWave *ext_wavefield_p = nullptr;
        SeisRecord *seis_record_p = nullptr;
        SubGeometry *sub_geometry_p = nullptr;
        //**********FOR PARA SHIFT**********
        Field<float, MemType::paged_device> *lambda_p = nullptr;
        Field<float, MemType::paged_device> *mu_p = nullptr;
        Field<float, MemType::paged_device> *rho_p = nullptr;

    protected:
        //**********DIFF COEFF**********
        Field<float, MemType::paged_device> phy_dc;
        Field<float, MemType::paged_device> phy_dc_list;
        //**********CUDA GRID**********
        Padding sub_pad;
        Frame grid_inside;
        Frame grid_top;
        Frame grid_bottom;
        Frame grid_front;
        Frame grid_back;
        Frame grid_right;
        Frame grid_left;
        //**********MPI**********
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_vx;
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_vy;
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_vz;
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_sxx;
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_sxy;
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_sxz;
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_syy;
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_syz;
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_szz;
        //***********************
        mpiHaloWave_p mpi_halo_wave;

    public:
        mpicu_halo_comm_Operator<float, MemBlock::multiple> mpi_multi_halo_vel;
        mpicu_halo_comm_Operator<float, MemBlock::multiple> mpi_multi_halo_stress;

    public:
        //**********MODULE INITIALIZE**********
        void set_base_wave_mpi();
        void set_base_wave_mpi_halo_wave();
        //
        //**********MODULE LINK**********
        void link(elasticGeoModel<domain::local> *_geomodel_p, SeisRecord *_seisrecord_p, elasticWave *_sub_ext_wavefield_p, SubGeometry *_geometry_p);
        void link_geomodel_para(Field<float, MemType::paged_device> *_lambda_p, Field<float, MemType::paged_device> *_mu_p, Field<float, MemType::paged_device> *_rho_p);
        //**********RESET**********
        void set_zero();
        void clear();
        void check_is_stable();

    protected:
        void set_grid();
        void update_outer_vel(SimulateType _simulate_type, int _start_fdorder, float _dt);
        void update_outer_stress(SimulateType _simulate_type, int _start_fdorder, float _dt);
        void update_inner_vel(SimulateType _simulate_type, float _dt);
        void update_inner_stress(SimulateType _simulate_type, float _dt);
    };

    template <SimulateType simulate_type>
    class mpicuPhy
    {
    };
    template <>
    class mpicuPhy<SimulateType::pure_forward> : public mpicuPhy_Base
    {
    public:
        void initialize();
        void set_phy_frame();
        void cuda_update_and_mpi_exchange_edge_vel();
        void cuda_update_and_mpi_exchange_edge_stress();
        void cuda_update_inside_vel();
        void cuda_update_inside_stress();
        //**********ADD SOURCE********
        void cuda_add_stress_source(int _it, int _i_shot);
    };
    template <>
    class mpicuPhy<SimulateType::rtm_forward> : public mpicuPhy_Base
    {
    public:
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_sau;

    public:
        void initialize();
        void set_rtm_wave_mpi();
        void set_rtm_wave_mpi_halo_wave();
        void set_phy_frame();
        void cuda_update_and_mpi_exchange_edge_vel();
        void cuda_update_and_mpi_exchange_edge_stress();
        void cuda_update_inside_vel();
        void cuda_update_inside_stress();
        //**********ADD SOURCE********
        void cuda_add_stress_source(int _it, int _i_shot);
    };
    template <>
    class mpicuPhy<SimulateType::rtm_backward> : public mpicuPhy_Base
    {
    public:
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_sau;

    public:
        void initialize();
        void set_rtm_wave_mpi();
        void set_rtm_wave_mpi_halo_wave();
        void set_phy_frame();
        void cuda_update_and_mpi_exchange_edge_vel();
        void cuda_update_and_mpi_exchange_edge_stress();
        void cuda_update_inside_vel();
        void cuda_update_inside_stress();
        //
        void cuda_subtract_stress_source(int _it, int _i_shot);
    };
    template <>
    class mpicuPhy<SimulateType::rtm_reverse> : public mpicuPhy_Base
    {
    public:
        mpicu_halo_comm_Operator<float, MemBlock::puppet> mpi_halo_sau;

    public:
        void initialize();
        void set_rtm_wave_mpi();
        void set_rtm_wave_mpi_halo_wave();
        void set_phy_frame();
        void cuda_update_and_mpi_exchange_edge_vel();
        void cuda_update_and_mpi_exchange_edge_stress();
        void cuda_update_inside_vel();
        void cuda_update_inside_stress();
        //
        void cuda_rtm_reverse_add_vel_source(int it, SeisRecord *raw_seis);
    };
}
#endif