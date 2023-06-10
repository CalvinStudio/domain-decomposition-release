#pragma once
#ifndef _FDTD3D_BLEND_ELASTIC_MODEL_MPI_HPP
#define _FDTD3D_BLEND_ELASTIC_MODEL_MPI_HPP
#include ".fdtd3d_blend_elastic_model_mpi_header_in.h"
namespace jarvis
{
    class fdtd3d_blend_elastic_model_mpi_module_Base : public elastic_fdtd3d_module_Base
    {
    public:
        elasticGeoModel<domain::global> *glb_model_p = nullptr;
        elasticGeoModel<domain::local> *sub_model_p = nullptr;
        SubGeometry *sub_survey_p = nullptr;
        elasticWave *sub_wave_p = nullptr;
        mpicuCPML *sub_pml_p = nullptr;
        SeisRecord *seis_record_p = nullptr;
    };

    template <SimulateType simulate_type>
    class fdtd3d_blend_elastic_model_mpi_module
    {
    };
    template <>
    class fdtd3d_blend_elastic_model_mpi_module<SimulateType::pure_forward> : public fdtd3d_blend_elastic_model_mpi_module_Base
    {
    public:
        mpicuPhy<SimulateType::pure_forward> *sub_phy_p = nullptr;
        void link();
        void initialize();
        void mpi_sub_forward(int _i_shot);
        void reset_for_next_shot();
    };
    template <>
    class fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_forward> : public fdtd3d_blend_elastic_model_mpi_module_Base
    {
    public:
        mpicuPhy<SimulateType::rtm_forward> *sub_phy_p = nullptr;
        elasticWaveRecon *sub_wrc_p = nullptr;

    public:
        void link(elasticGeoModel<domain::local> *_sub_model_p, SubGeometry *_sub_survey_p);
        void initialize();
        void mpi_sub_forward_for_restruct(int _i_shot);
        void reset_for_next_shot();
    };
    using fdtd3d_for_restruct_M = fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_forward>;
    template <>
    class fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_backward> : public fdtd3d_blend_elastic_model_mpi_module_Base
    {
    public:
        mpicuPhy<SimulateType::rtm_backward> *sub_phy_p = nullptr;
        elasticWaveRecon *sub_wrc_p = nullptr;

    public:
        void link(fdtd3d_for_restruct_M *_fdtd3d_for_for_restruct_p);
        void initialize();
        void reset_for_next_shot();
    };
    template <>
    class fdtd3d_blend_elastic_model_mpi_module<SimulateType::rtm_reverse> : public fdtd3d_blend_elastic_model_mpi_module_Base
    {
    public:
        mpicuPhy<SimulateType::rtm_reverse> *sub_phy_p = nullptr;
        void link(elasticGeoModel<domain::local> *_sub_model_p, SubGeometry *_sub_survey_p);
        void initialize();
        void reset_for_next_shot();
    };
}
#endif