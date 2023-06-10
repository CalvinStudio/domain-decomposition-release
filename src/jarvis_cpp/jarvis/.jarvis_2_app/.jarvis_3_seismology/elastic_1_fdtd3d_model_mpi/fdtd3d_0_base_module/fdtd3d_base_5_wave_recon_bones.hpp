#pragma once
#ifndef _FDTD3D_WRC_BONES_HPP
#define _FDTD3D_WRC_BONES_HPP
#include "fdtd3d_base_4_seisrecord_meat.hpp"
namespace jarvis
{
    class elasticWaveRecon : public elastic_fdtd3d_module_Base
    {
    public:
        //*MODULE INPUT
        elasticGeoModel<domain::local> *geomodel_p = nullptr;
        SeisRecord *seis_record_p = nullptr;
        elasticWave *ext_wavefield_p = nullptr;
        //*CUDA STREAM
        //*WAVE FIELD FOR WAVERECONSTRUCT
        cu_halo_store_Operator<float> halo_vx;
        cu_halo_store_Operator<float> halo_vy;
        cu_halo_store_Operator<float> halo_vz;
        cu_halo_store_Operator<float> halo_sxx;
        cu_halo_store_Operator<float> halo_sxy;
        cu_halo_store_Operator<float> halo_sxz;
        cu_halo_store_Operator<float> halo_syy;
        cu_halo_store_Operator<float> halo_syz;
        cu_halo_store_Operator<float> halo_szz;
        cu_halo_store_Operator<float> halo_sau;
        //
        //*FUNC
        void link(elasticGeoModel<domain::local> *_geomodel_p, SeisRecord *_seisrecord_p, elasticWave *_sub_ext_wavefield_p);
        //
        void initialize();
        //
        void extract_boundary_wave(int it);
        void copy_into_host_storage(int it);
        //
        void copy_vel_from_host_storage(int it);
        void copy_stress_from_host_storage(int it);
        void inject_boundary_vel(int it);
        void inject_boundary_stress(int it);
    };
}
#endif