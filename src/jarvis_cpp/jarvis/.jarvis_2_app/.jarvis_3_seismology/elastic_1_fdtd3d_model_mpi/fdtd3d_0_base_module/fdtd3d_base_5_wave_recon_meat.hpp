#pragma once
#ifndef _FDTD3D_WRC_MEAT_HPP
#define _FDTD3D_WRC_MEAT_HPP
#include "fdtd3d_base_5_wave_recon_bones.hpp"
namespace jarvis
{
    inline void elasticWaveRecon::link(elasticGeoModel<domain::local> *_geomodel_p, SeisRecord *_seisrecord_p,
                                       elasticWave *_sub_ext_wavefield_p)
    {
        geomodel_p = _geomodel_p, seis_record_p = _seisrecord_p;
        ext_wavefield_p = _sub_ext_wavefield_p;
    }
    //
    inline void elasticWaveRecon::initialize()
    {
        //*Assignment of static variables
        halo_vx.jarvis_mpi_cuda_stream_p = jarvis_mpi_cuda_stream_p;
        int rec_width = geo_const::rec_width;
        halo_vx.set_margin_halo_grid(&ext_wavefield_p->vx, geomodel_p->gridphy, seis_record_p->ntime - 1, rec_width - 1, rec_width, rec_width - 1, rec_width, rec_width, rec_width - 1);
        halo_vy.set_margin_halo_grid(&ext_wavefield_p->vy, geomodel_p->gridphy, seis_record_p->ntime - 1, rec_width - 1, rec_width, rec_width, rec_width - 1, rec_width - 1, rec_width);
        halo_vz.set_margin_halo_grid(&ext_wavefield_p->vz, geomodel_p->gridphy, seis_record_p->ntime - 1, rec_width, rec_width - 1, rec_width - 1, rec_width, rec_width - 1, rec_width);
        //
        halo_sxx.set_margin_halo_grid(&ext_wavefield_p->sxx, geomodel_p->gridphy, seis_record_p->ntime - 1, 0, 0, 0, 0, rec_width - 1, rec_width);
        halo_syy.set_margin_halo_grid(&ext_wavefield_p->syy, geomodel_p->gridphy, seis_record_p->ntime - 1, 0, 0, rec_width - 1, rec_width, 0, 0);
        halo_szz.set_margin_halo_grid(&ext_wavefield_p->szz, geomodel_p->gridphy, seis_record_p->ntime - 1, rec_width - 1, rec_width, 0, 0, 0, 0);
        halo_sau.set_margin_halo_grid(&ext_wavefield_p->sau, geomodel_p->gridphy, seis_record_p->ntime - 1, rec_width - 1, rec_width, rec_width - 1, rec_width, rec_width - 1, rec_width);
        //
        halo_sxy.set_margin_halo_grid(&ext_wavefield_p->sxy, geomodel_p->gridphy, seis_record_p->ntime - 1, 0, 0, rec_width, rec_width - 1, rec_width, rec_width - 1);
        halo_sxz.set_margin_halo_grid(&ext_wavefield_p->sxz, geomodel_p->gridphy, seis_record_p->ntime - 1, rec_width, rec_width - 1, 0, 0, rec_width, rec_width - 1);
        halo_syz.set_margin_halo_grid(&ext_wavefield_p->syz, geomodel_p->gridphy, seis_record_p->ntime - 1, rec_width, rec_width - 1, rec_width, rec_width - 1, 0, 0);
        //
        halo_vx.alloc_halo();
        halo_vy.alloc_halo();
        halo_vz.alloc_halo();
        //
        halo_sxx.alloc_halo();
        halo_syy.alloc_halo();
        halo_szz.alloc_halo();
        halo_sau.alloc_halo();
        //
        halo_sxy.alloc_halo();
        halo_sxz.alloc_halo();
        halo_syz.alloc_halo();
    }
    //
    inline void elasticWaveRecon::copy_into_host_storage(int it)
    {
        if (it >= 0 && it <= seis_record_p->ntime - 2)
        {
            jarvis_mpi_cuda_stream_p->sync_cal_stream();
            //
            halo_vx.copy_halo_into_host(it);
            halo_vy.copy_halo_into_host(it);
            halo_vz.copy_halo_into_host(it);
            halo_sxx.copy_halo_into_host(it);
            halo_syy.copy_halo_into_host(it);
            halo_szz.copy_halo_into_host(it);
            halo_sau.copy_halo_into_host(it);
            //
            halo_sxy.copy_halo_into_host(it);
            halo_sxz.copy_halo_into_host(it);
            halo_syz.copy_halo_into_host(it);
        }
    }
    //
    inline void elasticWaveRecon::copy_vel_from_host_storage(int it)
    {
        if (it >= 0 && it <= seis_record_p->ntime - 2)
        {
            halo_vx.copy_halo_from_host(it);
            halo_vy.copy_halo_from_host(it);
            halo_vz.copy_halo_from_host(it);
        }
    }
    //
    inline void elasticWaveRecon::copy_stress_from_host_storage(int it)
    {
        if (it >= 0 && it <= seis_record_p->ntime - 2)
        {
            halo_sxx.copy_halo_from_host(it);
            halo_syy.copy_halo_from_host(it);
            halo_szz.copy_halo_from_host(it);
            halo_sau.copy_halo_from_host(it);
            halo_sxy.copy_halo_from_host(it);
            halo_sxz.copy_halo_from_host(it);
            halo_syz.copy_halo_from_host(it);
        }
    }
    inline void elasticWaveRecon::extract_boundary_wave(int it)
    {
        if (it >= 0 && it <= seis_record_p->ntime - 2)
        {
            jarvis_mpi_cuda_stream_p->sync_d2h_stream();
            halo_vx.extract_into_halo();
            halo_vy.extract_into_halo();
            halo_vz.extract_into_halo();
            halo_sxx.extract_into_halo();
            halo_syy.extract_into_halo();
            halo_szz.extract_into_halo();
            halo_sau.extract_into_halo();
            halo_sxy.extract_into_halo();
            halo_sxz.extract_into_halo();
            halo_syz.extract_into_halo();
        }
    }
    //
    inline void elasticWaveRecon::inject_boundary_vel(int it)
    {
        if (it >= 0 && it <= seis_record_p->ntime - 2)
        {
            jarvis_mpi_cuda_stream_p->sync_h2d_stream();
            halo_vx.inject_from_halo();
            halo_vy.inject_from_halo();
            halo_vz.inject_from_halo();
            jarvis_mpi_cuda_stream_p->sync_cal_stream();
        }
    }
    //
    inline void elasticWaveRecon::inject_boundary_stress(int it)
    {
        if (it >= 0 && it <= seis_record_p->ntime - 2)
        {
            jarvis_mpi_cuda_stream_p->sync_h2d_stream();
            halo_sxx.inject_from_halo();
            halo_syy.inject_from_halo();
            halo_szz.inject_from_halo();
            halo_sau.inject_from_halo();
            halo_sxy.inject_from_halo();
            halo_sxz.inject_from_halo();
            halo_syz.inject_from_halo();
            jarvis_mpi_cuda_stream_p->sync_cal_stream();
        }
    }
}
#endif