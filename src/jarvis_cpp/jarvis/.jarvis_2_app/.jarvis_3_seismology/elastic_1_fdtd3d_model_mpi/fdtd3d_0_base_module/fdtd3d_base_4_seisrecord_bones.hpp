#pragma once
#ifndef _FDTD3D_SEISRECORD_BONES_HPP
#define _FDTD3D_SEISRECORD_BONES_HPP
#include "fdtd3d_base_3_geomodel_meat.hpp"
namespace jarvis
{
    struct SeisRecord : public elastic_fdtd3d_module_Base
    {
    public:
        elasticGeoModel<domain::local> *geomodel_p = nullptr;
        SubGeometry *sub_geometry_p = 0;
        elasticWave *ext_wavefield_p = 0;
        // seisdata
        float fm;
        float dt = -1;
        float T;
        int ntime;
        Field<float, MemType::paged_device> wave_let;
        Field<float, MemType::device> vx_seis;
        Field<float, MemType::device> vy_seis;
        Field<float, MemType::device> vz_seis;
        Field<float, MemType::device> vxp_seis;
        Field<float, MemType::device> vyp_seis;
        Field<float, MemType::device> vzp_seis;
        //
        Field<Field<float, MemType::pinned>> vx_host_vec;
        Field<Field<float, MemType::pinned>> vy_host_vec;
        Field<Field<float, MemType::pinned>> vz_host_vec;
        //
        void link(elasticGeoModel<domain::local> *_geomodel_p,
                  elasticWave *_sub_ext_wavefield_p,
                  SubGeometry *_SubGeometry_ptr)
        {
            geomodel_p = _geomodel_p;
            sub_geometry_p = _SubGeometry_ptr;
            ext_wavefield_p = _sub_ext_wavefield_p;
        }
        //
        void initialize();
        //
        void cuda_get_seis_data(int _i_time);
        //
        void alloc_host_storage();
        void copy_seis_into_host(int _i_shot);
        void copy_seis_from_host(int _i_shot);
    };
    //
    // struct PoyntingVec
    // {
    //     elasticGeoModel *geomodel_p;
    //     SubGeometry *survey_p;
    //     SeisRecord *seis_record_p;
    //     elasticWave *sub_ext_wavefield_p;
    //     cudaStream_t *cal_stream_p;
    //     Field<uvec3> pyt_vec;
    //     Field<uvec3> pyt_vec_rece;
    //     //
    //     void link_mpi_cuda_stream(cudaStream_t *_calc_stream_p)
    //     {
    //         cal_stream_p = _calc_stream_p;
    //     }
    //     //
    //     template <class GeoModelType>
    //     void link(GeoModelType *_geomodel_p,
    //                      SeisRecord *_seis_record_p,
    //                      elasticWave *_sub_ext_wavefield_p,
    //                      SubGeometry *_survey_p)
    //     {
    //         geomodel_p = _geomodel_p;
    //         survey_p = _survey_p;
    //         seis_record_p = _seis_record_p;
    //         sub_ext_wavefield_p = _sub_ext_wavefield_p;
    //     }
    //     //
    //     void initialize()
    //     {
    //         pyt_vec.alloc(MemPage::npin, geomodel_p->gridphy);
    //         pyt_vec_rece.alloc(MemPage::npin, seis_record_p->ntime, survey_p->sub_recv.frame.n_elem);
    //     }
    //     void cuda_cal_poynting_vec();
    //     void cuda_get_rece_pyt(int _i_time);
    //     void cudax_fliter_rece_pyt(Position pos_flag, Field<float,MemType::paged_device> *_seis);
    // };
}
#endif