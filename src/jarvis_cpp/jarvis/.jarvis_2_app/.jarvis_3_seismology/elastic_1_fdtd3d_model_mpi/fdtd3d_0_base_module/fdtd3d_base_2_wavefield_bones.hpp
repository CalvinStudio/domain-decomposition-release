#pragma once
#ifndef _FDTD3D_WAVEFIELD_BONES_HPP
#define _FDTD3D_WAVEFIELD_BONES_HPP
#include "fdtd3d_base_1_interface.hpp"
#include <unordered_map>
namespace jarvis
{
    class elasticWave : public elastic_fdtd3d_module_Base
    {
    public:
        struct orgWave
        {
            float *vx, *vy, *vz, *sxx, *sxy, *sxz, *syy, *syz, *szz;
        };
        struct rtmWave
        {
            float *vxp, *vyp, *vzp, *sau;
        };
        struct fwiWave
        {
            float *vxrt, *vyrt, *vzrt, *ux, *uy, *uz;
        };
        //*org 9
        Field<float, MemType::device> vx, vy, vz;
        Field<float, MemType::device> sxx, syy, szz;
        Field<float, MemType::device> sxy, sxz, syz;
        orgWave org_wave;
        //*rtm 4
        Field<float, MemType::device> vxp, vyp, vzp, sau;
        rtmWave rtm_wave;
        //*fwi 6
        Field<float, MemType::device> vxrt, vyrt, vzrt;
        Field<float, MemType::device> ux, uy, uz;
        fwiWave fwi_wave;
        //
        std::unordered_map<string, Field<float, MemType::device> *> wave_container;
        //
        void initialize_for_org(const Frame &_gridext)
        {
            _init_container();
            alloc_org_wavefield(_gridext);
        }
        void initialize_for_rtm(const Frame &_gridext)
        {
            _init_container();
            alloc_org_wavefield(_gridext);
            alloc_rtm_wavefield(_gridext);
        }
        //
        void alloc_org_wavefield(const Frame &grid);
        void alloc_rtm_wavefield(const Frame &grid);
        void alloc_fwi_wavefield(const Frame &grid);
        void set_zero();
        void clear();
        //
    protected:
        void _init_container();
    };
    //
    struct MemoryWave
    {
        Field<float, MemType::device> vx_x, vx_y, vx_z;
        Field<float, MemType::device> vy_x, vy_y, vy_z;
        Field<float, MemType::device> vz_x, vz_y, vz_z;
        Field<float, MemType::device> sxx_x, sxx_y, sxx_z;
        Field<float, MemType::device> syy_x, syy_y, syy_z;
        Field<float, MemType::device> szz_x, szz_y, szz_z;
        Field<float, MemType::device> sxy_x, sxy_y;
        Field<float, MemType::device> sxz_x, sxz_z;
        Field<float, MemType::device> syz_y, syz_z;
        std::unordered_map<string, Field<float, MemType::device> *> wave_container;
        void alloc_for_x_dir(const Frame &_frame);
        void alloc_for_y_dir(const Frame &_frame);
        void alloc_for_z_dir(const Frame &_frame);
        void set_zero();
        void clear();

    protected:
        void _init_container();
    };
}
#endif