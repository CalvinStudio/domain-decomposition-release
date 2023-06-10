#pragma once
#ifndef _FDTD3D_WAVEFIELD_MEAT_HPP
#define _FDTD3D_WAVEFIELD_MEAT_HPP
#include "fdtd3d_base_2_wavefield_bones.hpp"
namespace jarvis
{
    inline void elasticWave::_init_container()
    {
        wave_container["vx"] = &vx;
        wave_container["vy"] = &vy;
        wave_container["vz"] = &vz;
        wave_container["sxx"] = &sxx;
        wave_container["syy"] = &syy;
        wave_container["szz"] = &szz;
        wave_container["sxy"] = &sxy;
        wave_container["sxz"] = &sxz;
        wave_container["syz"] = &syz;
        wave_container["vxp"] = &vxp;
        wave_container["vyp"] = &vyp;
        wave_container["vzp"] = &vzp;
        wave_container["sau"] = &sau;
    }
    //
    inline void elasticWave::alloc_org_wavefield(const Frame &grid)
    {
        vx.alloc(grid);
        vy.alloc(grid);
        vz.alloc(grid);
        sxx.alloc(grid);
        sxy.alloc(grid);
        sxz.alloc(grid);
        syy.alloc(grid);
        syz.alloc(grid);
        szz.alloc(grid);
        sau.alloc(grid);
        org_wave.vx = vx.ptr();
        org_wave.vy = vy.ptr();
        org_wave.vz = vz.ptr();
        org_wave.sxx = sxx.ptr();
        org_wave.sxy = sxy.ptr();
        org_wave.sxz = sxz.ptr();
        org_wave.syy = syy.ptr();
        org_wave.syz = syz.ptr();
        org_wave.szz = szz.ptr();
        rtm_wave.sau = sau.ptr();
    }
    //
    inline void elasticWave::alloc_rtm_wavefield(const Frame &grid)
    {
        vxp.alloc(grid);
        vyp.alloc(grid);
        vzp.alloc(grid);
        rtm_wave.vxp = vxp.ptr();
        rtm_wave.vyp = vyp.ptr();
        rtm_wave.vzp = vzp.ptr();
    }
    //
    inline void elasticWave::alloc_fwi_wavefield(const Frame &grid)
    {
        vxrt.alloc(grid);
        vyrt.alloc(grid);
        vzrt.alloc(grid);
        ux.alloc(grid);
        uy.alloc(grid);
        uz.alloc(grid);
        fwi_wave.vxrt = vxrt.ptr();
        fwi_wave.vyrt = vyrt.ptr();
        fwi_wave.vzrt = vzrt.ptr();
        fwi_wave.ux = ux.ptr();
        fwi_wave.uy = uy.ptr();
        fwi_wave.uz = uz.ptr();
    }
    //
    inline void elasticWave::set_zero()
    {
        if (wave_container.size() > 0)
        {
            for (auto iter = wave_container.begin(); iter != wave_container.end(); ++iter)
            {
                if (!iter->second->Frame::is_empty())
                {
                    iter->second->set_zero();
                }
            }
        }
        else
        {
            printf("wave_container is empty!\n");
            std::abort();
        }
    }
    //
    inline void elasticWave::clear()
    {
        if (wave_container.size() > 0)
        {
            for (auto iter = wave_container.begin(); iter != wave_container.end(); ++iter)
            {
                if (!iter->second->Frame::is_empty())
                {
                    iter->second->clear();
                }
            }
        }
        else
        {
            printf("wave_container is empty!\n");
            std::abort();
        }
    }
    //
    inline void MemoryWave::_init_container()
    {
        wave_container["vx_x"] = &vx_x;
        wave_container["vy_x"] = &vy_x;
        wave_container["vz_x"] = &vz_x;
        wave_container["vx_y"] = &vx_y;
        wave_container["vy_y"] = &vy_y;
        wave_container["vz_y"] = &vz_y;
        wave_container["vx_z"] = &vx_z;
        wave_container["vy_z"] = &vy_z;
        wave_container["vz_z"] = &vz_z;
        wave_container["sxx_x"] = &sxx_x;
        wave_container["syy_x"] = &syy_x;
        wave_container["szz_x"] = &szz_x;
        wave_container["sxz_x"] = &sxz_x;
        wave_container["sxy_x"] = &sxy_x;
        wave_container["sxx_y"] = &sxx_y;
        wave_container["syy_y"] = &syy_y;
        wave_container["szz_y"] = &szz_y;
        wave_container["syz_y"] = &syz_y;
        wave_container["sxy_y"] = &sxy_y;
        wave_container["sxx_z"] = &sxx_z;
        wave_container["syy_z"] = &syy_z;
        wave_container["szz_z"] = &szz_z;
        wave_container["syz_z"] = &syz_z;
        wave_container["sxz_z"] = &sxz_z;
    }
    inline void MemoryWave::alloc_for_x_dir(const Frame &_frame)
    {
        _init_container();
        vx_x.alloc(_frame);
        vy_x.alloc(_frame);
        vz_x.alloc(_frame);
        sxx_x.alloc(_frame);
        sxy_x.alloc(_frame);
        sxz_x.alloc(_frame);
    }
    inline void MemoryWave::alloc_for_y_dir(const Frame &_frame)
    {
        _init_container();
        vx_x.alloc(_frame);
        vx_y.alloc(_frame);
        vy_x.alloc(_frame);
        vy_y.alloc(_frame);
        vz_x.alloc(_frame);
        vz_y.alloc(_frame);
        sxx_x.alloc(_frame);
        sxy_x.alloc(_frame);
        sxz_x.alloc(_frame);
        syy_y.alloc(_frame);
        sxy_y.alloc(_frame);
        syz_y.alloc(_frame);
    }
    inline void MemoryWave::alloc_for_z_dir(const Frame &_frame)
    {
        _init_container();
        for (auto iter = wave_container.begin(); iter != wave_container.end(); ++iter)
        {
            iter->second->alloc(_frame);
        }
    }
    //
    inline void MemoryWave::set_zero()
    {
        if (wave_container.size() > 0)
        {
            for (auto iter = wave_container.begin(); iter != wave_container.end(); ++iter)
            {
                if (!iter->second->Frame::is_empty())
                {
                    iter->second->set_zero();
                }
            }
        }
        else
        {
            printf("wave_container is empty!\n");
            std::abort();
        }
    }
    //
    inline void MemoryWave::clear()
    {
        if (wave_container.size() > 0)
        {
            for (auto iter = wave_container.begin(); iter != wave_container.end(); ++iter)
            {
                if (!iter->second->Frame::is_empty())
                {
                    iter->second->clear();
                }
            }
        }
        else
        {
            printf("wave_container is empty!\n");
            std::abort();
        }
    }
}
#endif