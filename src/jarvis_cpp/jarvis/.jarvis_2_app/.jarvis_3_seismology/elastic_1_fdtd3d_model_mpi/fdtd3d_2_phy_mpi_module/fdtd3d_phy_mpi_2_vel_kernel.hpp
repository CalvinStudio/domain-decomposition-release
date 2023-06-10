#pragma once
#ifndef _FDTD3D_PHY_MPI_VEL_KERNEL_HPP
#define _FDTD3D_PHY_MPI_VEL_KERNEL_HPP
#include "fdtd3d_phy_mpi_1_host_meat.hpp"
namespace jarvis
{
    inline __global__ void
    sub_update_vel_phy_domain_inside(SimulateType simulate_type,
                                     Frame sub_gridext, Frame grid_inside,
                                     float dt, float *phy_dc, float *rho,
                                     elasticWave::orgWave org_wave,
                                     elasticWave::rtmWave rtm_wave,
                                     elasticWave::fwiWave fwi_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_inside);
        float sxx_x_tmp = 0.f, sxy_y_tmp = 0.f, sxz_z_tmp = 0.f,
              syx_x_tmp = 0.f, syy_y_tmp = 0.f, syz_z_tmp = 0.f,
              szx_x_tmp = 0.f, szy_y_tmp = 0.f, szz_z_tmp = 0.f;
        if (idx < grid_inside.n_elem)
        {
            float dt_rho = dt / rho[ext_idx];
            float a = dt_rho / sub_gridext.d_rows;
            float b = dt_rho / sub_gridext.d_cols;
            float c = dt_rho / sub_gridext.d_slices;
            int ir;
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                szz_z_tmp += phy_dc[ir] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
            }
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
            }
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                sxz_z_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
            }
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                syz_z_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
            }

            org_wave.vx[ext_idx] += a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp;
            org_wave.vy[ext_idx] += a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp;
            org_wave.vz[ext_idx] += a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp;

            if (simulate_type == SimulateType::rtm_backward || simulate_type == SimulateType::rtm_reverse)
            {
                float sau_x_tmp = 0.f, sau_y_tmp = 0.f, sau_z_tmp = 0.f;
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                    sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                    sau_z_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                }
                rtm_wave.vxp[ext_idx] += a * sau_x_tmp;
                rtm_wave.vyp[ext_idx] += b * sau_y_tmp;
                rtm_wave.vzp[ext_idx] += c * sau_z_tmp;
            }
            else if (simulate_type == SimulateType::fwi)
            {
                fwi_wave.ux[ext_idx] += dt * org_wave.vx[ext_idx];
                fwi_wave.uy[ext_idx] += dt * org_wave.vy[ext_idx];
                fwi_wave.uz[ext_idx] += dt * org_wave.vz[ext_idx];
                fwi_wave.vxrt[ext_idx] = dt * (a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp);
                fwi_wave.vyrt[ext_idx] = dt * (a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp);
                fwi_wave.vzrt[ext_idx] = dt * (a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp);
            }
        }
    }
    //
    inline __global__ void
    sub_update_vel_phy_domain_edge_top(SimulateType simulate_type,
                                       Frame sub_gridext, Frame sub_gridphy, Frame grid_top,
                                       Padding sub_pad,
                                       isPosition is_pos,
                                       int start_fdorder,
                                       float dt, float *phy_dc, float *phy_dc_list, float *rho,
                                       //
                                       elasticWave::orgWave org_wave,
                                       elasticWave::rtmWave rtm_wave,
                                       elasticWave::fwiWave fwi_wave,
                                       mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_top);
        int phy_n_cols = sub_gridphy.n_cols;
        float sxx_x_tmp = 0.f, sxy_y_tmp = 0.f, sxz_z_tmp = 0.f,
              syx_x_tmp = 0.f, syy_y_tmp = 0.f, syz_z_tmp = 0.f,
              szx_x_tmp = 0.f, szy_y_tmp = 0.f, szz_z_tmp = 0.f;
        if (idx < grid_top.n_elem)
        {
            float dt_rho = dt / rho[ext_idx];
            float a = dt_rho / sub_gridext.d_rows;
            float b = dt_rho / sub_gridext.d_cols;
            float c = dt_rho / sub_gridext.d_slices;
            int ir;
            if (is_pos.is_top)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir <= k)
                    {
                        szz_z_tmp += phy_dc[ir] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
                    }
                    else
                    {
                        szz_z_tmp += phy_dc[ir] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - halo_wave.top_szz[i + j * n_rows + (k - ir - 1 + geo_const::phy_fdorder_half) * n_elem_slice]);
                    }
                    if (ir < k)
                    {
                        sxz_z_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                        syz_z_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                    }
                    else
                    {
                        sxz_z_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - halo_wave.top_sxz[i + j * n_rows + (k - ir - 1 + geo_const::phy_fdorder_half) * n_elem_slice]);
                        syz_z_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - halo_wave.top_syz[i + j * n_rows + (k - ir - 1 + geo_const::phy_fdorder_half) * n_elem_slice]);
                    }
                }
            }
            else if (!is_pos.is_top)
            {
                if (k <= start_fdorder - 1)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        sxz_z_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                        syz_z_tmp += phy_dc_list[ir] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                        szz_z_tmp += phy_dc_list[ir] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
                    }
                }
                else if (k >= start_fdorder && k < geo_const::phy_fdorder_half)
                {
                    for (ir = 0; ir <= k; ++ir)
                    {
                        sxz_z_tmp += phy_dc_list[ir + (k - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                        syz_z_tmp += phy_dc_list[ir + (k - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                        szz_z_tmp += phy_dc_list[ir + (k - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
                    }
                }
            }
            //*******X
            if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                    syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                    szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                }
            }
            else if (i < sub_pad.pad_right)
            {
                if (is_pos.is_right)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir <= i)
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        }
                        else
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - halo_wave.right_sxx[i - ir - 1 + geo_const::phy_fdorder_half + j * (geo_const::phy_fdorder_half - 1) + k * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                        }
                        if (ir < i)
                        {
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                        else
                        {
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - halo_wave.right_sxy[i - ir - 1 + geo_const::phy_fdorder_half + j * geo_const::phy_fdorder_half + k * geo_const::phy_fdorder_half * phy_n_cols]);
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - halo_wave.right_sxz[i - ir - 1 + geo_const::phy_fdorder_half + j * geo_const::phy_fdorder_half + k * geo_const::phy_fdorder_half * phy_n_cols]);
                        }
                    }
                }
                else if (!is_pos.is_right)
                {
                    if (i <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                    else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= i; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
            }
            else if (i > n_rows - sub_pad.pad_left - 1 && i < n_rows)
            {
                if (is_pos.is_left)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < n_rows - 1 - i)
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        }
                        else
                        {
                            sxx_x_tmp += phy_dc[ir] * (halo_wave.left_sxx[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + k * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.sxx[ext_idx - ir]);
                        }
                        if (ir <= n_rows - 1 - i)
                        {
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                        }
                        else
                        {
                            syx_x_tmp += phy_dc[ir] * (halo_wave.left_sxy[i + ir - n_rows + j * (geo_const::phy_fdorder_half - 1) + k * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc[ir] * (halo_wave.left_sxz[i + ir - n_rows + j * (geo_const::phy_fdorder_half - 1) + k * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
                else if (!is_pos.is_left)
                {
                    int _iw = i - n_rows + geo_const::phy_fdorder_half;
                    if (_iw >= geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                    else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
            }
            //********Y
            //
            if (j >= sub_pad.pad_front && j <= phy_n_cols - sub_pad.pad_back - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                    sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                    szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                }
            }
            else if (j < sub_pad.pad_front)
            {
                if (is_pos.is_front)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir <= j)
                        {
                            syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                        else
                        {
                            syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_syy[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + k * n_rows * (geo_const::phy_fdorder_half - 1)]);
                        }
                        if (ir < j)
                        {
                            szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                            sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                        else
                        {
                            sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - halo_wave.front_sxy[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + k * n_rows * geo_const::phy_fdorder_half]);
                            szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - halo_wave.front_syz[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + k * n_rows * geo_const::phy_fdorder_half]);
                        }
                    }
                }
                else if (!is_pos.is_front)
                {
                    if (j <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxy_y_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc_list[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                            syy_y_tmp += phy_dc_list[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                    }
                    else if (j >= start_fdorder && j < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= j; ++ir)
                        {
                            sxy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                            syy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                    }
                }
            }
            else if (j > phy_n_cols - sub_pad.pad_back - 1 && j < phy_n_cols)
            {
                if (is_pos.is_back)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < phy_n_cols - 1 - j)
                        {
                            syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                        else
                        {
                            syy_y_tmp += phy_dc[ir] * (halo_wave.back_syy[i + (j + ir + 1 - phy_n_cols) * n_rows + k * n_rows * geo_const::phy_fdorder_half] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                        if (ir <= phy_n_cols - 1 - j)
                        {
                            sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                        else
                        {
                            sxy_y_tmp += phy_dc[ir] * (halo_wave.back_sxy[i + (j + ir - phy_n_cols) * n_rows + k * n_rows * (geo_const::phy_fdorder_half - 1)] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc[ir] * (halo_wave.back_syz[i + (j + ir - phy_n_cols) * n_rows + k * n_rows * (geo_const::phy_fdorder_half - 1)] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                    }
                }
                else if (!is_pos.is_back)
                {
                    int _jw = j - phy_n_cols + geo_const::phy_fdorder_half;
                    //
                    if (_jw >= geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxy_y_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc_list[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                            syy_y_tmp += phy_dc_list[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                    }
                    else if (_jw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _jw; ++ir)
                        {
                            sxy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                            syy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                    }
                }
            }
            //********************
            org_wave.vx[ext_idx] += a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp;
            org_wave.vy[ext_idx] += a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp;
            org_wave.vz[ext_idx] += a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp;
            //
            if (simulate_type == SimulateType::rtm_backward || simulate_type == SimulateType::rtm_reverse)
            {
                float sau_x_tmp = 0.f, sau_y_tmp = 0.f, sau_z_tmp = 0.f;
                if (is_pos.is_top)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir <= k)
                        {
                            sau_z_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                        }
                        else
                        {
                            sau_z_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - halo_wave.top_sau[i + j * n_rows + (k - ir - 1 + geo_const::phy_fdorder_half) * n_elem_slice]);
                        }
                    }
                }
                else if (!is_pos.is_top)
                {
                    if (k <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sau_z_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                        }
                    }
                    else if (k >= start_fdorder && k < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= k; ++ir)
                        {
                            sau_z_tmp += phy_dc_list[ir + (k - start_fdorder + 1) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                        }
                    }
                }
                //*******X*************
                if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                    }
                }
                else if (i < sub_pad.pad_right)
                {
                    if (is_pos.is_right)
                    {
#pragma unroll
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir <= i)
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                            else
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - halo_wave.right_sau[i - ir - 1 + geo_const::phy_fdorder_half + j * (geo_const::phy_fdorder_half - 1) + k * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                            }
                        }
                    }
                    else if (!is_pos.is_right)
                    {
                        if (i <= start_fdorder - 1)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                        else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                        {
                            for (ir = 0; ir <= i; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                }
                else if (i > n_rows - sub_pad.pad_left - 1 && i < n_rows)
                {
                    if (is_pos.is_left)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir < n_rows - 1 - i)
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                            else
                            {
                                sau_x_tmp += phy_dc[ir] * (halo_wave.left_sau[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + k * geo_const::phy_fdorder_half * phy_n_cols] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                    else if (!is_pos.is_left)
                    {
                        int _iw = i - n_rows + geo_const::phy_fdorder_half;
                        if (_iw >= geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                        else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                }
                //********Y
                //
                if (j >= sub_pad.pad_front && j <= phy_n_cols - sub_pad.pad_back - 1)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                    }
                }
                else if (j < sub_pad.pad_front)
                {
                    if (is_pos.is_front)
                    {
#pragma unroll
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir <= j)
                            {
                                sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                            else
                            {
                                sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_sau[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + k * n_rows * (geo_const::phy_fdorder_half - 1)]);
                            }
                        }
                    }
                    else if (!is_pos.is_front)
                    {
                        if (j <= start_fdorder - 1)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_y_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                        }
                        else if (j >= start_fdorder && j < geo_const::phy_fdorder_half)
                        {
                            for (ir = 0; ir <= j; ++ir)
                            {
                                sau_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                        }
                    }
                }
                else if (j > phy_n_cols - sub_pad.pad_back - 1 && j < phy_n_cols)
                {
                    if (is_pos.is_back)
                    {
#pragma unroll
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir < phy_n_cols - 1 - j)
                            {
                                sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                            else
                            {
                                sau_y_tmp += phy_dc[ir] * (halo_wave.back_sau[i + (j + ir + 1 - phy_n_cols) * n_rows + k * n_rows * geo_const::phy_fdorder_half] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                        }
                    }
                    else if (!is_pos.is_back)
                    {
                        int _jw = j - phy_n_cols + geo_const::phy_fdorder_half;
                        //
                        if (_jw >= geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_y_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                        }
                        else if (_jw < geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < geo_const::phy_fdorder_half - _jw; ++ir)
                            {
                                sau_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                        }
                    }
                }
                rtm_wave.vxp[ext_idx] += a * sau_x_tmp;
                rtm_wave.vyp[ext_idx] += b * sau_y_tmp;
                rtm_wave.vzp[ext_idx] += c * sau_z_tmp;
            }
            if (simulate_type == SimulateType::fwi)
            {
                fwi_wave.ux[ext_idx] += dt * org_wave.vx[ext_idx];
                fwi_wave.uy[ext_idx] += dt * org_wave.vy[ext_idx];
                fwi_wave.uz[ext_idx] += dt * org_wave.vz[ext_idx];
                fwi_wave.vxrt[ext_idx] = dt * (a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp);
                fwi_wave.vyrt[ext_idx] = dt * (a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp);
                fwi_wave.vzrt[ext_idx] = dt * (a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp);
            }
        }
    }
    //
    inline __global__ void
    sub_update_vel_phy_domain_edge_bottom(SimulateType simulate_type,
                                          Frame sub_gridext, Frame sub_gridphy, Frame grid_bottom,
                                          Padding sub_pad,
                                          isPosition is_pos,
                                          int start_fdorder,
                                          float dt, float *phy_dc, float *phy_dc_list, float *rho,
                                          //
                                          elasticWave::orgWave org_wave,
                                          elasticWave::rtmWave rtm_wave,
                                          elasticWave::fwiWave fwi_wave,
                                          mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_bottom);
        int phy_n_cols = sub_gridphy.n_cols;
        int phy_n_slices = sub_gridphy.n_slices;
        float sxx_x_tmp = 0.f, sxy_y_tmp = 0.f, sxz_z_tmp = 0.f,
              syx_x_tmp = 0.f, syy_y_tmp = 0.f, syz_z_tmp = 0.f,
              szx_x_tmp = 0.f, szy_y_tmp = 0.f, szz_z_tmp = 0.f;
        if (idx < grid_bottom.n_elem)
        {
            float dt_rho = dt / rho[ext_idx];
            float a = dt_rho / sub_gridext.d_rows;
            float b = dt_rho / sub_gridext.d_cols;
            float c = dt_rho / sub_gridext.d_slices;
            int ir;
            if (is_pos.is_bottom)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir < geo_const::phy_fdorder_half - 1 - k)
                    {
                        szz_z_tmp += phy_dc[ir] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
                    }
                    else
                    {
                        szz_z_tmp += phy_dc[ir] * (halo_wave.bottom_szz[i + j * n_rows + (k + ir + 1 - geo_const::phy_fdorder_half) * n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
                    }
                    if (ir <= geo_const::phy_fdorder_half - 1 - k)
                    {
                        sxz_z_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                        syz_z_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                    }
                    else
                    {
                        sxz_z_tmp += phy_dc[ir] * (halo_wave.bottom_sxz[i + j * n_rows + (k + ir - geo_const::phy_fdorder_half) * n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                        syz_z_tmp += phy_dc[ir] * (halo_wave.bottom_syz[i + j * n_rows + (k + ir - geo_const::phy_fdorder_half) * n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                    }
                }
            }
            else if (!is_pos.is_bottom)
            {
                if (k >= geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        sxz_z_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                        syz_z_tmp += phy_dc_list[ir] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                        szz_z_tmp += phy_dc_list[ir] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
                    }
                }
                else if (k < geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < geo_const::phy_fdorder_half - k; ++ir)
                    {
                        sxz_z_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - k - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                        syz_z_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - k - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                        szz_z_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - k - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
                    }
                }
            }
            //*********X
            if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                    syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                    szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                }
            }
            else if (i < sub_pad.pad_right)
            {
                if (is_pos.is_right)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir <= i)
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        }
                        else
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - halo_wave.right_sxx[i - ir - 1 + geo_const::phy_fdorder_half + j * (geo_const::phy_fdorder_half - 1) + (k + phy_n_slices - sub_pad.pad_bottom) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                        }
                        if (ir < i)
                        {
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                        else
                        {
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - halo_wave.right_sxy[i - ir - 1 + geo_const::phy_fdorder_half + j * geo_const::phy_fdorder_half + (k + phy_n_slices - sub_pad.pad_bottom) * geo_const::phy_fdorder_half * phy_n_cols]);
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - halo_wave.right_sxz[i - ir - 1 + geo_const::phy_fdorder_half + j * geo_const::phy_fdorder_half + (k + phy_n_slices - sub_pad.pad_bottom) * geo_const::phy_fdorder_half * phy_n_cols]);
                        }
                    }
                }
                else if (!is_pos.is_right)
                {
                    if (i <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                    else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= i; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
            }
            else if (i > n_rows - sub_pad.pad_left - 1 && i < n_rows)
            {
                if (is_pos.is_left)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < n_rows - 1 - i)
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        }
                        else
                        {
                            sxx_x_tmp += phy_dc[ir] * (halo_wave.left_sxx[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + (k + phy_n_slices - sub_pad.pad_bottom) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.sxx[ext_idx - ir]);
                        }
                        if (ir <= n_rows - 1 - i)
                        {
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                        else
                        {
                            syx_x_tmp += phy_dc[ir] * (halo_wave.left_sxy[i + ir - n_rows + j * (geo_const::phy_fdorder_half - 1) + (k + phy_n_slices - sub_pad.pad_bottom) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc[ir] * (halo_wave.left_sxz[i + ir - n_rows + j * (geo_const::phy_fdorder_half - 1) + (k + phy_n_slices - sub_pad.pad_bottom) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
                else if (!is_pos.is_left)
                {
                    int _iw = i - n_rows + geo_const::phy_fdorder_half;
                    if (_iw >= geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                    else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
            }
            //*********Y
            //
            if (j >= sub_pad.pad_front && j <= phy_n_cols - sub_pad.pad_back - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                    sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                    szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                }
            }
            else if (j < sub_pad.pad_front)
            {
                if (is_pos.is_front)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir <= j)
                        {
                            syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                        else
                        {
                            syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_syy[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * (geo_const::phy_fdorder_half - 1)]);
                        }
                        if (ir < j)
                        {
                            sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                        else
                        {
                            sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - halo_wave.front_sxy[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * geo_const::phy_fdorder_half]);
                            szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - halo_wave.front_syz[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * geo_const::phy_fdorder_half]);
                        }
                    }
                }
                else if (!is_pos.is_front)
                {
                    if (j <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxy_y_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc_list[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                            syy_y_tmp += phy_dc_list[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                    }
                    else if (j >= start_fdorder && j < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= j; ++ir)
                        {
                            sxy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                            syy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                    }
                }
            }
            else if (j > phy_n_cols - sub_pad.pad_back - 1 && j < phy_n_cols)
            {
                if (is_pos.is_back)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < phy_n_cols - 1 - j)
                        {
                            syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                        else
                        {
                            syy_y_tmp += phy_dc[ir] * (halo_wave.back_syy[i + (j + ir + 1 - phy_n_cols) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * geo_const::phy_fdorder_half] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                        if (ir <= phy_n_cols - 1 - j)
                        {
                            sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                        else
                        {
                            sxy_y_tmp += phy_dc[ir] * (halo_wave.back_sxy[i + (j + ir - phy_n_cols) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * (geo_const::phy_fdorder_half - 1)] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc[ir] * (halo_wave.back_syz[i + (j + ir - phy_n_cols) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * (geo_const::phy_fdorder_half - 1)] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                    }
                }
                else if (!is_pos.is_back)
                {
                    int _jw = j - phy_n_cols + geo_const::phy_fdorder_half;
                    //
                    if (_jw >= geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxy_y_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc_list[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                            syy_y_tmp += phy_dc_list[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                    }
                    else if (_jw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _jw; ++ir)
                        {
                            sxy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                            szy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                            syy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        }
                    }
                }
            }
            //*********
            org_wave.vx[ext_idx] += a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp;
            org_wave.vy[ext_idx] += a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp;
            org_wave.vz[ext_idx] += a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp;
            //
            if (simulate_type == SimulateType::rtm_backward || simulate_type == SimulateType::rtm_reverse)
            {
                float sau_x_tmp = 0.f, sau_y_tmp = 0.f, sau_z_tmp = 0.f;
                if (is_pos.is_bottom)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < geo_const::phy_fdorder_half - 1 - k)
                        {
                            sau_z_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                        }
                        else
                        {
                            sau_z_tmp += phy_dc[ir] * (halo_wave.bottom_sau[i + j * n_rows + (k + ir + 1 - geo_const::phy_fdorder_half) * n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                        }
                    }
                }
                else if (!is_pos.is_bottom)
                {
                    if (k >= geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sau_z_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                        }
                    }
                    else if (k < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - k; ++ir)
                        {
                            sau_z_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - k - start_fdorder) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                        }
                    }
                }
                //*********X
                if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                    }
                }
                else if (i < sub_pad.pad_right)
                {
                    if (is_pos.is_right)
                    {
#pragma unroll
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir <= i)
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                            else
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - halo_wave.right_sau[i - ir - 1 + geo_const::phy_fdorder_half + j * (geo_const::phy_fdorder_half - 1) + (k + phy_n_slices - sub_pad.pad_bottom) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                            }
                        }
                    }
                    else if (!is_pos.is_right)
                    {
                        if (i <= start_fdorder - 1)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                        else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                        {
                            for (ir = 0; ir <= i; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                }
                else if (i > n_rows - sub_pad.pad_left - 1 && i < n_rows)
                {
                    if (is_pos.is_left)
                    {
#pragma unroll
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir < n_rows - 1 - i)
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                            else
                            {
                                sau_x_tmp += phy_dc[ir] * (halo_wave.left_sau[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + (k + phy_n_slices - sub_pad.pad_bottom) * geo_const::phy_fdorder_half * phy_n_cols] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                    else if (!is_pos.is_left)
                    {
                        int _iw = i - n_rows + geo_const::phy_fdorder_half;
                        if (_iw >= geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                        else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                }
                //*********Y
                //
                if (j >= sub_pad.pad_front && j <= phy_n_cols - sub_pad.pad_back - 1)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                    }
                }
                else if (j < sub_pad.pad_front)
                {
                    if (is_pos.is_front)
                    {
#pragma unroll
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir <= j)
                            {
                                sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                            else
                            {
                                sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_sau[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * (geo_const::phy_fdorder_half - 1)]);
                            }
                        }
                    }
                    else if (!is_pos.is_front)
                    {
                        if (j <= start_fdorder - 1)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_y_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                        }
                        else if (j >= start_fdorder && j < geo_const::phy_fdorder_half)
                        {
                            for (ir = 0; ir <= j; ++ir)
                            {
                                sau_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                        }
                    }
                }
                else if (j > phy_n_cols - sub_pad.pad_back - 1 && j < phy_n_cols)
                {
                    if (is_pos.is_back)
                    {
#pragma unroll
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir < phy_n_cols - 1 - j)
                            {
                                sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                            else
                            {
                                sau_y_tmp += phy_dc[ir] * (halo_wave.back_sau[i + (j + ir + 1 - phy_n_cols) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * geo_const::phy_fdorder_half] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                        }
                    }
                    else if (!is_pos.is_back)
                    {
                        int _jw = j - phy_n_cols + geo_const::phy_fdorder_half;
                        //
                        if (_jw >= geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_y_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                        }
                        else if (_jw < geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < geo_const::phy_fdorder_half - _jw; ++ir)
                            {
                                sau_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                            }
                        }
                    }
                }

                rtm_wave.vxp[ext_idx] += a * sau_x_tmp;
                rtm_wave.vyp[ext_idx] += b * sau_y_tmp;
                rtm_wave.vzp[ext_idx] += c * sau_z_tmp;
            }
            if (simulate_type == SimulateType::fwi)
            {
                fwi_wave.ux[ext_idx] += dt * org_wave.vx[ext_idx];
                fwi_wave.uy[ext_idx] += dt * org_wave.vy[ext_idx];
                fwi_wave.uz[ext_idx] += dt * org_wave.vz[ext_idx];
                fwi_wave.vxrt[ext_idx] = dt * (a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp);
                fwi_wave.vyrt[ext_idx] = dt * (a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp);
                fwi_wave.vzrt[ext_idx] = dt * (a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp);
            }
        }
    }
    //
    inline __global__ void
    sub_update_vel_phy_domain_edge_front(SimulateType simulate_type,
                                         Frame sub_gridext, Frame sub_gridphy, Frame grid_front,
                                         Padding sub_pad,
                                         isPosition is_pos,
                                         int start_fdorder,
                                         float dt, float *phy_dc, float *phy_dc_list, float *rho,
                                         //
                                         elasticWave::orgWave org_wave,
                                         elasticWave::rtmWave rtm_wave,
                                         elasticWave::fwiWave fwi_wave,
                                         mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_front);
        int phy_n_cols = sub_gridphy.n_cols;
        float sxx_x_tmp = 0.f, sxy_y_tmp = 0.f, sxz_z_tmp = 0.f,
              syx_x_tmp = 0.f, syy_y_tmp = 0.f, syz_z_tmp = 0.f,
              szx_x_tmp = 0.f, szy_y_tmp = 0.f, szz_z_tmp = 0.f;
        if (idx < grid_front.n_elem)
        {
            float dt_rho = dt / rho[ext_idx];
            float a = dt_rho / sub_gridext.d_rows;
            float b = dt_rho / sub_gridext.d_cols;
            float c = dt_rho / sub_gridext.d_slices;
            int ir;
            //********************
            if (is_pos.is_front)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir <= j)
                    {
                        syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                    }
                    else
                    {
                        syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_syy[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * (geo_const::phy_fdorder_half - 1)]);
                    }
                    if (ir < j)
                    {
                        sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                        szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                    }
                    else
                    {
                        sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - halo_wave.front_sxy[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * geo_const::phy_fdorder_half]);
                        szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - halo_wave.front_syz[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * geo_const::phy_fdorder_half]);
                    }
                }
            }
            else if (!is_pos.is_front)
            {
                if (j <= start_fdorder - 1)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        sxy_y_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                        szy_y_tmp += phy_dc_list[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                        syy_y_tmp += phy_dc_list[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                    }
                }
                else if (j >= start_fdorder && j < geo_const::phy_fdorder_half)
                {
                    for (ir = 0; ir <= j; ++ir)
                    {
                        sxy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                        szy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                        syy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                    }
                }
            }
            //************Z
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                szz_z_tmp += phy_dc[ir] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
                sxz_z_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                syz_z_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
            }
            //***********X
            //
            if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                    syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                    szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                }
            }
            else if (i < sub_pad.pad_right)
            {
                if (is_pos.is_right)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir <= i)
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        }
                        else
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - halo_wave.right_sxx[i - ir - 1 + geo_const::phy_fdorder_half + j * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                        }
                        if (ir < i)
                        {
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                        else
                        {
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - halo_wave.right_sxy[i - ir - 1 + geo_const::phy_fdorder_half + j * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols]);
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - halo_wave.right_sxz[i - ir - 1 + geo_const::phy_fdorder_half + j * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols]);
                        }
                    }
                }
                else if (!is_pos.is_right)
                {
                    if (i <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                    else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= i; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
            }
            else if (i > n_rows - sub_pad.pad_left - 1 && i < n_rows)
            {
                if (is_pos.is_left)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < n_rows - 1 - i)
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        }
                        else
                        {
                            sxx_x_tmp += phy_dc[ir] * (halo_wave.left_sxx[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.sxx[ext_idx - ir]);
                        }
                        if (ir <= n_rows - 1 - i)
                        {
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                        else
                        {
                            syx_x_tmp += phy_dc[ir] * (halo_wave.left_sxy[i + ir - n_rows + j * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc[ir] * (halo_wave.left_sxz[i + ir - n_rows + j * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
                else if (!is_pos.is_left)
                {
                    int _iw = i - n_rows + geo_const::phy_fdorder_half;
                    //
                    if (_iw >= geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                    else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
            }
            //********************
            org_wave.vx[ext_idx] += a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp;
            org_wave.vy[ext_idx] += a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp;
            org_wave.vz[ext_idx] += a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp;
            //
            if (simulate_type == SimulateType::rtm_backward || simulate_type == SimulateType::rtm_reverse)
            {
                float sau_x_tmp = 0.f, sau_y_tmp = 0.f, sau_z_tmp = 0.f;
                if (is_pos.is_front)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir <= j)
                        {
                            sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                        }
                        else
                        {
                            sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_sau[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * (geo_const::phy_fdorder_half - 1)]);
                        }
                    }
                }
                else if (!is_pos.is_front)
                {
                    if (j <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sau_y_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                        }
                    }
                    else if (j >= start_fdorder && j < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= j; ++ir)
                        {
                            sau_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                        }
                    }
                }
                //************Z
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    sau_z_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                }
                //***********X
                //
                if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                    }
                }
                else if (i < sub_pad.pad_right)
                {
                    if (is_pos.is_right)
                    {
#pragma unroll
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir <= i)
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                            else
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - halo_wave.right_sau[i - ir - 1 + geo_const::phy_fdorder_half + j * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                            }
                        }
                    }
                    else if (!is_pos.is_right)
                    {
                        if (i <= start_fdorder - 1)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                        else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                        {
                            for (ir = 0; ir <= i; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                }
                else if (i > n_rows - sub_pad.pad_left - 1 && i < n_rows)
                {
                    if (is_pos.is_left)
                    {
#pragma unroll
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir < n_rows - 1 - i)
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                            else
                            {
                                sau_x_tmp += phy_dc[ir] * (halo_wave.left_sau[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                    else if (!is_pos.is_left)
                    {
                        int _iw = i - n_rows + geo_const::phy_fdorder_half;
                        //
                        if (_iw >= geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                        else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                }
                rtm_wave.vxp[ext_idx] += a * sau_x_tmp;
                rtm_wave.vyp[ext_idx] += b * sau_y_tmp;
                rtm_wave.vzp[ext_idx] += c * sau_z_tmp;
            }
            if (simulate_type == SimulateType::fwi)
            {
                fwi_wave.ux[ext_idx] += dt * org_wave.vx[ext_idx];
                fwi_wave.uy[ext_idx] += dt * org_wave.vy[ext_idx];
                fwi_wave.uz[ext_idx] += dt * org_wave.vz[ext_idx];
                fwi_wave.vxrt[ext_idx] = dt * (a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp);
                fwi_wave.vyrt[ext_idx] = dt * (a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp);
                fwi_wave.vzrt[ext_idx] = dt * (a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp);
            }
        }
    }
    //
    inline __global__ void
    sub_update_vel_phy_domain_edge_back(SimulateType simulate_type,
                                        Frame sub_gridext, Frame sub_gridphy, Frame grid_back,
                                        Padding sub_pad,
                                        isPosition is_pos,
                                        int start_fdorder,
                                        float dt, float *phy_dc, float *phy_dc_list, float *rho,
                                        //
                                        elasticWave::orgWave org_wave,
                                        elasticWave::rtmWave rtm_wave,
                                        elasticWave::fwiWave fwi_wave,
                                        mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_back);
        int phy_n_cols = sub_gridphy.n_cols;
        float sxx_x_tmp = 0.f, sxy_y_tmp = 0.f, sxz_z_tmp = 0.f,
              syx_x_tmp = 0.f, syy_y_tmp = 0.f, syz_z_tmp = 0.f,
              szx_x_tmp = 0.f, szy_y_tmp = 0.f, szz_z_tmp = 0.f;
        if (idx < grid_back.n_elem)
        {
            float dt_rho = dt / rho[ext_idx];
            float a = dt_rho / sub_gridext.d_rows;
            float b = dt_rho / sub_gridext.d_cols;
            float c = dt_rho / sub_gridext.d_slices;
            int ir;
            //
            if (is_pos.is_back)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir < geo_const::phy_fdorder_half - 1 - j)
                    {
                        syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                    }
                    else
                    {
                        syy_y_tmp += phy_dc[ir] * (halo_wave.back_syy[i + (j + ir + 1 - geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * geo_const::phy_fdorder_half] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                    }
                    if (ir <= geo_const::phy_fdorder_half - 1 - j)
                    {
                        sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                        szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                    }
                    else
                    {
                        sxy_y_tmp += phy_dc[ir] * (halo_wave.back_sxy[i + (j + ir - geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * (geo_const::phy_fdorder_half - 1)] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                        szy_y_tmp += phy_dc[ir] * (halo_wave.back_syz[i + (j + ir - geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * (geo_const::phy_fdorder_half - 1)] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                    }
                }
            }
            else if (!is_pos.is_back)
            {
                if (j >= geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        syy_y_tmp += phy_dc_list[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        sxy_y_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                        szy_y_tmp += phy_dc_list[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                    }
                }
                else if (j < geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < geo_const::phy_fdorder_half - j; ++ir)
                    {
                        syy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - j - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                        sxy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - j - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                        szy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - j - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
                    }
                }
            }
            //********************
            //
            if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                    syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                    szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                }
            }
            else if (i < sub_pad.pad_right)
            {
                if (is_pos.is_right)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir <= i)
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        }
                        else
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - halo_wave.right_sxx[i - ir - 1 + geo_const::phy_fdorder_half + (j + phy_n_cols - sub_pad.pad_back) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                        }
                        if (ir < i)
                        {
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                        else
                        {
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - halo_wave.right_sxy[i - ir - 1 + geo_const::phy_fdorder_half + (j + phy_n_cols - sub_pad.pad_back) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols]);
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - halo_wave.right_sxz[i - ir - 1 + geo_const::phy_fdorder_half + (j + phy_n_cols - sub_pad.pad_back) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols]);
                        }
                    }
                }
                else if (!is_pos.is_right)
                {
                    if (i <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                    else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= i; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
            }
            else if (i > n_rows - sub_pad.pad_left - 1 && i < n_rows)
            {
                if (is_pos.is_left)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < n_rows - 1 - i)
                        {
                            sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        }
                        else
                        {
                            sxx_x_tmp += phy_dc[ir] * (halo_wave.left_sxx[i + ir + 1 - n_rows + (j + phy_n_cols - sub_pad.pad_back) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.sxx[ext_idx - ir]);
                        }
                        if (ir <= n_rows - 1 - i)
                        {
                            syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                        else
                        {
                            syx_x_tmp += phy_dc[ir] * (halo_wave.left_sxy[i + ir - n_rows + (j + phy_n_cols - sub_pad.pad_back) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc[ir] * (halo_wave.left_sxz[i + ir - n_rows + (j + phy_n_cols - sub_pad.pad_back) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
                else if (!is_pos.is_left)
                {
                    int _iw = i - n_rows + geo_const::phy_fdorder_half;
                    //
                    if (_iw >= geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                    else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                        {
                            sxx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                            syx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                            szx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                        }
                    }
                }
            }
            //********************
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                sxz_z_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                syz_z_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                szz_z_tmp += phy_dc[ir] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
            }
            //********************
            org_wave.vx[ext_idx] += a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp;
            org_wave.vy[ext_idx] += a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp;
            org_wave.vz[ext_idx] += a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp;
            //
            if (simulate_type == SimulateType::rtm_backward || simulate_type == SimulateType::rtm_reverse)
            {
                float sau_x_tmp = 0.f, sau_y_tmp = 0.f, sau_z_tmp = 0.f;
                if (is_pos.is_back)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < geo_const::phy_fdorder_half - 1 - j)
                        {
                            sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                        }
                        else
                        {
                            sau_y_tmp += phy_dc[ir] * (halo_wave.back_sau[i + (j + ir + 1 - geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * geo_const::phy_fdorder_half] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                        }
                    }
                }
                else if (!is_pos.is_back)
                {
                    if (j >= geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sau_y_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                        }
                    }
                    else if (j < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - j; ++ir)
                        {
                            sau_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - j - start_fdorder) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                        }
                    }
                }
                //********************
                //
                if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                    }
                }
                else if (i < sub_pad.pad_right)
                {
                    if (is_pos.is_right)
                    {
#pragma unroll
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir <= i)
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                            else
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - halo_wave.right_sau[i - ir - 1 + geo_const::phy_fdorder_half + (j + phy_n_cols - sub_pad.pad_back) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                            }
                        }
                    }
                    else if (!is_pos.is_right)
                    {
                        if (i <= start_fdorder - 1)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                        else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                        {
                            for (ir = 0; ir <= i; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                }
                else if (i > n_rows - sub_pad.pad_left - 1 && i < n_rows)
                {
                    if (is_pos.is_left)
                    {
#pragma unroll
                        for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                        {
                            if (ir < n_rows - 1 - i)
                            {
                                sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                            else
                            {
                                sau_x_tmp += phy_dc[ir] * (halo_wave.left_sau[i + ir + 1 - n_rows + (j + phy_n_cols - sub_pad.pad_back) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                    else if (!is_pos.is_left)
                    {
                        int _iw = i - n_rows + geo_const::phy_fdorder_half;
                        //
                        if (_iw >= geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < start_fdorder; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                        else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                        {
                            for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                            {
                                sau_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                            }
                        }
                    }
                }
                //********************
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    sau_z_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                }
                rtm_wave.vxp[ext_idx] += a * sau_x_tmp;
                rtm_wave.vyp[ext_idx] += b * sau_y_tmp;
                rtm_wave.vzp[ext_idx] += c * sau_z_tmp;
            }
            if (simulate_type == SimulateType::fwi)
            {
                fwi_wave.ux[ext_idx] += dt * org_wave.vx[ext_idx];
                fwi_wave.uy[ext_idx] += dt * org_wave.vy[ext_idx];
                fwi_wave.uz[ext_idx] += dt * org_wave.vz[ext_idx];
                fwi_wave.vxrt[ext_idx] = dt * (a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp);
                fwi_wave.vyrt[ext_idx] = dt * (a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp);
                fwi_wave.vzrt[ext_idx] = dt * (a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp);
            }
        }
    }
    //
    inline __global__ void
    sub_update_vel_phy_domain_edge_right(SimulateType simulate_type,
                                         Frame sub_gridext, Frame sub_gridphy, Frame grid_right,
                                         Padding sub_pad,
                                         isPosition is_pos,
                                         int start_fdorder,
                                         float dt, float *phy_dc, float *phy_dc_list, float *rho,
                                         //
                                         elasticWave::orgWave org_wave,
                                         elasticWave::rtmWave rtm_wave,
                                         elasticWave::fwiWave fwi_wave,
                                         mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_right);
        int phy_n_cols = sub_gridphy.n_cols;
        float sxx_x_tmp = 0.f, sxy_y_tmp = 0.f, sxz_z_tmp = 0.f,
              syx_x_tmp = 0.f, syy_y_tmp = 0.f, syz_z_tmp = 0.f,
              szx_x_tmp = 0.f, szy_y_tmp = 0.f, szz_z_tmp = 0.f;
        if (idx < grid_right.n_elem)
        {
            float dt_rho = dt / rho[ext_idx];
            float a = dt_rho / sub_gridext.d_rows;
            float b = dt_rho / sub_gridext.d_cols;
            float c = dt_rho / sub_gridext.d_slices;
            int ir;
            if (is_pos.is_right)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir <= i)
                    {
                        sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                    }
                    else
                    {
                        sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - halo_wave.right_sxx[i - ir - 1 + geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                    }
                    if (ir < i)
                    {
                        syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                        szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                    }
                    else
                    {
                        syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - halo_wave.right_sxy[i - ir - 1 + geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols]);
                        szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - halo_wave.right_sxz[i - ir - 1 + geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols]);
                    }
                }
            }
            else if (!is_pos.is_right)
            {
                if (i <= start_fdorder - 1)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        sxx_x_tmp += phy_dc_list[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        syx_x_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                        szx_x_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                    }
                }
                else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                {
                    for (ir = 0; ir <= i; ++ir)
                    {
                        sxx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        syx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                        szx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                    }
                }
            }
            //**************Y
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
            }
            //**************Z
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                sxz_z_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                szz_z_tmp += phy_dc[ir] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
                syz_z_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
            }
            //********************
            org_wave.vx[ext_idx] += a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp;
            org_wave.vy[ext_idx] += a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp;
            org_wave.vz[ext_idx] += a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp;
            //
            if (simulate_type == SimulateType::rtm_backward || simulate_type == SimulateType::rtm_reverse)
            {
                float sau_x_tmp = 0.f, sau_y_tmp = 0.f, sau_z_tmp = 0.f;
                if (is_pos.is_right)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir <= i)
                        {
                            sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                        }
                        else
                        {
                            sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - halo_wave.right_sau[i - ir - 1 + geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                        }
                    }
                }
                else if (!is_pos.is_right)
                {
                    if (i <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sau_x_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                        }
                    }
                    else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= i; ++ir)
                        {
                            sau_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                        }
                    }
                }
                //**************Y
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                }
                //**************Z
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    sau_z_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                }
                rtm_wave.vxp[ext_idx] += a * sau_x_tmp;
                rtm_wave.vyp[ext_idx] += b * sau_y_tmp;
                rtm_wave.vzp[ext_idx] += c * sau_z_tmp;
            }
            if (simulate_type == SimulateType::fwi)
            {
                fwi_wave.ux[ext_idx] += dt * org_wave.vx[ext_idx];
                fwi_wave.uy[ext_idx] += dt * org_wave.vy[ext_idx];
                fwi_wave.uz[ext_idx] += dt * org_wave.vz[ext_idx];
                fwi_wave.vxrt[ext_idx] = dt * (a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp);
                fwi_wave.vyrt[ext_idx] = dt * (a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp);
                fwi_wave.vzrt[ext_idx] = dt * (a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp);
            }
        }
    }
    //
    inline __global__ void
    sub_update_vel_phy_domain_edge_left(SimulateType simulate_type,
                                        Frame sub_gridext, Frame sub_gridphy, Frame grid_left,
                                        Padding sub_pad,
                                        isPosition is_pos,
                                        int start_fdorder,
                                        float dt, float *phy_dc, float *phy_dc_list, float *rho,
                                        //
                                        elasticWave::orgWave org_wave,
                                        elasticWave::rtmWave rtm_wave,
                                        elasticWave::fwiWave fwi_wave,
                                        mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_left);
        int phy_n_cols = sub_gridphy.n_cols;
        float sxx_x_tmp = 0.f, sxy_y_tmp = 0.f, sxz_z_tmp = 0.f,
              syx_x_tmp = 0.f, syy_y_tmp = 0.f, syz_z_tmp = 0.f,
              szx_x_tmp = 0.f, szy_y_tmp = 0.f, szz_z_tmp = 0.f;
        if (idx < grid_left.n_elem)
        {
            float dt_rho = dt / rho[ext_idx];
            float a = dt_rho / sub_gridext.d_rows;
            float b = dt_rho / sub_gridext.d_cols;
            float c = dt_rho / sub_gridext.d_slices;
            int ir;
            if (is_pos.is_left)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir < geo_const::phy_fdorder_half - 1 - i)
                    {
                        sxx_x_tmp += phy_dc[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                    }
                    else
                    {
                        sxx_x_tmp += phy_dc[ir] * (halo_wave.left_sxx[i + ir + 1 - geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.sxx[ext_idx - ir]);
                    }
                    if (ir <= geo_const::phy_fdorder_half - 1 - i)
                    {
                        syx_x_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                        szx_x_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                    }
                    else
                    {
                        syx_x_tmp += phy_dc[ir] * (halo_wave.left_sxy[i + ir - geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.sxy[ext_idx - ir - 1]);
                        szx_x_tmp += phy_dc[ir] * (halo_wave.left_sxz[i + ir - geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.sxz[ext_idx - ir - 1]);
                    }
                }
            }
            else
            {
                if (i >= geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        sxx_x_tmp += phy_dc_list[ir] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        syx_x_tmp += phy_dc_list[ir] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                        szx_x_tmp += phy_dc_list[ir] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                    }
                }
                else if (i < geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < geo_const::phy_fdorder_half - i; ++ir)
                    {
                        sxx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - i - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxx[ext_idx + ir + 1] - org_wave.sxx[ext_idx - ir]);
                        syx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - i - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxy[ext_idx + ir] - org_wave.sxy[ext_idx - ir - 1]);
                        szx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - i - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.sxz[ext_idx + ir] - org_wave.sxz[ext_idx - ir - 1]);
                    }
                }
            }
            //********************Y
            //*
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                sxy_y_tmp += phy_dc[ir] * (org_wave.sxy[ext_idx + ir * ext_n_rows] - org_wave.sxy[ext_idx - (ir + 1) * ext_n_rows]);
                syy_y_tmp += phy_dc[ir] * (org_wave.syy[ext_idx + (ir + 1) * ext_n_rows] - org_wave.syy[ext_idx - ir * ext_n_rows]);
                szy_y_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_rows] - org_wave.syz[ext_idx - (ir + 1) * ext_n_rows]);
            }
            //*
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                sxz_z_tmp += phy_dc[ir] * (org_wave.sxz[ext_idx + ir * ext_n_elem_slice] - org_wave.sxz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                syz_z_tmp += phy_dc[ir] * (org_wave.syz[ext_idx + ir * ext_n_elem_slice] - org_wave.syz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                szz_z_tmp += phy_dc[ir] * (org_wave.szz[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.szz[ext_idx - ir * ext_n_elem_slice]);
            }
            //********************
            org_wave.vx[ext_idx] += a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp;
            org_wave.vy[ext_idx] += a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp;
            org_wave.vz[ext_idx] += a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp;
            //
            if (simulate_type == SimulateType::rtm_backward || simulate_type == SimulateType::rtm_reverse)
            {
                float sau_x_tmp = 0.f, sau_y_tmp = 0.f, sau_z_tmp = 0.f;
                if (is_pos.is_left)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < geo_const::phy_fdorder_half - 1 - i)
                        {
                            sau_x_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                        }
                        else
                        {
                            sau_x_tmp += phy_dc[ir] * (halo_wave.left_sau[i + ir + 1 - geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - rtm_wave.sau[ext_idx - ir]);
                        }
                    }
                }
                else if (!is_pos.is_left)
                {
                    if (i >= geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            sau_x_tmp += phy_dc_list[ir] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                        }
                    }
                    else if (i < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - i; ++ir)
                        {
                            sau_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - i - start_fdorder) * geo_const::phy_fdorder_half] * (rtm_wave.sau[ext_idx + ir + 1] - rtm_wave.sau[ext_idx - ir]);
                        }
                    }
                }
                //********************Y
                //*
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    sau_y_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_rows] - rtm_wave.sau[ext_idx - ir * ext_n_rows]);
                }
                //*
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    sau_z_tmp += phy_dc[ir] * (rtm_wave.sau[ext_idx + (ir + 1) * ext_n_elem_slice] - rtm_wave.sau[ext_idx - ir * ext_n_elem_slice]);
                }
                rtm_wave.vxp[ext_idx] += a * sau_x_tmp;
                rtm_wave.vyp[ext_idx] += b * sau_y_tmp;
                rtm_wave.vzp[ext_idx] += c * sau_z_tmp;
            }
            if (simulate_type == SimulateType::fwi)
            {
                fwi_wave.ux[ext_idx] += dt * org_wave.vx[ext_idx];
                fwi_wave.uy[ext_idx] += dt * org_wave.vy[ext_idx];
                fwi_wave.uz[ext_idx] += dt * org_wave.vz[ext_idx];
                fwi_wave.vxrt[ext_idx] = dt * (a * sxx_x_tmp + b * sxy_y_tmp + c * sxz_z_tmp);
                fwi_wave.vyrt[ext_idx] = dt * (a * syx_x_tmp + b * syy_y_tmp + c * syz_z_tmp);
                fwi_wave.vzrt[ext_idx] = dt * (a * szx_x_tmp + b * szy_y_tmp + c * szz_z_tmp);
            }
        }
    }
}
#endif