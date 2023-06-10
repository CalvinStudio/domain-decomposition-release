#pragma once
#ifndef _FDTD3D_PHY_MPI_STRESS_KERNEL_HPP
#define _FDTD3D_PHY_MPI_STRESS_KERNEL_HPP
#include "fdtd3d_phy_mpi_2_vel_kernel.hpp"
namespace jarvis
{
    inline __global__ void
    sub_update_stress_phy_domain_inside(SimulateType simulate_type, Frame sub_gridext, Frame grid_inside,
                                        float dt, float *phy_dc, float *lambda, float *mu,
                                        elasticWave::orgWave org_wave, elasticWave::rtmWave rtm_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_inside);
        float a = dt / sub_gridext.d_rows;
        float b = dt / sub_gridext.d_cols;
        float c = dt / sub_gridext.d_slices;
        float vx_x_tmp = 0.f, vy_y_tmp = 0.f, vz_z_tmp = 0.f,
              vx_y_tmp = 0.f, vx_z_tmp = 0.f,
              vy_x_tmp = 0.f, vy_z_tmp = 0.f,
              vz_x_tmp = 0.f, vz_y_tmp = 0.f;
        if (idx < grid_inside.n_elem)
        {
            float lambda_tmp = lambda[ext_idx];
            float mu_tmp = mu[ext_idx];
            float lambda_2mu_tmp = lambda_tmp + 2.0f * mu_tmp;
            int ir;
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                vx_z_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
            }
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                vy_z_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
            }
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                vz_z_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
            }
            org_wave.sxx[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.syy[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.szz[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            org_wave.sxy[ext_idx] += a * mu_tmp * vy_x_tmp + b * mu_tmp * vx_y_tmp;
            org_wave.sxz[ext_idx] += a * mu_tmp * vz_x_tmp + c * mu_tmp * vx_z_tmp;
            org_wave.syz[ext_idx] += b * mu_tmp * vz_y_tmp + c * mu_tmp * vy_z_tmp;
            if (simulate_type == SimulateType::rtm_forward ||
                simulate_type == SimulateType::rtm_backward ||
                simulate_type == SimulateType::rtm_reverse)
            {
                rtm_wave.sau[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            }
        }
    }
    //
    inline __global__ void
    sub_update_stress_phy_domain_edge_top(SimulateType simulate_type,
                                          Frame sub_gridext, Frame sub_gridphy, Frame grid_top,
                                          Padding sub_pad,
                                          isPosition is_pos,
                                          int start_fdorder,
                                          float dt, float *phy_dc, float *phy_dc_list,
                                          float *lambda, float *mu,
                                          elasticWave::orgWave org_wave,
                                          elasticWave::rtmWave rtm_wave,
                                          mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_top);
        int phy_n_cols = sub_gridphy.n_cols;
        float a = dt / sub_gridext.d_rows;
        float b = dt / sub_gridext.d_cols;
        float c = dt / sub_gridext.d_slices;
        float lambda_tmp, mu_tmp, lambda_2mu_tmp;
        float vx_x_tmp = 0, vy_y_tmp = 0, vz_z_tmp = 0,
              vy_x_tmp = 0, vx_y_tmp = 0,
              vz_x_tmp = 0, vx_z_tmp = 0,
              vz_y_tmp = 0, vy_z_tmp = 0;
        if (idx < grid_top.n_elem)
        {
            lambda_tmp = lambda[ext_idx];
            mu_tmp = mu[ext_idx];
            lambda_2mu_tmp = lambda_tmp + 2.0 * mu_tmp;
            int ir;
            if (is_pos.is_top)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir <= k)
                    {
                        vx_z_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
                        vy_z_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
                    }
                    else
                    {
                        vx_z_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - halo_wave.top_vx[i + j * n_rows + (k - ir - 1 + geo_const::phy_fdorder_half) * n_elem_slice]);
                        vy_z_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - halo_wave.top_vy[i + j * n_rows + (k - ir - 1 + geo_const::phy_fdorder_half) * n_elem_slice]);
                    }
                    if (ir < k)
                    {
                        vz_z_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                    }
                    else
                    {
                        vz_z_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - halo_wave.top_vz[i + j * n_rows + (k - ir - 1 + geo_const::phy_fdorder_half) * n_elem_slice]);
                    }
                }
            }
            else if (!is_pos.is_top)
            {
                if (k <= start_fdorder - 1)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        vx_z_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
                        vy_z_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
                        vz_z_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                    }
                }
                else if (k >= start_fdorder && k < geo_const::phy_fdorder_half)
                {
                    for (ir = 0; ir <= k; ++ir)
                    {
                        vx_z_tmp += phy_dc_list[ir + (k - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
                        vy_z_tmp += phy_dc_list[ir + (k - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
                        vz_z_tmp += phy_dc_list[ir + (k - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                    }
                }
            }
            //*X
            if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                    vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                    vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                }
            }
            else if (i < sub_pad.pad_right)
            {
                if (is_pos.is_right)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < i)
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        else
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - halo_wave.right_vx[i - ir - 1 + geo_const::phy_fdorder_half + j * geo_const::phy_fdorder_half + k * geo_const::phy_fdorder_half * phy_n_cols]);
                        }
                        if (ir <= i)
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                        else
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - halo_wave.right_vy[i - ir - 1 + geo_const::phy_fdorder_half + j * (geo_const::phy_fdorder_half - 1) + k * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - halo_wave.right_vz[i - ir - 1 + geo_const::phy_fdorder_half + j * (geo_const::phy_fdorder_half - 1) + k * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                        }
                    }
                }
                else if (!is_pos.is_right)
                {
                    if (i <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            vx_x_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                    }
                    else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= i; ++ir)
                        {
                            vx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
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
                        if (ir <= n_rows - 1 - i)
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        else
                        {
                            vx_x_tmp += phy_dc[ir] * (halo_wave.left_vx[i + ir - n_rows + j * (geo_const::phy_fdorder_half - 1) + k * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        if (ir < n_rows - 1 - i)
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                        else
                        {
                            vy_x_tmp += phy_dc[ir] * (halo_wave.left_vy[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + k * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (halo_wave.left_vz[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + k * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.vz[ext_idx - ir]);
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
                            vy_x_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                            vx_x_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        }
                    }
                    else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                        {
                            vy_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                            vx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        }
                    }
                }
            }
            //*Y
            if (j >= sub_pad.pad_front && j <= phy_n_cols - sub_pad.pad_back - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                    vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                    vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                }
            }
            else if (j < sub_pad.pad_front)
            {
                if (is_pos.is_front)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < j)
                        {
                            vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                        else
                        {
                            vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - halo_wave.front_vy[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + k * n_rows * geo_const::phy_fdorder_half]);
                        }
                        if (ir <= j)
                        {
                            vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                        }
                        else
                        {
                            vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_vx[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + k * n_rows * (geo_const::phy_fdorder_half - 1)]);
                            vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_vz[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + k * n_rows * (geo_const::phy_fdorder_half - 1)]);
                        }
                    }
                }
                else if (!is_pos.is_front)
                {
                    if (j <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            vy_y_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                            vx_y_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                        }
                    }
                    else if (j >= start_fdorder && j < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= j; ++ir)
                        {
                            vy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                            vx_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
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
                        if (ir <= phy_n_cols - 1 - j)
                        {
                            vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                        else
                        {
                            vy_y_tmp += phy_dc[ir] * (halo_wave.back_vy[i + (j + ir - phy_n_cols) * n_rows + k * n_rows * (geo_const::phy_fdorder_half - 1)] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                        if (ir < phy_n_cols - 1 - j)
                        {
                            vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                        }
                        else
                        {
                            vx_y_tmp += phy_dc[ir] * (halo_wave.back_vx[i + (j + ir + 1 - phy_n_cols) * n_rows + k * n_rows * geo_const::phy_fdorder_half] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc[ir] * (halo_wave.back_vz[i + (j + ir + 1 - phy_n_cols) * n_rows + k * n_rows * geo_const::phy_fdorder_half] - org_wave.vz[ext_idx - ir * ext_n_rows]);
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
                            vy_y_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                            vx_y_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                        }
                    }
                    else if (_jw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _jw; ++ir)
                        {
                            vy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                            vx_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                        }
                    }
                }
            }
            //*********************
            org_wave.sxx[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.syy[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.szz[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            org_wave.sxy[ext_idx] += a * mu_tmp * vy_x_tmp + b * mu_tmp * vx_y_tmp;
            org_wave.sxz[ext_idx] += a * mu_tmp * vz_x_tmp + c * mu_tmp * vx_z_tmp;
            org_wave.syz[ext_idx] += b * mu_tmp * vz_y_tmp + c * mu_tmp * vy_z_tmp;
            if (simulate_type == SimulateType::rtm_forward ||
                simulate_type == SimulateType::rtm_backward ||
                simulate_type == SimulateType::rtm_reverse)
            {
                rtm_wave.sau[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            }
        }
    }
    //
    inline __global__ void
    sub_update_stress_phy_domain_edge_bottom(SimulateType simulate_type,
                                             Frame sub_gridext, Frame sub_gridphy, Frame grid_bottom,
                                             Padding sub_pad,
                                             isPosition is_pos,
                                             int start_fdorder,
                                             float dt, float *phy_dc, float *phy_dc_list,
                                             float *lambda, float *mu,
                                             elasticWave::orgWave org_wave,
                                             elasticWave::rtmWave rtm_wave,
                                             mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_bottom);
        int phy_n_cols = sub_gridphy.n_cols;
        int phy_n_slices = sub_gridphy.n_slices;
        float a = dt / sub_gridext.d_rows;
        float b = dt / sub_gridext.d_cols;
        float c = dt / sub_gridext.d_slices;
        float lambda_tmp, mu_tmp, lambda_2mu_tmp;
        float vx_x_tmp = 0.f, vy_y_tmp = 0.f, vz_z_tmp = 0.f,
              vy_x_tmp = 0.f, vx_y_tmp = 0.f,
              vz_x_tmp = 0.f, vx_z_tmp = 0.f,
              vz_y_tmp = 0.f, vy_z_tmp = 0.f;
        if (idx < grid_bottom.n_elem)
        {
            lambda_tmp = lambda[ext_idx];
            mu_tmp = mu[ext_idx];
            lambda_2mu_tmp = lambda_tmp + 2.0 * mu_tmp;
            int ir;
            if (is_pos.is_bottom)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir <= geo_const::phy_fdorder_half - 1 - k)
                    {
                        vz_z_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                    }
                    else
                    {
                        vz_z_tmp += phy_dc[ir] * (halo_wave.bottom_vz[i + j * n_rows + (k + ir - geo_const::phy_fdorder_half) * n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                    }
                    if (ir < geo_const::phy_fdorder_half - 1 - k)
                    {
                        vx_z_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
                        vy_z_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
                    }
                    else
                    {
                        vx_z_tmp += phy_dc[ir] * (halo_wave.bottom_vx[i + j * n_rows + (k + ir + 1 - geo_const::phy_fdorder_half) * n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
                        vy_z_tmp += phy_dc[ir] * (halo_wave.bottom_vy[i + j * n_rows + (k + ir + 1 - geo_const::phy_fdorder_half) * n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
                    }
                }
            }
            else if (!is_pos.is_bottom)
            {
                if (k >= geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        vx_z_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
                        vy_z_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
                        vz_z_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                    }
                }
                else if (k < geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < geo_const::phy_fdorder_half - k; ++ir)
                    {
                        vx_z_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - k - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
                        vy_z_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - k - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
                        vz_z_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - k - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                    }
                }
            }
            //*X
            if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                    vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                    vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                }
            }
            else if (i < sub_pad.pad_right)
            {
                if (is_pos.is_right)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < i)
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        else
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - halo_wave.right_vx[(i - ir - 1 + geo_const::phy_fdorder_half) + j * geo_const::phy_fdorder_half + (k + phy_n_slices - sub_pad.pad_bottom) * geo_const::phy_fdorder_half * phy_n_cols]);
                        }
                        if (ir <= i)
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                        else
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - halo_wave.right_vy[(i - ir - 1 + geo_const::phy_fdorder_half) + j * (geo_const::phy_fdorder_half - 1) + (k + phy_n_slices - sub_pad.pad_bottom) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - halo_wave.right_vz[(i - ir - 1 + geo_const::phy_fdorder_half) + j * (geo_const::phy_fdorder_half - 1) + (k + phy_n_slices - sub_pad.pad_bottom) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                        }
                    }
                }
                else if (!is_pos.is_right)
                {
                    if (i <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            vx_x_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                    }
                    else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= i; ++ir)
                        {
                            vx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
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
                        if (ir <= n_rows - 1 - i)
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        else
                        {
                            vx_x_tmp += phy_dc[ir] * (halo_wave.left_vx[i + ir - n_rows + j * (geo_const::phy_fdorder_half - 1) + (k + phy_n_slices - sub_pad.pad_bottom) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        if (ir < n_rows - 1 - i)
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                        else
                        {
                            vy_x_tmp += phy_dc[ir] * (halo_wave.left_vy[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + (k + phy_n_slices - sub_pad.pad_bottom) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (halo_wave.left_vz[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + (k + phy_n_slices - sub_pad.pad_bottom) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.vz[ext_idx - ir]);
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
                            vx_x_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                    }
                    else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                        {
                            vx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                    }
                }
            }
            //*Y
            if (j >= sub_pad.pad_front && j <= phy_n_cols - sub_pad.pad_back - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                    vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                    vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                }
            }
            else if (j < sub_pad.pad_front)
            {
                if (is_pos.is_front)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < j)
                        {
                            vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                        else
                        {
                            vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - halo_wave.front_vy[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * geo_const::phy_fdorder_half]);
                        }
                        if (ir <= j)
                        {
                            vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                        }
                        else
                        {
                            vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_vx[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * (geo_const::phy_fdorder_half - 1)]);
                            vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_vz[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * (geo_const::phy_fdorder_half - 1)]);
                        }
                    }
                }
                else if (!is_pos.is_front)
                {
                    if (j <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            vy_y_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                            vx_y_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                        }
                    }
                    else if (j >= start_fdorder && j < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= j; ++ir)
                        {
                            vy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                            vx_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
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
                        if (ir <= phy_n_cols - 1 - j)
                        {
                            vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                        else
                        {
                            vy_y_tmp += phy_dc[ir] * (halo_wave.back_vy[i + (j + ir - phy_n_cols) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * (geo_const::phy_fdorder_half - 1)] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                        }
                        if (ir < phy_n_cols - 1 - j)
                        {
                            vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                        }
                        else
                        {
                            vx_y_tmp += phy_dc[ir] * (halo_wave.back_vx[i + (j + ir + 1 - phy_n_cols) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * geo_const::phy_fdorder_half] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc[ir] * (halo_wave.back_vz[i + (j + ir + 1 - phy_n_cols) * n_rows + (k + phy_n_slices - sub_pad.pad_bottom) * n_rows * geo_const::phy_fdorder_half] - org_wave.vz[ext_idx - ir * ext_n_rows]);
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
                            vy_y_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                            vx_y_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                        }
                    }
                    else if (_jw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _jw; ++ir)
                        {
                            vy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                            vx_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                            vz_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _jw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                        }
                    }
                }
            }
            org_wave.sxx[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.syy[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.szz[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            org_wave.sxy[ext_idx] += a * mu_tmp * vy_x_tmp + b * mu_tmp * vx_y_tmp;
            org_wave.sxz[ext_idx] += a * mu_tmp * vz_x_tmp + c * mu_tmp * vx_z_tmp;
            org_wave.syz[ext_idx] += b * mu_tmp * vz_y_tmp + c * mu_tmp * vy_z_tmp;
            if (simulate_type == SimulateType::rtm_forward ||
                simulate_type == SimulateType::rtm_backward ||
                simulate_type == SimulateType::rtm_reverse)
            {
                rtm_wave.sau[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            }
        }
    }
    //
    //
    inline __global__ void
    sub_update_stress_phy_domain_edge_front(SimulateType simulate_type,
                                            Frame sub_gridext, Frame sub_gridphy, Frame grid_front,
                                            Padding sub_pad,
                                            isPosition is_pos,
                                            int start_fdorder,
                                            float dt, float *phy_dc, float *phy_dc_list,
                                            float *lambda, float *mu,
                                            elasticWave::orgWave org_wave,
                                            elasticWave::rtmWave rtm_wave,
                                            mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_front);
        int phy_n_cols = sub_gridphy.n_cols;
        float a = dt / sub_gridext.d_rows;
        float b = dt / sub_gridext.d_cols;
        float c = dt / sub_gridext.d_slices;
        float lambda_tmp, mu_tmp, lambda_2mu_tmp;
        float vx_x_tmp = 0.f, vy_y_tmp = 0.f, vz_z_tmp = 0.f,
              vy_x_tmp = 0.f, vx_y_tmp = 0.f,
              vz_x_tmp = 0.f, vx_z_tmp = 0.f,
              vz_y_tmp = 0.f, vy_z_tmp = 0.f;
        if (idx < grid_front.n_elem)
        {
            lambda_tmp = lambda[ext_idx];
            mu_tmp = mu[ext_idx];
            lambda_2mu_tmp = lambda_tmp + 2.0 * mu_tmp;
            int ir;
            if (is_pos.is_front)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir < j)
                    {
                        vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                    }
                    else
                    {
                        vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - halo_wave.front_vy[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * geo_const::phy_fdorder_half]);
                    }
                    if (ir <= j)
                    {
                        vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                        vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                    }
                    else
                    {
                        vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_vx[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * (geo_const::phy_fdorder_half - 1)]);
                        vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.front_vz[i + (j - ir - 1 + geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * (geo_const::phy_fdorder_half - 1)]);
                    }
                }
            }
            else if (!is_pos.is_front)
            {
                if (j <= start_fdorder - 1)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        vy_y_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                        vx_y_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                        vz_y_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                    }
                }
                else if (j >= start_fdorder && j < geo_const::phy_fdorder_half)
                {
                    for (ir = 0; ir <= j; ++ir)
                    {
                        vy_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                        vx_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                        vz_y_tmp += phy_dc_list[ir + (j - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                    }
                }
            }
            //*X
            if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                    vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                    vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                }
            }
            else if (i < sub_pad.pad_right)
            {
                if (is_pos.is_right)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < i)
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        else
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - halo_wave.right_vx[i - ir - 1 + geo_const::phy_fdorder_half + j * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols]);
                        }
                        if (ir <= i)
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                        else
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.right_vy[i - ir - 1 + geo_const::phy_fdorder_half + j * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - halo_wave.right_vz[i - ir - 1 + geo_const::phy_fdorder_half + j * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                        }
                    }
                }
                else if (!is_pos.is_right)
                {
                    if (i <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            vx_x_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                    }
                    else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= i; ++ir)
                        {
                            vx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
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
                        if (ir <= n_rows - 1 - i)
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        else
                        {
                            vx_x_tmp += phy_dc[ir] * (halo_wave.left_vx[i + ir - n_rows + j * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        if (ir < n_rows - 1 - i)
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                        else
                        {
                            vy_x_tmp += phy_dc[ir] * (halo_wave.left_vy[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (halo_wave.left_vz[i + ir + 1 - n_rows + j * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.vz[ext_idx - ir]);
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
                            vx_x_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                    }
                    else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                        {
                            vx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                    }
                }
            }
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                vz_z_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                vx_z_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
                vy_z_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
            }
            org_wave.sxx[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.syy[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.szz[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            org_wave.sxy[ext_idx] += a * mu_tmp * vy_x_tmp + b * mu_tmp * vx_y_tmp;
            org_wave.sxz[ext_idx] += a * mu_tmp * vz_x_tmp + c * mu_tmp * vx_z_tmp;
            org_wave.syz[ext_idx] += b * mu_tmp * vz_y_tmp + c * mu_tmp * vy_z_tmp;
            if (simulate_type == SimulateType::rtm_forward ||
                simulate_type == SimulateType::rtm_backward ||
                simulate_type == SimulateType::rtm_reverse)
            {
                rtm_wave.sau[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            }
        }
    }
    //
    inline __global__ void
    sub_update_stress_phy_domain_edge_back(SimulateType simulate_type,
                                           Frame sub_gridext, Frame sub_gridphy, Frame grid_back,
                                           Padding sub_pad,
                                           isPosition is_pos,
                                           int start_fdorder,
                                           float dt, float *phy_dc, float *phy_dc_list, float *lambda, float *mu,
                                           elasticWave::orgWave org_wave,
                                           elasticWave::rtmWave rtm_wave,
                                           mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_back);
        int phy_n_cols = sub_gridphy.n_cols;
        float a = dt / sub_gridext.d_rows;
        float b = dt / sub_gridext.d_cols;
        float c = dt / sub_gridext.d_slices;
        float lambda_tmp, mu_tmp, lambda_2mu_tmp;
        float vx_x_tmp = 0.f, vy_y_tmp = 0.f, vz_z_tmp = 0.f,
              vy_x_tmp = 0.f, vx_y_tmp = 0.f,
              vz_x_tmp = 0.f, vx_z_tmp = 0.f,
              vz_y_tmp = 0.f, vy_z_tmp = 0.f;
        if (idx < grid_back.n_elem)
        {
            lambda_tmp = lambda[ext_idx];
            mu_tmp = mu[ext_idx];
            lambda_2mu_tmp = lambda_tmp + 2.0 * mu_tmp;
            int ir;
            if (is_pos.is_back)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir < geo_const::phy_fdorder_half - 1 - j)
                    {
                        vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                        vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                    }
                    else
                    {
                        vx_y_tmp += phy_dc[ir] * (halo_wave.back_vx[i + (j + ir + 1 - geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * geo_const::phy_fdorder_half] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                        vz_y_tmp += phy_dc[ir] * (halo_wave.back_vz[i + (j + ir + 1 - geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * geo_const::phy_fdorder_half] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                    }
                    if (ir <= geo_const::phy_fdorder_half - 1 - j)
                    {
                        vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                    }
                    else
                    {
                        vy_y_tmp += phy_dc[ir] * (halo_wave.back_vy[i + (j + ir - geo_const::phy_fdorder_half) * n_rows + (k + sub_pad.pad_top) * n_rows * (geo_const::phy_fdorder_half - 1)] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                    }
                }
            }
            else if (!is_pos.is_back)
            {
                if (j >= geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        vy_y_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                        vx_y_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                        vz_y_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                    }
                }
                else if (j < geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < geo_const::phy_fdorder_half - j; ++ir)
                    {
                        vy_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - j - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                        vx_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - j - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                        vz_y_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - j - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
                    }
                }
            }
            //*X
            if (i >= sub_pad.pad_right && i <= n_rows - sub_pad.pad_left - 1)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                    vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                    vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                }
            }
            else if (i < sub_pad.pad_right)
            {
                if (is_pos.is_right)
                {
#pragma unroll
                    for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                    {
                        if (ir < i)
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        else
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - halo_wave.right_vx[i - ir - 1 + geo_const::phy_fdorder_half + (j + phy_n_cols - sub_pad.pad_back) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols]);
                        }
                        if (ir <= i)
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                        else
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - halo_wave.right_vy[i - ir - 1 + geo_const::phy_fdorder_half + (j + phy_n_cols - sub_pad.pad_back) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - halo_wave.right_vz[i - ir - 1 + geo_const::phy_fdorder_half + (j + phy_n_cols - sub_pad.pad_back) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                        }
                    }
                }
                else if (!is_pos.is_right)
                {
                    if (i <= start_fdorder - 1)
                    {
                        for (ir = 0; ir < start_fdorder; ++ir)
                        {
                            vx_x_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                    }
                    else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                    {
                        for (ir = 0; ir <= i; ++ir)
                        {
                            vx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
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
                        if (ir <= n_rows - 1 - i)
                        {
                            vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        else
                        {
                            vx_x_tmp += phy_dc[ir] * (halo_wave.left_vx[i + ir - n_rows + (j + phy_n_cols - sub_pad.pad_back) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.vx[ext_idx - ir - 1]);
                        }
                        if (ir < n_rows - 1 - i)
                        {
                            vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                        else
                        {
                            vy_x_tmp += phy_dc[ir] * (halo_wave.left_vy[i + ir + 1 - n_rows + (j + phy_n_cols - sub_pad.pad_back) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc[ir] * (halo_wave.left_vz[i + ir + 1 - n_rows + (j + phy_n_cols - sub_pad.pad_back) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.vz[ext_idx - ir]);
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
                            vx_x_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                    }
                    else if (_iw < geo_const::phy_fdorder_half - start_fdorder)
                    {
                        for (ir = 0; ir < geo_const::phy_fdorder_half - _iw; ++ir)
                        {
                            vx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                            vy_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                            vz_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - _iw - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                        }
                    }
                }
            }
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                vz_z_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
                vx_z_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
                vy_z_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
            }
            org_wave.sxx[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.syy[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.szz[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            org_wave.sxy[ext_idx] += a * mu_tmp * vy_x_tmp + b * mu_tmp * vx_y_tmp;
            org_wave.sxz[ext_idx] += a * mu_tmp * vz_x_tmp + c * mu_tmp * vx_z_tmp;
            org_wave.syz[ext_idx] += b * mu_tmp * vz_y_tmp + c * mu_tmp * vy_z_tmp;
            if (simulate_type == SimulateType::rtm_forward ||
                simulate_type == SimulateType::rtm_backward ||
                simulate_type == SimulateType::rtm_reverse)
            {
                rtm_wave.sau[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            }
        }
    }
    //
    inline __global__ void
    sub_update_stress_phy_domain_edge_right(SimulateType simulate_type,
                                            Frame sub_gridext, Frame sub_gridphy, Frame grid_right,
                                            Padding sub_pad,
                                            isPosition is_pos,
                                            int start_fdorder,
                                            float dt, float *phy_dc, float *phy_dc_list,
                                            float *lambda, float *mu,
                                            elasticWave::orgWave org_wave,
                                            elasticWave::rtmWave rtm_wave,
                                            mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_right);
        int phy_n_cols = sub_gridphy.n_cols;
        float a = dt / sub_gridext.d_rows;
        float b = dt / sub_gridext.d_cols;
        float c = dt / sub_gridext.d_slices;
        float lambda_tmp, mu_tmp, lambda_2mu_tmp;
        float vx_x_tmp = 0.f, vy_y_tmp = 0.f, vz_z_tmp = 0.f,
              vy_x_tmp = 0.f, vx_y_tmp = 0.f,
              vz_x_tmp = 0.f, vx_z_tmp = 0.f,
              vz_y_tmp = 0.f, vy_z_tmp = 0.f;
        if (idx < grid_right.n_elem)
        {
            lambda_tmp = lambda[ext_idx];
            mu_tmp = mu[ext_idx];
            lambda_2mu_tmp = lambda_tmp + 2.0 * mu_tmp;
            int ir;
            if (is_pos.is_right)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir < i)
                    {
                        vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                    }
                    else
                    {
                        vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - halo_wave.right_vx[i - ir - 1 + geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols]);
                    }
                    if (ir <= i)
                    {
                        vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                        vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                    }
                    else
                    {
                        vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - halo_wave.right_vy[i - ir - 1 + geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                        vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - halo_wave.right_vz[i - ir - 1 + geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols]);
                    }
                }
            }
            else if (!is_pos.is_right)
            {
                if (i <= start_fdorder - 1)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        vx_x_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        vy_x_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                        vz_x_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                    }
                }
                else if (i >= start_fdorder && i < geo_const::phy_fdorder_half)
                {
                    for (ir = 0; ir <= i; ++ir)
                    {
                        vx_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        vy_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                        vz_x_tmp += phy_dc_list[ir + (i - start_fdorder + 1) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                    }
                }
            }
            //*Y
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
            }
            //*Z
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                vx_z_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
                vy_z_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
                vz_z_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
            }
            org_wave.sxx[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.syy[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.szz[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            org_wave.sxy[ext_idx] += a * mu_tmp * vy_x_tmp + b * mu_tmp * vx_y_tmp;
            org_wave.sxz[ext_idx] += a * mu_tmp * vz_x_tmp + c * mu_tmp * vx_z_tmp;
            org_wave.syz[ext_idx] += b * mu_tmp * vz_y_tmp + c * mu_tmp * vy_z_tmp;
            if (simulate_type == SimulateType::rtm_forward ||
                simulate_type == SimulateType::rtm_backward ||
                simulate_type == SimulateType::rtm_reverse)
            {
                rtm_wave.sau[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            }
        }
    }
    //
    inline __global__ void
    sub_update_stress_phy_domain_edge_left(SimulateType simulate_type,
                                           Frame sub_gridext, Frame sub_gridphy, Frame grid_left,
                                           Padding sub_pad,
                                           isPosition is_pos,
                                           int start_fdorder,
                                           float dt, float *phy_dc, float *phy_dc_list,
                                           float *lambda, float *mu,
                                           elasticWave::orgWave org_wave,
                                           elasticWave::rtmWave rtm_wave,
                                           mpiHaloWave_p halo_wave)
    {
        set_cufield_grid_3d_idx(sub_gridext, grid_left);
        int phy_n_cols = sub_gridphy.n_cols;
        float a = dt / sub_gridext.d_rows;
        float b = dt / sub_gridext.d_cols;
        float c = dt / sub_gridext.d_slices;
        float lambda_tmp, mu_tmp, lambda_2mu_tmp;
        float vx_x_tmp = 0.f, vy_y_tmp = 0.f, vz_z_tmp = 0.f,
              vy_x_tmp = 0.f, vx_y_tmp = 0.f,
              vz_x_tmp = 0.f, vx_z_tmp = 0.f,
              vz_y_tmp = 0.f, vy_z_tmp = 0.f;
        if (idx < grid_left.n_elem)
        {
            lambda_tmp = lambda[ext_idx];
            mu_tmp = mu[ext_idx];
            lambda_2mu_tmp = lambda_tmp + 2.0 * mu_tmp;
            int ir;
            if (is_pos.is_left)
            {
#pragma unroll
                for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
                {
                    if (ir <= geo_const::phy_fdorder_half - 1 - i)
                    {
                        vx_x_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                    }
                    else
                    {
                        vx_x_tmp += phy_dc[ir] * (halo_wave.left_vx[i + ir - geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * (geo_const::phy_fdorder_half - 1) + (k + sub_pad.pad_top) * (geo_const::phy_fdorder_half - 1) * phy_n_cols] - org_wave.vx[ext_idx - ir - 1]);
                    }
                    if (ir < geo_const::phy_fdorder_half - 1 - i)
                    {
                        vy_x_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                        vz_x_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                    }
                    else
                    {
                        vy_x_tmp += phy_dc[ir] * (halo_wave.left_vy[i + ir + 1 - geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.vy[ext_idx - ir]);
                        vz_x_tmp += phy_dc[ir] * (halo_wave.left_vz[i + ir + 1 - geo_const::phy_fdorder_half + (j + sub_pad.pad_front) * geo_const::phy_fdorder_half + (k + sub_pad.pad_top) * geo_const::phy_fdorder_half * phy_n_cols] - org_wave.vz[ext_idx - ir]);
                    }
                }
            }
            else if (!is_pos.is_left)
            {
                if (i >= geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < start_fdorder; ++ir)
                    {
                        vx_x_tmp += phy_dc_list[ir] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        vy_x_tmp += phy_dc_list[ir] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                        vz_x_tmp += phy_dc_list[ir] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                    }
                }
                else if (i < geo_const::phy_fdorder_half - start_fdorder)
                {
                    for (ir = 0; ir < geo_const::phy_fdorder_half - i; ++ir)
                    {
                        vx_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - i - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vx[ext_idx + ir] - org_wave.vx[ext_idx - ir - 1]);
                        vy_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - i - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vy[ext_idx + ir + 1] - org_wave.vy[ext_idx - ir]);
                        vz_x_tmp += phy_dc_list[ir + (geo_const::phy_fdorder_half - i - start_fdorder) * geo_const::phy_fdorder_half] * (org_wave.vz[ext_idx + ir + 1] - org_wave.vz[ext_idx - ir]);
                    }
                }
            }
            //*Y
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                vx_y_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vx[ext_idx - ir * ext_n_rows]);
                vy_y_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + ir * ext_n_rows] - org_wave.vy[ext_idx - (ir + 1) * ext_n_rows]);
                vz_y_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + (ir + 1) * ext_n_rows] - org_wave.vz[ext_idx - ir * ext_n_rows]);
            }
#pragma unroll
            for (ir = 0; ir < geo_const::phy_fdorder_half; ++ir)
            {
                vx_z_tmp += phy_dc[ir] * (org_wave.vx[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vx[ext_idx - ir * ext_n_elem_slice]);
                vy_z_tmp += phy_dc[ir] * (org_wave.vy[ext_idx + (ir + 1) * ext_n_elem_slice] - org_wave.vy[ext_idx - ir * ext_n_elem_slice]);
                vz_z_tmp += phy_dc[ir] * (org_wave.vz[ext_idx + ir * ext_n_elem_slice] - org_wave.vz[ext_idx - (ir + 1) * ext_n_elem_slice]);
            }
            org_wave.sxx[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.syy[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_tmp * vz_z_tmp;
            org_wave.szz[ext_idx] += a * lambda_tmp * vx_x_tmp + b * lambda_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            org_wave.sxy[ext_idx] += a * mu_tmp * vy_x_tmp + b * mu_tmp * vx_y_tmp;
            org_wave.sxz[ext_idx] += a * mu_tmp * vz_x_tmp + c * mu_tmp * vx_z_tmp;
            org_wave.syz[ext_idx] += b * mu_tmp * vz_y_tmp + c * mu_tmp * vy_z_tmp;
            if (simulate_type == SimulateType::rtm_forward ||
                simulate_type == SimulateType::rtm_backward ||
                simulate_type == SimulateType::rtm_reverse)
            {
                rtm_wave.sau[ext_idx] += a * lambda_2mu_tmp * vx_x_tmp + b * lambda_2mu_tmp * vy_y_tmp + c * lambda_2mu_tmp * vz_z_tmp;
            }
        }
    }
    //
    inline __global__ void
    sub_add_stress_source(int it, int iShot, Frame sub_gridext, Point3D *shot, float *wave,
                          float *sxx, float *syy, float *szz)
    {
        set_frame_nd(sub_gridext);
        float l_rows, l_cols, l_slices, r_rows, r_cols, r_slices;
        l_rows = sub_gridext.l_rows;
        l_cols = sub_gridext.l_cols;
        l_slices = sub_gridext.l_slices;
        r_rows = sub_gridext.r_rows;
        r_cols = sub_gridext.r_cols;
        r_slices = sub_gridext.r_slices;
        float shot_x = shot[iShot].x;
        float shot_y = shot[iShot].y;
        float shot_z = shot[iShot].z;

        if (shot_x >= l_rows && shot_y >= l_cols && shot_z >= l_slices &&
            shot_x <= r_rows && shot_y <= r_cols && shot_z <= r_slices)
        {
            int i_loc = get_gap_num(l_rows, shot_x, d_rows);
            int j_loc = get_gap_num(l_cols, shot_y, d_cols);
            int k_loc = get_gap_num(l_slices, shot_z, d_slices);
            int idx = i_loc + j_loc * n_rows + k_loc * n_elem_slice;
            sxx[idx] += wave[it];
            syy[idx] += wave[it];
            szz[idx] += wave[it];
        }
    }
    //
    inline __global__ void
    sub_subtract_stress_source(int it, int iShot, Frame sub_gridext, Point3D *shot, float *wave,
                               float *sxx, float *syy, float *szz)
    {
        set_frame_nd(sub_gridext);
        float l_rows, l_cols, l_slices, r_rows, r_cols, r_slices;
        l_rows = sub_gridext.l_rows;
        l_cols = sub_gridext.l_cols;
        l_slices = sub_gridext.l_slices;
        r_rows = sub_gridext.r_rows;
        r_cols = sub_gridext.r_cols;
        r_slices = sub_gridext.r_slices;
        float shot_x = shot[iShot].x;
        float shot_y = shot[iShot].y;
        float shot_z = shot[iShot].z;

        if (shot_x >= l_rows && shot_y >= l_cols && shot_z >= l_slices &&
            shot_x <= r_rows && shot_y <= r_cols && shot_z <= r_slices)
        {
            int i_loc = get_gap_num(l_rows, shot_x, d_rows);
            int j_loc = get_gap_num(l_cols, shot_y, d_cols);
            int k_loc = get_gap_num(l_slices, shot_z, d_slices);
            int idx = i_loc + j_loc * n_rows + k_loc * n_elem_slice;
            sxx[idx] -= wave[it];
            syy[idx] -= wave[it];
            szz[idx] -= wave[it];
        }
    }
}
#endif