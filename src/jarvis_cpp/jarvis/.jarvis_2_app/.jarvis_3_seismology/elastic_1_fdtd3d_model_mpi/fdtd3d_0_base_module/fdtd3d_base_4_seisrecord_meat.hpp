#pragma once
#ifndef _FDTD3D_SEISRECORD_MEAT_HPP
#define _FDTD3D_SEISRECORD_MEAT_HPP
#include "fdtd3d_base_4_seisrecord_bones.hpp"
namespace jarvis
{
    inline __global__ void
    cal_sub_seis_data_aux(Frame gridpml, Margin margin,
                          int it, int ntime, int ReceNum, Point3D *rece,
                          float *vxp, float *vyp, float *vzp,
                          float *vxp_seis, float *vyp_seis, float *vzp_seis)
    {
        set_frame_ndl(gridpml);
        set_model_padding(margin);
        int iRece = threadIdx.x;
        int idx_seis = it + iRece * ntime;
        int Rx, Ry, Rz, idx_rece;
        if (iRece < ReceNum)
        {
            Rx = get_gap_num(l_rows, rece[iRece].x, d_rows);
            Ry = get_gap_num(l_cols, rece[iRece].y, d_cols);
            Rz = get_gap_num(l_slices, rece[iRece].z, d_slices);
            idx_rece = Rx + Ry * n_rows + Rz * n_elem_slice;
            vxp_seis[idx_seis] = vxp[idx_rece];
            vyp_seis[idx_seis] = vyp[idx_rece];
            vzp_seis[idx_seis] = vzp[idx_rece];
        }
    }

    inline __global__ void
    sub_cal_seis_data_rtm(Frame gridpml, Margin margin,
                          int it, int ntime, int ReceNum, Point3D *rece,
                          float *vx, float *vy, float *vz,
                          float *vx_seis, float *vy_seis, float *vz_seis,
                          float *vxp, float *vyp, float *vzp,
                          float *vxp_seis, float *vyp_seis, float *vzp_seis)
    {
        set_frame_ndl(gridpml);
        int iRece = threadIdx.x;
        int idx_seis = it + iRece * ntime;
        int Rx, Ry, Rz, idx_rece;
        if (iRece < ReceNum)
        {
            Rx = get_gap_num(l_rows, rece[iRece].x, d_rows);
            Ry = get_gap_num(l_cols, rece[iRece].y, d_cols);
            Rz = get_gap_num(l_slices, rece[iRece].z, d_slices);
            idx_rece = Rx + Ry * n_rows + Rz * n_elem_slice;
            vx_seis[idx_seis] = vx[idx_rece];
            vy_seis[idx_seis] = vy[idx_rece];
            vz_seis[idx_seis] = vz[idx_rece];
            vxp_seis[idx_seis] = vxp[idx_rece];
            vyp_seis[idx_seis] = vyp[idx_rece];
            vzp_seis[idx_seis] = vzp[idx_rece];
        }
    }

    inline __global__ void
    sub_cal_seis_data_base(Frame gridpml, Margin margin,
                           int it, int ntime, int ReceNum, Point3D *rece,
                           float *vx, float *vy, float *vz,
                           float *vx_seis, float *vy_seis, float *vz_seis)
    {
        set_frame_ndl(gridpml);
        int iRece = threadIdx.x;
        int idx_seis = it + iRece * ntime;
        int Rx, Ry, Rz, idx_rece;
        if (iRece < ReceNum)
        {
            Rx = get_gap_num(l_rows, rece[iRece].x, d_rows);
            Ry = get_gap_num(l_cols, rece[iRece].y, d_cols);
            Rz = get_gap_num(l_slices, rece[iRece].z, d_slices);
            idx_rece = Rx + Ry * n_rows + Rz * n_elem_slice;
            vx_seis[idx_seis] = vx[idx_rece];
            vy_seis[idx_seis] = vy[idx_rece];
            vz_seis[idx_seis] = vz[idx_rece];
        }
    }

    inline __global__ void
    sub_rtm_reverse_add_vel_source(Frame sub_gridext,
                                   int it, int ntime, int ReceNum, Point3D *rece,
                                   float *vx, float *vy, float *vz,
                                   float *vx_seis, float *vy_seis, float *vz_seis)
    {
        set_frame_nd(sub_gridext);
        int iRece = threadIdx.x + blockDim.x * blockIdx.x;
        int iRx, iRy, iRz, idx_rece, idx_seis;
        double l_rows, l_cols, l_slices, r_rows, r_cols, r_slices;
        l_rows = sub_gridext.l_rows;
        l_cols = sub_gridext.l_cols;
        l_slices = sub_gridext.l_slices;
        r_rows = sub_gridext.r_rows;
        r_cols = sub_gridext.r_cols;
        r_slices = sub_gridext.r_slices;
        //
        if (iRece < ReceNum)
        {
            iRx = get_gap_num(l_rows, rece[iRece].x, d_rows);
            iRy = get_gap_num(l_cols, rece[iRece].y, d_cols);
            iRz = get_gap_num(l_slices, rece[iRece].z, d_slices);
            idx_seis = it + iRece * ntime;
            idx_rece = iRx + iRy * n_rows + iRz * n_elem_slice;
            vx[idx_rece] += vx_seis[idx_seis];
            vy[idx_rece] += vy_seis[idx_seis];
            vz[idx_rece] += vz_seis[idx_seis];
        }
    }

    inline __global__ void
    sub_add_reverse_vel_source(SimulateType simulate_type, Frame sub_gridext,
                               int it, int ntime, int ReceNum, Point3D *rece,
                               float *vx, float *vy, float *vz,
                               float *vx_seis, float *vy_seis, float *vz_seis,
                               float *vx_fd, float *vy_fd, float *vz_fd)
    {
        set_frame_nd(sub_gridext);
        int iRece = threadIdx.x + blockDim.x * blockIdx.x;
        int iRx, iRy, iRz, idx_rece, idx_seis;
        double l_rows, l_cols, l_slices, r_rows, r_cols, r_slices;
        l_rows = sub_gridext.l_rows;
        l_cols = sub_gridext.l_cols;
        l_slices = sub_gridext.l_slices;
        r_rows = sub_gridext.r_rows;
        r_cols = sub_gridext.r_cols;
        r_slices = sub_gridext.r_slices;
        //
        if (iRece < ReceNum)
        {
            iRx = get_gap_num(l_rows, rece[iRece].x, d_rows);
            iRy = get_gap_num(l_cols, rece[iRece].y, d_cols);
            iRz = get_gap_num(l_slices, rece[iRece].z, d_slices);
            idx_seis = it + iRece * ntime;
            idx_rece = iRx + iRy * n_rows + iRz * n_elem_slice;
            if (simulate_type == SimulateType::rtm_reverse)
            {
                vx[idx_rece] += vx_seis[idx_seis];
                vy[idx_rece] += vy_seis[idx_seis];
                vz[idx_rece] += vz_seis[idx_seis];
            }
            if (simulate_type == SimulateType::fwi)
            {
                vx[idx_rece] += (vx_seis[idx_seis] - vx_fd[idx_seis]);
                vy[idx_rece] += (vy_seis[idx_seis] - vy_fd[idx_seis]);
                vz[idx_rece] += (vz_seis[idx_seis] - vz_fd[idx_seis]);
            }
        }
    }
    //
    inline __global__ void mute_first_wave(int iShot, int ReceNum, Point3D *shot, Point3D *rece,
                                           float v_first, float dt, float fm, int ntime,
                                           float *vx_seis, float *vy_seis, float *vz_seis)
    {
        int i = threadIdx.x + blockDim.x * blockIdx.x;
        if (i < ReceNum)
        {
            if (abs(shot[iShot].y - rece[i].y) < 20)
            {
                float t = abs(rece[i].x - shot[iShot].x) / v_first + 0.04;
                int idx_l = int(t / dt);
                int idx_r = int(t / dt + 1.0 / (fm * dt)) + 200;
                if (idx_l > ntime - 1)
                    idx_l = ntime - 1;

                if (idx_r > ntime - 1)
                    idx_r = ntime - 1;
                for (int it = 0; it <= idx_r; it++)
                {
                    vx_seis[it + i * ntime] = 0;
                    vy_seis[it + i * ntime] = 0;
                    vz_seis[it + i * ntime] = 0;
                }
            }
        }
    }
    inline __global__ void cal_poynting_vec(Frame gridext, Frame gridphy, Margin sub_margin,
                                            float *vx, float *vy, float *vz,
                                            float *sxx, float *sxy, float *sxz,
                                            float *syy, float *syz, float *szz, uvec3 *pyt_vec)
    {
        set_cufield_3d_idx(gridphy, _idx, _i, _j, _k);
        if (_idx < gridphy.n_elem)
        {
            int ext_i = _i + sub_margin.right_margin;
            int ext_j = _j + sub_margin.front_margin;
            int ext_k = _k + sub_margin.top_margin;
            int ext_idx = ext_i + ext_j * gridext.n_rows + ext_k * gridext.n_elem_slice;

            double pyt1 = -(sxx[ext_idx] * vx[ext_idx] +
                            sxy[ext_idx] * vy[ext_idx] +
                            sxz[ext_idx] * vz[ext_idx]);

            double pyt2 = -(sxy[ext_idx] * vx[ext_idx] +
                            syy[ext_idx] * vy[ext_idx] +
                            syz[ext_idx] * vz[ext_idx]);

            double pyt3 = -(sxz[ext_idx] * vx[ext_idx] +
                            syz[ext_idx] * vy[ext_idx] +
                            szz[ext_idx] * vz[ext_idx]);

            double pyt_L = sqrt(pyt1 * pyt1 + pyt2 * pyt2 + pyt3 * pyt3);

            if (pyt_L < 1e-60) //*Prevent pyt_L too small to avoid division by zero
            {
                pyt_vec[_idx].x = 0;
                pyt_vec[_idx].y = 0;
                pyt_vec[_idx].z = 0;
            }
            else
            {
                pyt_vec[_idx].x = pyt1 / pyt_L;
                pyt_vec[_idx].y = pyt2 / pyt_L;
                pyt_vec[_idx].z = pyt3 / pyt_L;
            }

            // F_(pyt_vec, i, j, k).x = pyt1;
            // F_(pyt_vec, i, j, k).y = pyt2;
            // F_(pyt_vec, i, j, k).z = pyt3;
        }
    }
    //
    inline __global__ void
    get_rece_pyt(Frame gridext, Margin sub_margin,
                 int it, int ntime, int ReceNum, Point3D *rece,
                 float *vx, float *vy, float *vz,
                 float *sxx, float *sxy, float *sxz,
                 float *syy, float *syz, float *szz, uvec3 *pyt_vec_rece)
    {
        set_cufield_3d_idx(gridext, idx, i, j, k);
        int iRece = threadIdx.x;
        int idx_seis = it + iRece * ntime;
        int Rx, Ry, Rz, idx_rece;
        double l_rows, l_cols, l_slices, r_rows, r_cols, r_slices;
        l_rows = gridext.l_rows;
        l_cols = gridext.l_cols;
        l_slices = gridext.l_slices;
        r_rows = gridext.r_rows;
        r_cols = gridext.r_cols;
        r_slices = gridext.r_slices;
        if (iRece < ReceNum)
        {
            Rx = get_gap_num(l_rows, rece[iRece].x, gridext.d_rows);
            Ry = get_gap_num(l_cols, rece[iRece].y, gridext.d_cols);
            Rz = get_gap_num(l_slices, rece[iRece].z, gridext.d_slices);
            idx_rece = Rx + Ry * n_rows + Rz * n_elem_slice;
            double pyt1 = -(sxx[idx_rece] * vx[idx_rece] +
                            sxy[idx_rece] * vy[idx_rece] +
                            sxz[idx_rece] * vz[idx_rece]);

            double pyt2 = -(sxy[idx_rece] * vx[idx_rece] +
                            syy[idx_rece] * vy[idx_rece] +
                            syz[idx_rece] * vz[idx_rece]);

            double pyt3 = -(sxz[idx_rece] * vx[idx_rece] +
                            syz[idx_rece] * vy[idx_rece] +
                            szz[idx_rece] * vz[idx_rece]);

            double pyt_L = sqrt(pyt1 * pyt1 + pyt2 * pyt2 + pyt3 * pyt3);

            if (pyt_L < 1e-60) //*Prevent pyt_L too small to avoid division by zero
            {
                pyt_vec_rece[idx_seis].x = 0;
                pyt_vec_rece[idx_seis].y = 0;
                pyt_vec_rece[idx_seis].z = 0;
            }
            else
            {
                pyt_vec_rece[idx_seis].x = pyt1 / pyt_L;
                pyt_vec_rece[idx_seis].y = pyt2 / pyt_L;
                pyt_vec_rece[idx_seis].z = pyt3 / pyt_L;
            }
        }
    }
    inline void SeisRecord::initialize()
    {
        fm = arg_list_p->fm;
        dt = arg_list_p->dt;
        T = arg_list_p->T;
        ntime = ceil(T / dt);
        float wavelet_dura = 1.3f / fm;
        int wavelet_nodes = (int)floor(wavelet_dura / dt + 0.5f);
        wave_let.alloc(ntime);
        for (int i = 0; i < ntime; i++)
        {
            wave_let(i) = 1.0e4f * (1.0f - 2.0f * pow(jarvis_const::pi * fm * ((i - wavelet_nodes) * dt), 2)) * exp((-1) * pow(jarvis_const::pi * fm * ((i - wavelet_nodes) * dt), 2));
        }
        wave_let.copy_h2d();
        //
        if (sub_geometry_p->is_have_recv)
        {
            vx_seis.alloc(Frame(ntime, sub_geometry_p->sub_recv.n_elem));
            vy_seis.alloc(Frame(ntime, sub_geometry_p->sub_recv.n_elem));
            vz_seis.alloc(Frame(ntime, sub_geometry_p->sub_recv.n_elem));
            vxp_seis.alloc(Frame(ntime, sub_geometry_p->sub_recv.n_elem));
            vyp_seis.alloc(Frame(ntime, sub_geometry_p->sub_recv.n_elem));
            vzp_seis.alloc(Frame(ntime, sub_geometry_p->sub_recv.n_elem));
        }
    }
    //
    inline void SeisRecord::alloc_host_storage()
    {
        if (sub_geometry_p->is_have_recv)
        {
            vx_host_vec.alloc(sub_geometry_p->shot.n_elem);
            vy_host_vec.alloc(sub_geometry_p->shot.n_elem);
            vz_host_vec.alloc(sub_geometry_p->shot.n_elem);
            for (int ishot = 0; ishot < sub_geometry_p->shot.n_elem; ishot++)
            {
                vx_host_vec(ishot).alloc(Frame(ntime, sub_geometry_p->sub_recv.n_elem));
                vy_host_vec(ishot).alloc(Frame(ntime, sub_geometry_p->sub_recv.n_elem));
                vz_host_vec(ishot).alloc(Frame(ntime, sub_geometry_p->sub_recv.n_elem));
            }
        }
    }
    //
    inline void SeisRecord::copy_seis_into_host(int _i_shot)
    {
        if (sub_geometry_p->is_have_recv)
        {
            vx_seis.copy_into_host(vx_host_vec(_i_shot), jarvis_mpi_cuda_stream_p->d2h_stream());
            vy_seis.copy_into_host(vy_host_vec(_i_shot), jarvis_mpi_cuda_stream_p->d2h_stream());
            vz_seis.copy_into_host(vz_host_vec(_i_shot), jarvis_mpi_cuda_stream_p->d2h_stream());
            jarvis_mpi_cuda_stream_p->sync_d2h_stream();
        }
    }
    //
    inline void SeisRecord::copy_seis_from_host(int _i_shot)
    {
        if (sub_geometry_p->is_have_recv)
        {
            vx_seis.copy_from_host(vx_host_vec(_i_shot), jarvis_mpi_cuda_stream_p->h2d_stream());
            vy_seis.copy_from_host(vy_host_vec(_i_shot), jarvis_mpi_cuda_stream_p->h2d_stream());
            vz_seis.copy_from_host(vz_host_vec(_i_shot), jarvis_mpi_cuda_stream_p->h2d_stream());
            jarvis_mpi_cuda_stream_p->sync_h2d_stream();
        }
    }

    inline void SeisRecord::cuda_get_seis_data(int it)
    {
        if (sub_geometry_p->is_have_recv)
        {
            void *args_list[] = {&geomodel_p->gridext,
                                 &geomodel_p->margin,
                                 &it,
                                 &ntime,
                                 &sub_geometry_p->sub_recv.n_elem,
                                 &sub_geometry_p->sub_recv.vector<Point3D, MemType::device>::ptr(),
                                 &ext_wavefield_p->vx.ptr(),
                                 &ext_wavefield_p->vy.ptr(),
                                 &ext_wavefield_p->vz.ptr(),
                                 &vx_seis.ptr(),
                                 &vy_seis.ptr(),
                                 &vz_seis.ptr()};
            cudaLaunchKernel((void *)sub_cal_seis_data_base, 1, sub_geometry_p->sub_recv.n_elem, args_list, 0, jarvis_mpi_cuda_stream_p->cal_stream());
        }
    }

    // //
    // inline __global__ void
    // fliter_rece_pyt(Position pos_flag, int ntime, float *seis, uvec3 *pyt_vec_rece)
    // {
    //     int idx = JARVIS_CUDA_IDX_X;
    //     int it = idx % ntime;
    //     int iRece = idx / ntime;
    //     float sum = 0;
    //     int width = 200;
    //     if (pos_flag == Position::top)
    //     {
    //         if (it > width && it < ntime - width)
    //         {
    //             for (int i = idx - width; i < idx + width; i++)
    //             {
    //                 sum += pyt_vec_rece[i].z;
    //             }
    //         }
    //         if (sum >= 0)
    //             seis[idx] = 0;
    //     }
    //     else if (pos_flag == Position::bottom)
    //     {
    //         if (it > width && it < ntime - width)
    //         {
    //             for (int i = idx - width; i < idx + width; i++)
    //             {
    //                 sum += pyt_vec_rece[i].z;
    //             }
    //         }
    //         if (sum < 0)
    //             seis[idx] = 0;
    //     }
    //     else if (pos_flag == Position::left)
    //     {
    //         if (it > width && it < ntime - width)
    //         {
    //             for (int i = idx - width; i < idx + width; i++)
    //             {
    //                 sum += pyt_vec_rece[i].x;
    //             }
    //         }
    //         if (sum < 0)
    //             seis[idx] = 0;
    //     }
    //     else if (pos_flag == Position::right)
    //     {
    //         if (it > width && it < ntime - width)
    //         {
    //             for (int i = idx - width; i < idx + width; i++)
    //             {
    //                 sum += pyt_vec_rece[i].x;
    //             }
    //         }
    //         if (sum >= 0)
    //             seis[idx] = 0;
    //     }
    //     else if (pos_flag == Position::front)
    //     {
    //         if (it > width && it < ntime - width)
    //         {
    //             for (int i = idx - width; i < idx + width; i++)
    //             {
    //                 sum += pyt_vec_rece[i].y;
    //             }
    //         }
    //         if (sum >= 0)
    //             seis[idx] = 0;
    //     }
    //     else if (pos_flag == Position::back)
    //     {
    //         if (it > width && it < ntime - width)
    //         {
    //             for (int i = idx - width; i < idx + width; i++)
    //             {
    //                 sum += pyt_vec_rece[i].y;
    //             }
    //         }
    //         if (sum < 0)
    //             seis[idx] = 0;
    //     }
    // }
}
#endif