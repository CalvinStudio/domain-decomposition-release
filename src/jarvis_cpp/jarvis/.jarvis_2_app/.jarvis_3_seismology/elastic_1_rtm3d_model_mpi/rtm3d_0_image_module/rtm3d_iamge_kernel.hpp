#ifndef _RTM3D_BASE_MODULE_KERNEL_HPP
#define _RTM3D_BASE_MODULE_KERNEL_HPP
#include "rtm3d_image_bones.hpp"
namespace jarvis
{
    inline __global__ void
    cal_correlation_from_wave(Frame sub_gridext, Frame sub_gridphy,
                              float *vx_fd, float *vy_fd, float *vz_fd,
                              float *vxp_fd, float *vyp_fd, float *vzp_fd,
                              float *vx, float *vy, float *vz,
                              float *vxp, float *vyp, float *vzp,
                              float *Ip_fd,
                              float *image_pp, float *image_ps, float *image_sp, float *image_ss)
    {
        set_cufield_grid_3d_idx(sub_gridext, sub_gridphy);
        if (idx < sub_gridphy.n_elem)
        {
            float vxp_fd_tmp = vxp_fd[ext_idx];
            float vyp_fd_tmp = vyp_fd[ext_idx];
            float vzp_fd_tmp = vzp_fd[ext_idx];
            float vxs_fd_tmp = vx_fd[ext_idx] - vxp_fd_tmp;
            float vys_fd_tmp = vy_fd[ext_idx] - vyp_fd_tmp;
            float vzs_fd_tmp = vz_fd[ext_idx] - vzp_fd_tmp;
            float vxp_tmp = vxp[ext_idx];
            float vyp_tmp = vyp[ext_idx];
            float vzp_tmp = vzp[ext_idx];
            float vxs_tmp = vx[ext_idx] - vxp_tmp;
            float vys_tmp = vy[ext_idx] - vyp_tmp;
            float vzs_tmp = vz[ext_idx] - vzp_tmp;

            Ip_fd[idx] += vxp_fd_tmp * vxp_fd_tmp +
                          vyp_fd_tmp * vyp_fd_tmp +
                          vzp_fd_tmp * vzp_fd_tmp;

            image_pp[idx] += vxp_fd_tmp * vxp_tmp +
                             vyp_fd_tmp * vyp_tmp +
                             vzp_fd_tmp * vzp_tmp;

            image_ps[idx] += vxp_fd_tmp * vxs_tmp +
                             vyp_fd_tmp * vys_tmp +
                             vzp_fd_tmp * vzs_tmp;

            image_sp[idx] += vxs_fd_tmp * vxp_tmp +
                             vys_fd_tmp * vyp_tmp +
                             vzs_fd_tmp * vzp_tmp;

            image_ss[idx] += vxs_fd_tmp * vxs_tmp +
                             vys_fd_tmp * vys_tmp +
                             vzs_fd_tmp * vzs_tmp;
        }
    }
    //
    inline __global__ void
    cal_correlation_from_tmp_image(Frame sub_gridphy,
                                   float *Ip_fd_tmp,
                                   float *Ipp_tmp, float *Ips_tmp, float *Isp_tmp, float *Iss_tmp,
                                   float *Ip_fd,
                                   float *image_pp, float *image_ps, float *image_sp, float *image_ss)
    {
        set_cufield_3d_idx(sub_gridphy, idx, _i, _j, _k);
        if (idx < sub_gridphy.n_elem)
        {
            Ip_fd[idx] += Ip_fd_tmp[idx];
            image_pp[idx] += Ipp_tmp[idx];
            image_ps[idx] += Ips_tmp[idx];
            image_sp[idx] += Isp_tmp[idx];
            image_ss[idx] += Iss_tmp[idx];
        }
    }
    //
    inline __global__ void
    cal_correlation_tmp(Frame sub_gridext, Frame sub_gridphy,
                        Margin margin,
                        float *vx, float *vy, float *vz,
                        float *vxp, float *vyp, float *vzp,
                        float *vx_fd, float *vy_fd, float *vz_fd,
                        float *vxp_fd, float *vyp_fd, float *vzp_fd,
                        float *Ip_fd_tmp,
                        float *Ipp_tmp, float *Ips_tmp, float *Isp_tmp, float *Iss_tmp)
    {
        set_cufield_3d_idx(sub_gridphy, idx, _i, _j, _k);
        if (idx < sub_gridphy.n_elem)
        {
            int ext_i = _i + margin.right_margin;
            int ext_j = _j + margin.front_margin;
            int ext_k = _k + margin.top_margin;
            int ext_idx = ext_i + ext_j * sub_gridext.n_rows + ext_k * sub_gridext.n_elem_slice;

            float vxs = vx[ext_idx] - vxp[ext_idx];
            float vys = vy[ext_idx] - vyp[ext_idx];
            float vzs = vz[ext_idx] - vzp[ext_idx];
            float vxs_forward = vx_fd[ext_idx] - vxp_fd[ext_idx];
            float vys_forward = vy_fd[ext_idx] - vyp_fd[ext_idx];
            float vzs_forward = vz_fd[ext_idx] - vzp_fd[ext_idx];

            Ip_fd_tmp[idx] = vxp_fd[ext_idx] * vxp_fd[ext_idx] +
                             vyp_fd[ext_idx] * vyp_fd[ext_idx] +
                             vzp_fd[ext_idx] * vzp_fd[ext_idx];

            Ipp_tmp[idx] = vxp_fd[ext_idx] * vxp[ext_idx] +
                           vyp_fd[ext_idx] * vyp[ext_idx] +
                           vzp_fd[ext_idx] * vzp[ext_idx];

            Ips_tmp[idx] = vxp_fd[ext_idx] * vxs +
                           vyp_fd[ext_idx] * vys +
                           vzp_fd[ext_idx] * vzs;

            Isp_tmp[idx] = vxs_forward * vxp[ext_idx] +
                           vys_forward * vyp[ext_idx] +
                           vzs_forward * vzp[ext_idx];

            Iss_tmp[idx] = vxs_forward * vxs +
                           vys_forward * vys +
                           vzs_forward * vzs;
        }
    }
    //
    inline __global__ void
    source_compensation(Frame sub_gridphy,
                        float *Ip_fd,
                        float *image_pp, float *image_ps, float *image_sp, float *image_ss)
    {
        set_cufield_3d_idx(sub_gridphy, idx, _i, _j, _k);
        if (idx < sub_gridphy.n_elem)
        {
            float Ip_fd_tmp = Ip_fd[idx];
            if (Ip_fd_tmp > 3.4E-38)
            {
                image_pp[idx] = image_pp[idx] / Ip_fd_tmp;
                image_ps[idx] = image_ps[idx] / Ip_fd_tmp;
                image_sp[idx] = image_sp[idx] / Ip_fd_tmp;
                image_ss[idx] = image_ss[idx] / Ip_fd_tmp;
            }
        }
    }
    //
    inline __global__ void
    cal_poynting_vec_angle(Frame gridphy, uvec3 *p1, uvec3 *p2, float *angle)
    {
        set_cufield_3d_idx(gridphy, idx, i, j, k);
        float angle_tmp;
        if (idx < gridphy.n_elem)
        {
            angle_tmp = (p1[idx].x * p2[idx].x +
                         p1[idx].y * p2[idx].y +
                         p1[idx].z * p2[idx].z);
            if (angle_tmp > 1)
                angle_tmp = 1;
            if (angle_tmp < -1)
                angle_tmp = -1;

            angle[idx] = acos(angle_tmp) * 90.0 / jarvis_const::pi;
        }
    }
    //
    inline __global__ void
    cal_adcig_image_one_angle(Frame sub_gridext, Frame sub_gridphy,
                              Margin margin,
                              float *_pyt_angle, float _angle,
                              float *Ipp_tmp, float *Ips_tmp, float *Isp_tmp, float *Iss_tmp,
                              float *Ipp_angle, float *Ips_angle, float *Isp_angle, float *Iss_angle)
    {
        set_cufield_3d_idx(sub_gridphy, idx, _i, _j, _k);
        if (idx < sub_gridphy.n_elem)
        {
            int ext_i = _i + margin.right_margin;
            int ext_j = _j + margin.front_margin;
            int ext_k = _k + margin.top_margin;
            int ext_idx = ext_i + ext_j * sub_gridext.n_rows + ext_k * sub_gridext.n_elem_slice;

            float delta = 2;
            float angle_weight = exp(-((_pyt_angle[idx] - _angle) * (_pyt_angle[idx] - _angle) / (2.0 * delta * delta)));

            Ipp_angle[idx] += Ipp_tmp[idx] * angle_weight;

            Ips_angle[idx] += Ips_tmp[idx] * angle_weight;

            Isp_angle[idx] += Isp_tmp[idx] * angle_weight;

            Iss_angle[idx] += Iss_tmp[idx] * angle_weight;
        }
    }

    inline __global__ void
    laplace_filter_cal_for_image(Frame gridphy,
                                 float *image_pp, float *image_ps, float *image_sp, float *image_ss,
                                 float *Ipp_laplace, float *Ips_laplace,
                                 float *Isp_laplace, float *Iss_laplace)
    {
        set_cufield_3d_idx(gridphy, idx, i, j, k);
        if (i >= 2 && i <= n_rows - 3)
            if (j >= 2 && j <= gridphy.n_cols - 3)
                if (k >= 2 && k <= gridphy.n_slices - 3)
                {
                    Ipp_laplace[idx] = image_pp[idx + 1] + image_pp[idx - 1] + image_pp[idx + n_rows] + image_pp[idx - n_rows] + image_pp[idx + n_elem_slice] + image_pp[idx - n_elem_slice] - 6.f * image_pp[idx];
                    Ips_laplace[idx] = image_ps[idx + 1] + image_ps[idx - 1] + image_ps[idx + n_rows] + image_ps[idx - n_rows] + image_ps[idx + n_elem_slice] + image_ps[idx - n_elem_slice] - 6.f * image_ps[idx];
                    Isp_laplace[idx] = image_sp[idx + 1] + image_sp[idx - 1] + image_sp[idx + n_rows] + image_sp[idx - n_rows] + image_sp[idx + n_elem_slice] + image_sp[idx - n_elem_slice] - 6.f * image_sp[idx];
                    Iss_laplace[idx] = image_ss[idx + 1] + image_ss[idx - 1] + image_ss[idx + n_rows] + image_ss[idx - n_rows] + image_ss[idx + n_elem_slice] + image_ss[idx - n_elem_slice] - 6.f * image_ss[idx];
                }
    }
}
#endif