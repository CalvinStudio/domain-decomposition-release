#pragma once
#ifndef _FDTD3D_GEOMODEL_MEAT_HPP
#define _FDTD3D_GEOMODEL_MEAT_HPP
#include "fdtd3d_base_3_geomodel_bones.hpp"
namespace jarvis
{
    inline __global__ void smooth_model_x_dir(Frame grid, int x_smooth_size, float *model_data_org, float *model_data_smo)
    {
        set_cufield_3d_idx(grid, idx, i, j, k);
        int lbx, rbx;
        int is;
        if (i < grid.n_rows && j < grid.n_cols && k < grid.n_slices)
        {
            lbx = -x_smooth_size;
            rbx = x_smooth_size;
            //
            if (i < x_smooth_size)
            {
                lbx = -i;
            }
            else if (i >= grid.n_rows - x_smooth_size)
            {
                rbx = grid.n_rows - i - 1;
            }
            double Sum = 0;
            int sm_grid_size = (rbx - lbx + 1);
            for (is = lbx; is <= rbx; is++)
            {
                Sum += model_data_org[idx + is];
            }
            model_data_smo[idx] = Sum / sm_grid_size;
        }
    }
    //

    inline __global__ void smooth_model_y_dir(Frame grid,
                                              int y_smooth_size,
                                              float *model_data_org, float *model_data_smo)
    {
        set_cufield_3d_idx(grid, idx, i, j, k);
        int lby, rby;
        int js;
        if (i < grid.n_rows && j < grid.n_cols && k < grid.n_slices)
        {
            lby = -y_smooth_size;
            rby = y_smooth_size;
            //
            if (j < y_smooth_size)
            {
                lby = -j;
            }
            else if (j >= grid.n_cols - y_smooth_size)
            {
                rby = grid.n_cols - j - 1;
            }
            double Sum = 0;
            int sm_grid_size = (rby - lby + 1);
            for (js = lby; js <= rby; js++)
            {
                Sum += model_data_org[idx + n_rows];
            }
            model_data_smo[idx] = Sum / sm_grid_size;
        }
    }
    //

    inline __global__ void smooth_model_z_dir(Frame grid,
                                              int z_smooth_size,
                                              float *model_data_org, float *model_data_smo)
    {
        set_cufield_3d_idx(grid, idx, i, j, k);
        int lbz, rbz;
        int ks;
        if (i < grid.n_rows && j < grid.n_cols && k < grid.n_slices)
        {
            lbz = -z_smooth_size;
            rbz = z_smooth_size;
            //
            if (k < z_smooth_size)
                lbz = -k;
            else if (k >= grid.n_slices - z_smooth_size)
                rbz = grid.n_slices - k - 1;
            float Sum = 0;
            int sm_grid_size = (rbz - lbz + 1);
            for (ks = lbz; ks <= rbz; ks++)
            {
                Sum += model_data_org[idx + ks * n_elem_slice];
            }
            model_data_smo[idx] = Sum / sm_grid_size;
        }
    }

    inline __global__ void adjust_para(float adjust, float *lambda, float *mu, float *rho)
    {
        int idx = threadIdx.x + blockDim.x * blockIdx.x;
        mu[idx] = mu[idx] * (1.0f + adjust) * (1.0f + adjust);
        lambda[idx] = lambda[idx] * (1.0f + adjust) * (1.0f + adjust);
    }

    template <class MarginType>
    inline void extend_model(Frame gridext, MarginType margin, Field<float, MemType::paged_device> &_a_ext, Field<float> &_a_phy)
    {
        // gridext.print_info("gridext:");
        _a_ext.alloc(gridext);
        /* Conversion */
        for (int k = margin.top_margin; k < (gridext.n_slices - margin.bottom_margin); k++)
            for (int j = margin.front_margin; j < (gridext.n_cols - margin.back_margin); j++)
                for (int i = margin.right_margin; i < (gridext.n_rows - margin.left_margin); i++)
                {
                    _a_ext(i, j, k) = _a_phy(i - margin.right_margin, j - margin.front_margin, k - margin.top_margin);
                }
        /*
         *    Parameterize the boundary layers.
         *    26 regions.
         */
        int left_margin = margin.left_margin;
        int right_margin = margin.right_margin;
        int front_margin = margin.front_margin;
        int back_margin = margin.back_margin;
        int top_margin = margin.top_margin;
        int bottom_margin = margin.bottom_margin;
        set_frame_nd(gridext);
        /* start x-axial-left */
        /* start y-axial-back */

        // 1 -- left-upper-back cube, is the corner.
        for (int k = 0; k < top_margin; k++)
            for (int j = 0; j < front_margin; j++)
                for (int i = 0; i < right_margin; i++)
                {
                    _a_ext(i, j, k) = _a_ext(right_margin, front_margin, top_margin);
                }
        // cout << miu(top_margin, front_margin, right_margin) << endl;
        // cout << top_margin << front_margin << right_margin << endl;
        // 2 -- left-back cube, is the edge.
        for (int k = top_margin; k < (n_slices - bottom_margin); k++)
            for (int j = 0; j < front_margin; j++)
                for (int i = 0; i < right_margin; i++)
                {
                    _a_ext(i, j, k) = _a_ext(right_margin, front_margin, k);
                }

        // 3 -- left-down-back cube, is the corner.
        for (int k = (n_slices - bottom_margin); k < n_slices; k++)
            for (int j = 0; j < front_margin; j++)
                for (int i = 0; i < right_margin; i++)
                {
                    _a_ext(i, j, k) = _a_ext(right_margin, front_margin, n_slices - bottom_margin - 1);
                }

        /* end y-direction-back */

        /* start y-direction-middle */
        // 4 -- left-upper cube, is the edge.
        for (int k = 0; k < top_margin; k++)
            for (int j = front_margin; j < (n_cols - back_margin); j++)
                for (int i = 0; i < right_margin; i++)
                {
                    _a_ext(i, j, k) = _a_ext(right_margin, j, top_margin);
                }

        // 5 -- left-cube, is the plane.
        for (int k = top_margin; k < (n_slices - bottom_margin); k++)
            for (int j = front_margin; j < (n_cols - back_margin); j++)
                for (int i = 0; i < right_margin; i++)
                {
                    _a_ext(i, j, k) = _a_ext(right_margin, j, k);
                }

        // 6 -- left-down cube, is the edge.
        for (int k = n_slices - bottom_margin; k < n_slices; k++)
            for (int j = front_margin; j < (n_cols - back_margin); j++)
                for (int i = 0; i < right_margin; i++)
                {
                    _a_ext(i, j, k) = _a_ext(right_margin, j, n_slices - bottom_margin - 1);
                }

        /* end y-axial-middle */

        /* start: y-axial-front */
        // 7 -- left-upper-front cube, is the corner.
        for (int k = 0; k < top_margin; k++)
            for (int j = n_cols - back_margin; j < n_cols; j++)
                for (int i = 0; i < right_margin; i++)
                {
                    _a_ext(i, j, k) = _a_ext(right_margin, n_cols - back_margin - 1, top_margin);
                }

        // 8 -- left-front cube, is the edge.
        for (int k = top_margin; k < (n_slices - bottom_margin); k++)
            for (int j = (n_cols - back_margin); j < n_cols; j++)
                for (int i = 0; i < right_margin; i++)
                {
                    _a_ext(i, j, k) = _a_ext(right_margin, n_cols - back_margin - 1, k);
                }

        // 9 -- left-down-front cube, is the corner.
        for (int k = (n_slices - bottom_margin); k < n_slices; k++)
            for (int j = (n_cols - back_margin); j < n_cols; j++)
                for (int i = 0; i < right_margin; i++)
                {
                    _a_ext(i, j, k) = _a_ext(right_margin, n_cols - back_margin - 1, n_slices - bottom_margin - 1);
                }

        /* end y-axial-front */
        /* end x-axial-left */

        /* start: x-axial-middle */
        /* start: y-axial-back. */
        // 10 -- back-upper cube, is the edge.
        for (int k = 0; k < top_margin; k++)
            for (int j = 0; j < front_margin; j++)
                for (int i = right_margin; i < (n_rows - left_margin); i++)
                {
                    _a_ext(i, j, k) = _a_ext(i, front_margin, top_margin);
                }

        // 11 -- back cube, is the plane.
        for (int k = top_margin; k < (n_slices - bottom_margin); k++)
            for (int j = 0; j < front_margin; j++)
                for (int i = right_margin; i < (n_rows - left_margin); i++)
                {
                    _a_ext(i, j, k) = _a_ext(i, front_margin, k);
                }

        // 12 -- back-down cube, is the edge.
        for (int k = (n_slices - bottom_margin); k < n_slices; k++)
            for (int j = 0; j < front_margin; j++)
                for (int i = right_margin; i < (n_rows - left_margin); i++)
                {
                    _a_ext(i, j, k) = _a_ext(i, front_margin, n_slices - bottom_margin - 1);
                }

        /* end y-axial-back */

        /* start: y-axial-middle */
        // 13 -- upper cube, is the plane.
        for (int k = 0; k < top_margin; k++)
            for (int j = front_margin; j < (n_cols - back_margin); j++)
                for (int i = right_margin; i < (n_rows - left_margin); i++)
                {
                    _a_ext(i, j, k) = _a_ext(i, j, top_margin);
                }

        // Here! Skip the physical domain.

        // 14 -- down cube, is the plane.
        for (int k = (n_slices - bottom_margin); k < n_slices; k++)
            for (int j = front_margin; j < (n_cols - back_margin); j++)
                for (int i = right_margin; i < (n_rows - left_margin); i++)
                {
                    _a_ext(i, j, k) = _a_ext(i, j, n_slices - bottom_margin - 1);
                }

        /* end y-axial-middle. */

        /* start: y-axial-front. */
        // 15 -- upper-front cube, is the edge.
        for (int k = 0; k < top_margin; k++)
            for (int j = (n_cols - back_margin); j < n_cols; j++)
                for (int i = right_margin; i < (n_rows - left_margin); i++)
                {
                    _a_ext(i, j, k) = _a_ext(i, n_cols - back_margin - 1, top_margin);
                }

        // 16 -- front cube, is the plane.
        for (int k = top_margin; k < (n_slices - bottom_margin); k++)
            for (int j = (n_cols - back_margin); j < n_cols; j++)
                for (int i = right_margin; i < (n_rows - left_margin); i++)
                {
                    _a_ext(i, j, k) = _a_ext(i, n_cols - back_margin - 1, k);
                }

        // 17 -- down-front cube, is the edge.
        for (int k = (n_slices - bottom_margin); k < n_slices; k++)
            for (int j = (n_cols - back_margin); j < n_cols; j++)
                for (int i = right_margin; i < (n_rows - left_margin); i++)
                {
                    _a_ext(i, j, k) = _a_ext(i, n_cols - back_margin - 1, n_slices - bottom_margin - 1);
                }

        /* end y-axial-front. */
        /* end x-axial-middle */

        /* start: x-axial-right */
        /* start: y-axial-back */
        // 18 -- right-upper-back cube, is the corner.
        for (int k = 0; k < top_margin; k++)
            for (int j = 0; j < front_margin; j++)
                for (int i = (n_rows - left_margin); i < n_rows; i++)
                {
                    _a_ext(i, j, k) = _a_ext(n_rows - left_margin - 1, front_margin, top_margin);
                }

        // 19 -- right-back cube, is the edge.
        for (int k = top_margin; k < (n_slices - bottom_margin); k++)
            for (int j = 0; j < front_margin; j++)
                for (int i = (n_rows - left_margin); i < n_rows; i++)
                {
                    _a_ext(i, j, k) = _a_ext(n_rows - left_margin - 1, front_margin, k);
                }

        // 20 -- right-down-back cube, is the corner.
        for (int k = (n_slices - bottom_margin); k < n_slices; k++)
            for (int j = 0; j < front_margin; j++)
                for (int i = (n_rows - left_margin); i < n_rows; i++)
                {
                    _a_ext(i, j, k) = _a_ext(n_rows - left_margin - 1, front_margin, n_slices - bottom_margin - 1);
                }

        /* end y-axial-back. */

        /* start: y-axial-middle */
        // 21 -- right-upper cube, is the edge.
        for (int k = 0; k < top_margin; k++)
            for (int j = front_margin; j < (n_cols - back_margin); j++)
                for (int i = (n_rows - left_margin); i < n_rows; i++)
                {
                    _a_ext(i, j, k) = _a_ext(n_rows - left_margin - 1, j, top_margin);
                }

        // 22 -- right cube, is the plane.
        for (int k = top_margin; k < (n_slices - bottom_margin); k++)
            for (int j = front_margin; j < (n_cols - back_margin); j++)
                for (int i = (n_rows - left_margin); i < n_rows; i++)
                {
                    _a_ext(i, j, k) = _a_ext(n_rows - left_margin - 1, j, k);
                }

        // 23 -- right-down cube, is the edge.
        for (int k = (n_slices - bottom_margin); k < n_slices; k++)
            for (int j = front_margin; j < (n_cols - back_margin); j++)
                for (int i = (n_rows - left_margin); i < n_rows; i++)
                {
                    _a_ext(i, j, k) = _a_ext(n_rows - left_margin - 1, j, n_slices - bottom_margin - 1);
                }

        /* end y-axial-middle */

        /* start: y-axial-front */
        // 24 -- right-upper-front cube, is the corner.
        for (int k = 0; k < top_margin; k++)
            for (int j = (n_cols - back_margin); j < n_cols; j++)
                for (int i = (n_rows - left_margin); i < n_rows; i++)
                {
                    _a_ext(i, j, k) = _a_ext(n_rows - left_margin - 1, n_cols - back_margin - 1, top_margin);
                }

        // 25 -- right-front cube, is the edge.
        for (int k = top_margin; k < (n_slices - bottom_margin); k++)
            for (int j = (n_cols - back_margin); j < n_cols; j++)
                for (int i = (n_rows - left_margin); i < n_rows; i++)
                {
                    _a_ext(i, j, k) = _a_ext(n_rows - left_margin - 1, n_cols - back_margin - 1, k);
                }

        // 26 -- right-down-front cube, is the corner.
        for (int k = (n_slices - bottom_margin); k < n_slices; k++)
            for (int j = (n_cols - back_margin); j < n_cols; j++)
                for (int i = (n_rows - left_margin); i < n_rows; i++)
                {
                    _a_ext(i, j, k) = _a_ext(n_rows - left_margin - 1, n_cols - back_margin - 1, n_slices - bottom_margin - 1);
                }
    } // end function.

    inline void smooth_model3d_host(Field<float, MemType::paged_device> &model_data_org, Field<float, MemType::paged_device> &model_data_smo, int smooth_wsize)
    {
        set_frame_nd(model_data_org);
        int lbx, rbx, lby, rby, lbz, rbz;
        int is, js, ks;

        lbx = lby = lbz = -smooth_wsize;
        rbx = rby = rbz = smooth_wsize;
        for (int k = 0; k < n_slices; k++)
            for (int j = 0; j < n_cols; j++)
                for (int i = 0; i < n_rows; i++)
                {
                    if (i < smooth_wsize)
                        lbx = -i;
                    else if (i >= n_rows - smooth_wsize)
                        rbx = n_rows - i - 1;

                    if (j < smooth_wsize)
                        lby = -j;
                    else if (j >= n_cols - smooth_wsize)
                        rby = n_cols - j - 1;

                    if (k < smooth_wsize)
                        lbz = -k;
                    else if (k >= n_slices - smooth_wsize)
                        rbz = n_slices - k - 1;

                    float Sum = 0;
                    int sm_grid_size = (rbx - lbx + 1) * (rby - lby + 1) * (rbz - lbz + 1);
                    for (is = lbx; is <= rbx; is++)
                        for (js = lby; js <= rby; js++)
                            for (ks = lbz; ks <= rbz; ks++)
                            {
                                Sum += model_data_org(i + is, j + js, k + ks);
                            }
                    model_data_smo(i, j, k) = Sum / sm_grid_size;
                }
    }
    //
    inline void elasticGeoModel<domain::global>::initialize()
    {
        model_path = arg_list_p->model_path;
        margin.top_margin = arg_list_p->pml_num;
        margin.bottom_margin = arg_list_p->pml_num;
        margin.left_margin = arg_list_p->pml_num;
        margin.right_margin = arg_list_p->pml_num;
        margin.front_margin = arg_list_p->pml_num;
        margin.back_margin = arg_list_p->pml_num;

        GmsReader gms(model_path);
        gms.read_gms_by_order_to_ffield(phy_vp);
        gms.read_gms_by_order_to_ffield(phy_vs);
        gms.read_gms_by_order_to_ffield(phy_rho);
        gms.clear();

        // phy_vp.alloc(Frame(800, 800, 800, 1, 1, 1, 0, 0, 0));
        // phy_vs.alloc(Frame(800, 800, 800, 1, 1, 1, 0, 0, 0));
        // phy_rho.alloc(Frame(800, 800, 800, 1, 1, 1, 0, 0, 0));
        // phy_vp.fill(3000, 0, phy_vp.size() - 1);
        // phy_vs.fill(2000, 0, phy_vs.size() - 1);
        // phy_rho.fill(2400, 0, phy_rho.size() - 1);

        // for (int k = 0; k < 200; k++)
        //     for (int j = 0; j < 800; j++)
        //         for (int i = 0; i < 800; i++)
        //         {
        //             phy_vp(i, j, k) = 2000;
        //             phy_vs(i, j, k) = 1400;
        //             phy_rho(i, j, k) = 1400;
        //         }

        gridphy = phy_vp;

        gridext.set_ndl(gridphy.n_rows + margin.left_margin + margin.right_margin,
                        gridphy.n_cols + margin.front_margin + margin.back_margin,
                        gridphy.n_slices + margin.top_margin + margin.bottom_margin,
                        gridphy.d_rows, gridphy.d_cols, gridphy.d_slices,
                        -margin.right_margin * gridphy.d_rows,
                        -margin.front_margin * gridphy.d_cols,
                        -margin.top_margin * gridphy.d_slices);
        //
        extend_model(gridext, margin, vp, phy_vp);
        extend_model(gridext, margin, vs, phy_vs);
        extend_model(gridext, margin, rho, phy_rho);

        lambda.alloc(gridext);
        mu.alloc(gridext);
        for (int i = 0; i < gridext.n_elem; i++)
        {
            lambda[i] = rho[i] * (vp[i] * vp[i] - 2.0f * vs[i] * vs[i]);
            mu[i] = rho[i] * (vs[i] * vs[i]);
        }
        //
        lambda.copy_h2d();
        mu.copy_h2d();
        rho.copy_h2d();
        cudaDeviceSynchronize();
        lambda_smo.alloc(gridext);
        mu_smo.alloc(gridext);
        rho_smo.alloc(gridext);

        lambda_inv.alloc(gridext);
        mu_inv.alloc(gridext);
        rho_inv.alloc(gridext);

        cuda_smooth_model(20, 4, 1);
        cudaDeviceSynchronize();
        for (int i = 0; i < gridext.n_elem; i++)
        {
            lambda_inv[i] = lambda[0];
            mu_inv[i] = mu[0];
            rho_inv[i] = rho[0];
        }
        lambda_inv.copy_h2d();
        mu_inv.copy_h2d();
        rho_inv.copy_h2d();
        lambda_smo.copy_d2h();
        mu_smo.copy_d2h();
        rho_smo.copy_d2h();
        cudaDeviceSynchronize();
    }
    inline void elasticGeoModel<domain::global>::cuda_smooth_model(int x_smo_size, int y_smo_size, int z_smo_size)
    {
        void *args_list0[] = {&gridext, &x_smo_size, &lambda.device::ptr(), &lambda_smo.device::ptr()};
        cudaLaunchKernel((void *)smooth_model_x_dir, jarvis_cuda_kernel_size(gridext.n_elem), args_list0, 0, jarvis_default_cuda_stream);
        //
        void *args_list1[] = {&gridext, &x_smo_size, &mu.device::ptr(), &mu_smo.device::ptr()};
        cudaLaunchKernel((void *)smooth_model_x_dir, jarvis_cuda_kernel_size(gridext.n_elem), args_list1, 0, jarvis_default_cuda_stream);
        //
        void *args_list2[] = {&gridext, &x_smo_size, &rho.device::ptr(), &rho_smo.device::ptr()};
        cudaLaunchKernel((void *)smooth_model_x_dir, jarvis_cuda_kernel_size(gridext.n_elem), args_list2, 0, jarvis_default_cuda_stream);
        //
        void *args_list3[] = {&gridext, &x_smo_size, &lambda.device::ptr(), &lambda.device::ptr()};
        cudaLaunchKernel((void *)smooth_model_y_dir, jarvis_cuda_kernel_size(gridext.n_elem), args_list3, 0, jarvis_default_cuda_stream);
        //
        void *args_list4[] = {&gridext, &x_smo_size, &mu.device::ptr(), &mu.device::ptr()};
        cudaLaunchKernel((void *)smooth_model_y_dir, jarvis_cuda_kernel_size(gridext.n_elem), args_list4, 0, jarvis_default_cuda_stream);
        //
        void *args_list5[] = {&gridext, &x_smo_size, &rho.device::ptr(), &rho.device::ptr()};
        cudaLaunchKernel((void *)smooth_model_y_dir, jarvis_cuda_kernel_size(gridext.n_elem), args_list5, 0, jarvis_default_cuda_stream);
    }

    inline void elasticGeoModel<domain::global>::clear()
    {
        phy_vp.clear();
        phy_vs.clear();
        phy_rho.clear();
        //
        vp.clear();
        vs.clear();
        lambda.clear();
        mu.clear();
        rho.clear();
        lambda_smo.clear();
        mu_smo.clear();
        rho_smo.clear();
        lambda_inv.clear();
        mu_inv.clear();
        rho_inv.clear();
    }
    //
    inline void elasticGeoModel<domain::local>::initialize(elasticGeoModel<domain::global> *_glb_model_p)
    {
        set_gridphy(*_glb_model_p);
        set_gridext(*_glb_model_p);
        glb_vp_min = _glb_model_p->phy_vp.min();
        glb_vs_min = _glb_model_p->phy_vs.min();
        glb_vp_max = _glb_model_p->phy_vp.max();
        glb_vs_max = _glb_model_p->phy_vs.max();
        //
        vp.alloc(gridext);
        lambda.alloc(gridext);
        mu.alloc(gridext);
        rho.alloc(gridext);
        //
        lambda_smo.alloc(gridext);
        mu_smo.alloc(gridext);
        rho_smo.alloc(gridext);
        //
        lambda_inv.alloc(gridext);
        mu_inv.alloc(gridext);
        rho_inv.alloc(gridext);
        //
        set_model_para(_glb_model_p->gridext, _glb_model_p->vp, vp);
        set_model_para(_glb_model_p->gridext, _glb_model_p->lambda, lambda);
        set_model_para(_glb_model_p->gridext, _glb_model_p->mu, mu);
        set_model_para(_glb_model_p->gridext, _glb_model_p->rho, rho);
        //
        set_model_para(_glb_model_p->gridext, _glb_model_p->lambda_smo, lambda_smo);
        set_model_para(_glb_model_p->gridext, _glb_model_p->mu_smo, mu_smo);
        set_model_para(_glb_model_p->gridext, _glb_model_p->rho_smo, rho_smo);
        //
        set_model_para(_glb_model_p->gridext, _glb_model_p->lambda_inv, lambda_inv);
        set_model_para(_glb_model_p->gridext, _glb_model_p->mu_inv, mu_inv);
        set_model_para(_glb_model_p->gridext, _glb_model_p->rho_inv, rho_inv);
        //
        lambda.copy_h2d();
        mu.copy_h2d();
        rho.copy_h2d();
        //
        lambda_smo.copy_h2d();
        mu_smo.copy_h2d();
        rho_smo.copy_h2d();
        //
        lambda_inv.copy_h2d();
        mu_inv.copy_h2d();
        rho_inv.copy_h2d();
        cudaDeviceSynchronize();
    }

    inline void elasticGeoModel<domain::local>::set_gridphy(elasticGeoModel<domain::global> &_glb_model)
    {
        mpiFrame &_model_rank = jarvis_mpi_cuda_stream_p->mpi_frame;
        if (_model_rank.n_rows == 1)
        {
            gridphy.n_rows = _glb_model.gridphy.n_rows;
        }
        else
        {
            if (_model_rank.i_rows == 0)
            {
                gridphy.n_rows = _glb_model.gridext.n_rows / _model_rank.n_rows - _glb_model.margin.right_margin;
            }
            else if (_model_rank.i_rows > 0 && _model_rank.i_rows < _model_rank.n_rows - 1)
            {
                gridphy.n_rows = _glb_model.gridext.n_rows / _model_rank.n_rows;
            }
            else if (_model_rank.i_rows == _model_rank.n_rows - 1)
            {
                gridphy.n_rows = _glb_model.gridext.n_rows / _model_rank.n_rows + _glb_model.gridext.n_rows % _model_rank.n_rows - _glb_model.margin.left_margin;
            }
        }
        //
        int width = 0;
        if (_model_rank.n_cols == 1)
        {
            gridphy.n_cols = _glb_model.gridphy.n_cols;
        }
        else if (_model_rank.n_cols == 2)
        {
            if (_model_rank.i_cols == 0)
            {
                gridphy.n_cols = _glb_model.gridext.n_cols / _model_rank.n_cols - _glb_model.margin.front_margin;
            }
            else if (_model_rank.i_cols == 1)
            {
                gridphy.n_cols = _glb_model.gridext.n_cols / _model_rank.n_cols + _glb_model.gridext.n_cols % _model_rank.n_cols - _glb_model.margin.back_margin;
            }
        }
        else if (_model_rank.n_cols > 2)
        {
            if (_model_rank.i_cols == 0)
            {
                gridphy.n_cols = (_glb_model.gridext.n_cols - 2 * width) / _model_rank.n_cols - (_glb_model.margin.front_margin - width);
            }
            if (_model_rank.i_cols > 0 && _model_rank.i_cols < _model_rank.n_cols - 1)
            {
                gridphy.n_cols = (_glb_model.gridext.n_cols - 2 * width) / _model_rank.n_cols;
            }
            else if (_model_rank.i_cols == _model_rank.n_cols - 1)
            {
                gridphy.n_cols = (_glb_model.gridext.n_cols - 2 * width) / _model_rank.n_cols + (_glb_model.gridext.n_cols - 2 * width) % _model_rank.n_cols - (_glb_model.margin.back_margin - width);
            }
        }
        //
        if (_model_rank.n_slices == 1)
        {
            gridphy.n_slices = _glb_model.gridphy.n_slices;
        }
        else
        {
            if (_model_rank.i_slices == 0)
            {
                gridphy.n_slices = _glb_model.gridext.n_slices / _model_rank.n_slices - _glb_model.margin.top_margin;
            }
            if (_model_rank.i_slices > 0 && _model_rank.i_slices < _model_rank.n_slices - 1)
            {
                gridphy.n_slices = _glb_model.gridext.n_slices / _model_rank.n_slices;
            }
            else if (_model_rank.i_slices == _model_rank.n_slices - 1)
            {
                gridphy.n_slices = _glb_model.gridext.n_slices / _model_rank.n_slices + _glb_model.gridext.n_slices % _model_rank.n_slices - _glb_model.margin.bottom_margin;
            }
        }
        if (_model_rank.i_rows == 0)
            gridphy.l_rows = 0;
        if (_model_rank.i_rows > 0)
            gridphy.l_rows = (_model_rank.i_rows * int((_glb_model.gridext.n_rows) / _model_rank.n_rows)) * _glb_model.gridphy.d_rows - _glb_model.margin.right_margin * _glb_model.gridphy.d_rows;
        //
        if (_model_rank.i_cols == 0)
            gridphy.l_cols = 0;
        if (_model_rank.i_cols > 0)
            gridphy.l_cols = (_model_rank.i_cols * int((_glb_model.gridext.n_cols - 2 * width) / _model_rank.n_cols)) * _glb_model.gridphy.d_cols - (_glb_model.margin.front_margin - width) * _glb_model.gridphy.d_cols;
        //
        if (_model_rank.i_slices == 0)
            gridphy.l_slices = 0;
        if (_model_rank.i_slices > 0)
            gridphy.l_slices = (_model_rank.i_slices * int(_glb_model.gridext.n_slices / _model_rank.n_slices)) * _glb_model.gridphy.d_slices - _glb_model.margin.top_margin * _glb_model.gridphy.d_slices;
        //
        gridphy.d_rows = _glb_model.gridphy.d_rows;
        gridphy.d_cols = _glb_model.gridphy.d_cols;
        gridphy.d_slices = _glb_model.gridphy.d_slices;
        //
        gridphy.re_setndl();
        gridphy.print_info("mpi_rank==" + to_string(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank));
    }

    inline void elasticGeoModel<domain::local>::set_gridext(elasticGeoModel<domain::global> &_glb_model)
    {
        mpiFrame &_model_rank = jarvis_mpi_cuda_stream_p->mpi_frame;
        if (!_model_rank.near.is_top)
            margin.top_margin = _glb_model.margin.top_margin;

        if (!_model_rank.near.is_bottom)
            margin.bottom_margin = _glb_model.margin.bottom_margin;

        if (!_model_rank.near.is_front)
            margin.front_margin = _glb_model.margin.front_margin;

        if (!_model_rank.near.is_back)
            margin.back_margin = _glb_model.margin.back_margin;

        if (!_model_rank.near.is_right)
            margin.right_margin = _glb_model.margin.right_margin;

        if (!_model_rank.near.is_left)
            margin.left_margin = _glb_model.margin.left_margin;

        gridext.set_ndl(gridphy.n_rows + margin.left_margin + margin.right_margin,
                        gridphy.n_cols + margin.front_margin + margin.back_margin,
                        gridphy.n_slices + margin.top_margin + margin.bottom_margin,
                        gridphy.d_rows, gridphy.d_cols, gridphy.d_slices,
                        -margin.right_margin * gridphy.d_rows + gridphy.l_rows,
                        -margin.front_margin * gridphy.d_cols + gridphy.l_cols,
                        -margin.top_margin * gridphy.d_slices + gridphy.l_slices);
    }

    inline void elasticGeoModel<domain::local>::set_model_para(const Frame &_glb_gridext, Field<float, MemType::paged_device> &_glb_model_gms_para, Field<float, MemType::paged_device> &_sub_model)
    {
        int ix_l = get_gap_num(_glb_gridext.l_rows, gridext.l_rows, gridext.d_rows);
        int iy_l = get_gap_num(_glb_gridext.l_cols, gridext.l_cols, gridext.d_cols);
        int iz_l = get_gap_num(_glb_gridext.l_slices, gridext.l_slices, gridext.d_slices);
        frame_for(gridext, i, j, k)
        {
            int gi = i + ix_l;
            int gj = j + iy_l;
            int gk = k + iz_l;
            _sub_model(i, j, k) = _glb_model_gms_para(gi, gj, gk);
        }
    }
    //
    inline void elasticGeoModel<domain::local>::adjust_model_para(float adjust)
    {
        using device = jarvis::vector<float, MemType::device>;
        void *args_list[] = {&adjust,
                             &lambda.device::ptr(),
                             &mu.device::ptr(),
                             &rho.device::ptr()};
        cudaLaunchKernel((void *)adjust_para, jarvis_cuda_kernel_size(gridext.n_elem), args_list, 0, jarvis_default_cuda_stream);
    }

    inline void elasticGeoModel<domain::local>::host_clear()
    {
        lambda.host::clear();
        mu.host::clear();
        rho.host::clear();
        //
        lambda_smo.host::clear();
        mu_smo.host::clear();
        rho_smo.host::clear();
        //
        lambda_inv.host::clear();
        mu_inv.host::clear();
        rho_inv.host::clear();

        vp.device::clear();
    }
    // //
    // void PoyntingVec::cuda_cal_poynting_vec()
    // {
    //     cal_poynting_vec<<<  jarvis_cuda_kernel_size(geomodel_p->gridphy.n_elem), 0, *cal_stream_p>>>(
    //         geomodel_p->gridext, geomodel_p->gridphy, geomodel_p->margin,
    //         sub_ext_wavefield_p->vx.device::ptr(), sub_ext_wavefield_p->vy.device::ptr(), sub_ext_wavefield_p->vz.device::ptr(),
    //         sub_ext_wavefield_p->sxx.device::ptr(), sub_ext_wavefield_p->sxy.device::ptr(), sub_ext_wavefield_p->sxz.device::ptr(),
    //         sub_ext_wavefield_p->syy.device::ptr(), sub_ext_wavefield_p->syz.device::ptr(), sub_ext_wavefield_p->szz.device::ptr(), pyt_vec.device::ptr());
    // }
    // //
    // void PoyntingVec::cuda_get_rece_pyt(int it)
    // {
    //     get_rece_pyt<<<1, survey_p->sub_recv.n_elem, 0, *cal_stream_p>>>(
    //         geomodel_p->gridext, geomodel_p->margin,
    //         it, seis_record_p->ntime, survey_p->sub_recv.n_elem, survey_p->sub_recv.device::ptr(),
    //         sub_ext_wavefield_p->vx.device::ptr(), sub_ext_wavefield_p->vy.device::ptr(), sub_ext_wavefield_p->vz.device::ptr(),
    //         sub_ext_wavefield_p->sxx.device::ptr(), sub_ext_wavefield_p->sxy.device::ptr(), sub_ext_wavefield_p->sxz.device::ptr(),
    //         sub_ext_wavefield_p->syy.device::ptr(), sub_ext_wavefield_p->syz.device::ptr(), sub_ext_wavefield_p->szz.device::ptr(), pyt_vec_rece.device::ptr());
    // }
    // //
    // void PoyntingVec::cudax_fliter_rece_pyt(Position pos_flag, Field<float,MemType::paged_device> *_seis)
    // {
    //     int block_num = (seis_record_p->ntime * survey_p->sub_recv.n_elem + jarvis_const::block_size - 1) / jarvis_const::block_size;
    //     int thread_num = jarvis_const::block_size;
    //     fliter_rece_pyt<<<block_num, thread_num, 0, *cal_stream_p>>>(pos_flag, seis_record_p->ntime, _seis->cu_mem, pyt_vec_rece.device::ptr());
    // }
    //
}
#endif