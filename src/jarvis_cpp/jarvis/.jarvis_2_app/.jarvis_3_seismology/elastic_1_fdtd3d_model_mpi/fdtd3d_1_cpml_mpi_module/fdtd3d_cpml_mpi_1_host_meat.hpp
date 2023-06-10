#pragma once
#ifndef _FDTD3D_CPML_MPI_HOST_MEAT_HPP
#define _FDTD3D_CPML_MPI_HOST_MEAT_HPP
#include "fdtd3d_cpml_mpi_0.hpp"
namespace jarvis
{
    inline void mpicuCPML::cu_alloc_damp_coeff() // 3
    {
        if (is_pml.is_top)
        {
            d_x.top.alloc(damp_host_memory, grid_top);
            alpha_x.top.alloc(damp_host_memory, grid_top);
            a_x.top.alloc(damp_device_memory, grid_top);
            b_x.top.alloc(damp_device_memory, grid_top);
            //
            d_y.top.alloc(damp_host_memory, grid_top);
            alpha_y.top.alloc(damp_host_memory, grid_top);
            a_y.top.alloc(damp_device_memory, grid_top);
            b_y.top.alloc(damp_device_memory, grid_top);
            //
            //
            d_z.top.alloc(damp_host_memory, grid_top);
            alpha_z.top.alloc(damp_host_memory, grid_top);
            a_z.top.alloc(damp_device_memory, grid_top);
            b_z.top.alloc(damp_device_memory, grid_top);
            //
            d_x_half.top.alloc(damp_host_memory, grid_top);
            alpha_x_half.top.alloc(damp_host_memory, grid_top);
            a_x_half.top.alloc(damp_device_memory, grid_top);
            b_x_half.top.alloc(damp_device_memory, grid_top);
            //
            d_y_half.top.alloc(damp_host_memory, grid_top);
            alpha_y_half.top.alloc(damp_host_memory, grid_top);
            a_y_half.top.alloc(damp_device_memory, grid_top);
            b_y_half.top.alloc(damp_device_memory, grid_top);
            //
            d_z_half.top.alloc(damp_host_memory, grid_top);
            alpha_z_half.top.alloc(damp_host_memory, grid_top);
            a_z_half.top.alloc(damp_device_memory, grid_top);
            b_z_half.top.alloc(damp_device_memory, grid_top);
        }
        if (is_pml.is_bottom)
        {
            d_x.bottom.alloc(damp_host_memory, grid_bottom);
            alpha_x.bottom.alloc(damp_host_memory, grid_bottom);
            a_x.bottom.alloc(damp_device_memory, grid_bottom);
            b_x.bottom.alloc(damp_device_memory, grid_bottom);
            //
            d_y.bottom.alloc(damp_host_memory, grid_bottom);
            alpha_y.bottom.alloc(damp_host_memory, grid_bottom);
            a_y.bottom.alloc(damp_device_memory, grid_bottom);
            b_y.bottom.alloc(damp_device_memory, grid_bottom);
            //
            //
            d_z.bottom.alloc(damp_host_memory, grid_bottom);
            alpha_z.bottom.alloc(damp_host_memory, grid_bottom);
            a_z.bottom.alloc(damp_device_memory, grid_bottom);
            b_z.bottom.alloc(damp_device_memory, grid_bottom);
            //
            d_x_half.bottom.alloc(damp_host_memory, grid_bottom);
            alpha_x_half.bottom.alloc(damp_host_memory, grid_bottom);
            a_x_half.bottom.alloc(damp_device_memory, grid_bottom);
            b_x_half.bottom.alloc(damp_device_memory, grid_bottom);
            //
            d_y_half.bottom.alloc(damp_host_memory, grid_bottom);
            alpha_y_half.bottom.alloc(damp_host_memory, grid_bottom);
            a_y_half.bottom.alloc(damp_device_memory, grid_bottom);
            b_y_half.bottom.alloc(damp_device_memory, grid_bottom);
            //
            d_z_half.bottom.alloc(damp_host_memory, grid_bottom);
            alpha_z_half.bottom.alloc(damp_host_memory, grid_bottom);
            a_z_half.bottom.alloc(damp_device_memory, grid_bottom);
            b_z_half.bottom.alloc(damp_device_memory, grid_bottom);
        }
        if (is_pml.is_front)
        {
            d_x.front.alloc(damp_host_memory, grid_front);
            alpha_x.front.alloc(damp_host_memory, grid_front);
            a_x.front.alloc(damp_device_memory, grid_front);
            b_x.front.alloc(damp_device_memory, grid_front);
            //
            d_y.front.alloc(damp_host_memory, grid_front);
            alpha_y.front.alloc(damp_host_memory, grid_front);
            a_y.front.alloc(damp_device_memory, grid_front);
            b_y.front.alloc(damp_device_memory, grid_front);
            //
            d_x_half.front.alloc(damp_host_memory, grid_front);
            alpha_x_half.front.alloc(damp_host_memory, grid_front);
            a_x_half.front.alloc(damp_device_memory, grid_front);
            b_x_half.front.alloc(damp_device_memory, grid_front);
            //
            d_y_half.front.alloc(damp_host_memory, grid_front);
            alpha_y_half.front.alloc(damp_host_memory, grid_front);
            a_y_half.front.alloc(damp_device_memory, grid_front);
            b_y_half.front.alloc(damp_device_memory, grid_front);
        }
        if (is_pml.is_back)
        {
            d_x.back.alloc(damp_host_memory, grid_back);
            alpha_x.back.alloc(damp_host_memory, grid_back);
            a_x.back.alloc(damp_device_memory, grid_back);
            b_x.back.alloc(damp_device_memory, grid_back);
            //
            d_y.back.alloc(damp_host_memory, grid_back);
            alpha_y.back.alloc(damp_host_memory, grid_back);
            a_y.back.alloc(damp_device_memory, grid_back);
            b_y.back.alloc(damp_device_memory, grid_back);
            //
            d_x_half.back.alloc(damp_host_memory, grid_back);
            alpha_x_half.back.alloc(damp_host_memory, grid_back);
            a_x_half.back.alloc(damp_device_memory, grid_back);
            b_x_half.back.alloc(damp_device_memory, grid_back);
            //
            d_y_half.back.alloc(damp_host_memory, grid_back);
            alpha_y_half.back.alloc(damp_host_memory, grid_back);
            a_y_half.back.alloc(damp_device_memory, grid_back);
            b_y_half.back.alloc(damp_device_memory, grid_back);
        }
        if (is_pml.is_left)
        {
            d_x.left.alloc(damp_host_memory, grid_left);
            alpha_x.left.alloc(damp_host_memory, grid_left);
            a_x.left.alloc(damp_device_memory, grid_left);
            b_x.left.alloc(damp_device_memory, grid_left);
            //
            d_x_half.left.alloc(damp_host_memory, grid_left);
            alpha_x_half.left.alloc(damp_host_memory, grid_left);
            a_x_half.left.alloc(damp_device_memory, grid_left);
            b_x_half.left.alloc(damp_device_memory, grid_left);
        }
        if (is_pml.is_right)
        {
            d_x.right.alloc(damp_host_memory, grid_right);
            alpha_x.right.alloc(damp_host_memory, grid_right);
            a_x.right.alloc(damp_device_memory, grid_right);
            b_x.right.alloc(damp_device_memory, grid_right);
            //
            d_x_half.right.alloc(damp_host_memory, grid_right);
            alpha_x_half.right.alloc(damp_host_memory, grid_right);
            a_x_half.right.alloc(damp_device_memory, grid_right);
            b_x_half.right.alloc(damp_device_memory, grid_right);
        }
        damp_host_memory.alloc();
        damp_device_memory.alloc();
    }
    //
    inline void mpicuCPML::calculate_damp_coeff() // 4th
    {
        elasticGeoModel<domain::local> &model = *geomodel_p; // link the elasticGeoModel
        float DAMP_COEF = 7.0f;
        int NUMB_POWER = 2;
        float ALPHA_MAX_PML = jarvis_const::pi * seis_record_p->fm;
        // Not include the count of the half of the spatial order.
        float top_pml_num = float(model.margin.top_margin - geo_const::pml_fdorder_half);
        float bottom_pml_num = float(model.margin.bottom_margin - geo_const::pml_fdorder_half);
        float front_pml_num = float(model.margin.front_margin - geo_const::pml_fdorder_half);
        float back_pml_num = float(model.margin.back_margin - geo_const::pml_fdorder_half);
        float right_pml_num = float(model.margin.right_margin - geo_const::pml_fdorder_half);
        float left_pml_num = float(model.margin.left_margin - geo_const::pml_fdorder_half);
        // top
        if (is_pml.is_top)
        {
            int pad_front = 0;
            if (model.margin.front_margin >= geo_const::pml_fdorder_half)
                pad_front = geo_const::pml_fdorder_half;

            int pad_back = 0;
            if (model.margin.back_margin >= geo_const::pml_fdorder_half)
                pad_back = geo_const::pml_fdorder_half;

            int pad_left = 0;
            if (model.margin.left_margin >= geo_const::pml_fdorder_half)
                pad_left = geo_const::pml_fdorder_half;

            int pad_right = 0;
            if (model.margin.right_margin >= geo_const::pml_fdorder_half)
                pad_right = geo_const::pml_fdorder_half;

            for (int k = geo_const::pml_fdorder_half; k < model.margin.top_margin; k++)
            {
                for (int i = pad_right; i < model.gridext.n_rows - pad_left; i++)
                {
                    for (int j = pad_front; j < model.gridext.n_cols - pad_back; j++)
                    {
                        float temp_z = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * top_pml_num * model.gridext.d_slices);
                        d_z.top(i, j, k - geo_const::pml_fdorder_half) = temp_z * pow((model.margin.top_margin - k) / top_pml_num, NUMB_POWER);
                        alpha_z.top(i, j, k - geo_const::pml_fdorder_half) = ALPHA_MAX_PML * (1.0f - (model.margin.top_margin - k) / top_pml_num);
                        d_z_half.top(i, j, k - geo_const::pml_fdorder_half) = temp_z * pow((model.margin.top_margin - k - 0.5f) / top_pml_num, NUMB_POWER);
                        alpha_z_half.top(i, j, k - geo_const::pml_fdorder_half) = ALPHA_MAX_PML * (1.0f - (model.margin.top_margin - k - 0.5f) / top_pml_num);
                    }
                    for (int j = geo_const::pml_fdorder_half; j < model.margin.front_margin; j++)
                    {
                        float temp_y = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * front_pml_num * model.gridext.d_cols);
                        d_y.top(i, j, k - geo_const::pml_fdorder_half) = temp_y * pow((model.margin.front_margin - j) / front_pml_num, NUMB_POWER);
                        alpha_y.top(i, j, k - geo_const::pml_fdorder_half) = ALPHA_MAX_PML * (1.0f - (model.margin.front_margin - j) / front_pml_num);
                        d_y_half.top(i, j, k - geo_const::pml_fdorder_half) = temp_y * pow((model.margin.front_margin - j - 0.5f) / front_pml_num, NUMB_POWER);
                        alpha_y_half.top(i, j, k - geo_const::pml_fdorder_half) = ALPHA_MAX_PML * (1.0f - (model.margin.front_margin - j - 0.5f) / front_pml_num);
                    }
                    for (int j = (model.gridext.n_cols - model.margin.back_margin); j < (model.gridext.n_cols - geo_const::pml_fdorder_half); j++)
                    {
                        float temp_y = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * back_pml_num * model.gridext.d_cols);
                        d_y.top(i, j, k - geo_const::pml_fdorder_half) = temp_y * pow((j - (model.gridext.n_cols - model.margin.back_margin - 1) - 0.5f) / back_pml_num, NUMB_POWER);
                        alpha_y.top(i, j, k - geo_const::pml_fdorder_half) = ALPHA_MAX_PML * (1.0f - (j - (model.gridext.n_cols - model.margin.back_margin - 1) - 0.5f) / back_pml_num);
                        d_y_half.top(i, j, k - geo_const::pml_fdorder_half) = temp_y * pow(((j - (model.gridext.n_cols - model.margin.back_margin - 1))) / back_pml_num, NUMB_POWER);
                        alpha_y_half.top(i, j, k - geo_const::pml_fdorder_half) = ALPHA_MAX_PML * (1.0f - ((j - (model.gridext.n_cols - model.margin.back_margin - 1))) / back_pml_num);
                    }
                }
                for (int j = pad_front; j < model.gridext.n_cols - pad_back; j++)
                {
                    for (int i = geo_const::pml_fdorder_half; i < model.margin.right_margin; i++)
                    {
                        float temp_x = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * right_pml_num * model.gridext.d_rows);
                        d_x.top(i, j, k - geo_const::pml_fdorder_half) = temp_x * pow((model.margin.right_margin - i) / right_pml_num, NUMB_POWER);
                        alpha_x.top(i, j, k - geo_const::pml_fdorder_half) = ALPHA_MAX_PML * (1.0f - (model.margin.right_margin - i) / right_pml_num);
                        d_x_half.top(i, j, k - geo_const::pml_fdorder_half) = temp_x * pow((model.margin.right_margin - i - 0.5f) / right_pml_num, NUMB_POWER);
                        alpha_x_half.top(i, j, k - geo_const::pml_fdorder_half) = ALPHA_MAX_PML * (1.0f - (model.margin.right_margin - i - 0.5f) / right_pml_num);
                    }
                    for (int i = (model.gridext.n_rows - model.margin.left_margin); i < (model.gridext.n_rows - geo_const::pml_fdorder_half); i++)
                    {
                        float temp_x = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * left_pml_num * model.gridext.d_rows);
                        d_x.top(i, j, k - geo_const::pml_fdorder_half) = temp_x * pow((i - (model.gridext.n_rows - model.margin.left_margin - 1) - 0.5f) / left_pml_num, NUMB_POWER);
                        alpha_x.top(i, j, k - geo_const::pml_fdorder_half) = ALPHA_MAX_PML * (1.0f - (i - (model.gridext.n_rows - model.margin.left_margin - 1) - 0.5f) / left_pml_num);
                        d_x_half.top(i, j, k - geo_const::pml_fdorder_half) = temp_x * pow((i - (model.gridext.n_rows - model.margin.left_margin - 1)) / left_pml_num, NUMB_POWER);
                        alpha_x_half.top(i, j, k - geo_const::pml_fdorder_half) = ALPHA_MAX_PML * (1.0f - (i - (model.gridext.n_rows - model.margin.left_margin - 1)) / left_pml_num);
                    }
                }
            }
            for (int i = 0; i < grid_top.n_elem; i++)
            {
                if (alpha_x.top[i] < 0)
                    alpha_x.top[i] = 0;
                if (alpha_x_half.top[i] < 0)
                    alpha_x_half.top[i] = 0;
                if (d_x.top[i] >= 0)
                {
                    b_x.top[i] = exp(-(d_x.top[i] + alpha_x.top[i]) * seis_record_p->dt);
                    b_x_half.top[i] = exp(-(d_x_half.top[i] + alpha_x_half.top[i]) * seis_record_p->dt);
                }
                if (abs(d_x.top[i]) >= 1e-8)
                {
                    a_x.top[i] = d_x.top[i] * (b_x.top[i] - 1.0f) / ((d_x.top[i] + alpha_x.top[i]));
                    a_x_half.top[i] = d_x_half.top[i] * (b_x_half.top[i] - 1.0f) / ((d_x_half.top[i] + alpha_x_half.top[i]));
                }
                if (alpha_y.top[i] < 0)
                    alpha_y.top[i] = 0;
                if (alpha_y_half.top[i] < 0)
                    alpha_y_half.top[i] = 0;

                if (d_y.top[i] >= 0)
                {
                    b_y.top[i] = exp(-(d_y.top[i] + alpha_y.top[i]) * seis_record_p->dt);
                    b_y_half.top[i] = exp(-(d_y_half.top[i] + alpha_y_half.top[i]) * seis_record_p->dt);
                }
                if (abs(d_y.top[i]) >= 1e-8)
                {
                    a_y.top[i] = d_y.top[i] * (b_y.top[i] - 1.0f) / ((d_y.top[i] + alpha_y.top[i]));
                    a_y_half.top[i] = d_y_half.top[i] * (b_y_half.top[i] - 1.0f) / ((d_y_half.top[i] + alpha_y_half.top[i]));
                }

                if (alpha_z.top[i] < 0)
                    alpha_z.top[i] = 0;
                if (alpha_z_half.top[i] < 0)
                    alpha_z_half.top[i] = 0;

                if (d_z.top[i] >= 0)
                {
                    b_z.top[i] = exp(-(d_z.top[i] + alpha_z.top[i]) * seis_record_p->dt);
                    b_z_half.top[i] = exp(-(d_z_half.top[i] + alpha_z_half.top[i]) * seis_record_p->dt);
                }
                if (abs(d_z.top[i]) >= 1e-8)
                {
                    a_z.top[i] = d_z.top[i] * (b_z.top[i] - 1.0f) / ((d_z.top[i] + alpha_z.top[i]));
                    a_z_half.top[i] = d_z_half.top[i] * (b_z_half.top[i] - 1.0f) / ((d_z_half.top[i] + alpha_z_half.top[i]));
                }
            }
        }
        // bottom
        if (is_pml.is_bottom)
        {
            int pad_right = 0;
            if (model.margin.right_margin >= geo_const::pml_fdorder_half)
                pad_right = geo_const::pml_fdorder_half;

            int pad_left = 0;
            if (model.margin.left_margin >= geo_const::pml_fdorder_half)
                pad_left = geo_const::pml_fdorder_half;

            int pad_front = 0;
            if (model.margin.front_margin >= geo_const::pml_fdorder_half)
                pad_front = geo_const::pml_fdorder_half;

            int pad_back = 0;
            if (model.margin.back_margin >= geo_const::pml_fdorder_half)
                pad_back = geo_const::pml_fdorder_half;

            for (int k = (model.gridext.n_slices - model.margin.bottom_margin); k < (model.gridext.n_slices - geo_const::pml_fdorder_half); k++)
            {
                for (int i = pad_right; i < model.gridext.n_rows - pad_left; i++)
                {
                    for (int j = pad_front; j < model.gridext.n_cols - pad_back; j++)
                    // bottom dampx
                    {
                        float temp_z = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * bottom_pml_num * model.gridext.d_slices);
                        d_z.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = temp_z * pow((k - (model.gridext.n_slices - model.margin.bottom_margin - 1) - 0.5f) / bottom_pml_num, NUMB_POWER);
                        alpha_z.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = ALPHA_MAX_PML * (1.0f - (k - (model.gridext.n_slices - model.margin.bottom_margin - 1) - 0.5f) / bottom_pml_num);
                        d_z_half.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = temp_z * pow((k - (model.gridext.n_slices - model.margin.bottom_margin - 1)) / bottom_pml_num, NUMB_POWER);
                        alpha_z_half.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = ALPHA_MAX_PML * (1.0f - (k - (model.gridext.n_slices - model.margin.bottom_margin - 1)) / bottom_pml_num);
                    }
                    // bottom dampy
                    for (int j = geo_const::pml_fdorder_half; j < model.margin.front_margin; j++)
                    {
                        float temp_y = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * front_pml_num * model.gridext.d_cols);

                        d_y.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = temp_y * pow((model.margin.front_margin - j) / front_pml_num, NUMB_POWER);
                        alpha_y.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = ALPHA_MAX_PML * (1.0f - (model.margin.front_margin - j) / front_pml_num);
                        d_y_half.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = temp_y * pow((model.margin.front_margin - j - 0.5f) / front_pml_num, NUMB_POWER);
                        alpha_y_half.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = ALPHA_MAX_PML * (1.0f - (model.margin.front_margin - j - 0.5f) / front_pml_num);
                    }
                    for (int j = (model.gridext.n_cols - model.margin.back_margin); j < (model.gridext.n_cols - geo_const::pml_fdorder_half); j++)
                    {
                        float temp_y = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * back_pml_num * model.gridext.d_cols);

                        d_y.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = temp_y * pow((j - (model.gridext.n_cols - model.margin.back_margin - 1) - 0.5f) / back_pml_num, NUMB_POWER);
                        alpha_y.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = ALPHA_MAX_PML * (1.0f - (j - (model.gridext.n_cols - model.margin.back_margin - 1) - 0.5f) / back_pml_num);
                        d_y_half.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = temp_y * pow((j - (model.gridext.n_cols - model.margin.back_margin - 1)) / back_pml_num, NUMB_POWER);
                        alpha_y_half.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = ALPHA_MAX_PML * (1.0f - (j - (model.gridext.n_cols - model.margin.back_margin - 1)) / back_pml_num);
                    }
                }

                for (int j = pad_front; j < model.gridext.n_cols - pad_back; j++)
                {
                    for (int i = geo_const::pml_fdorder_half; i < model.margin.right_margin; i++)
                    {
                        float temp_x = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * right_pml_num * model.gridext.d_rows);

                        d_x.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = temp_x * pow((model.margin.right_margin - i) / right_pml_num, NUMB_POWER);
                        alpha_x.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = ALPHA_MAX_PML * (1.0f - (model.margin.right_margin - i) / right_pml_num);
                        d_x_half.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = temp_x * pow((model.margin.right_margin - i - 0.5f) / right_pml_num, NUMB_POWER);
                        alpha_x_half.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = ALPHA_MAX_PML * (1.0f - (model.margin.right_margin - i - 0.5f) / right_pml_num);
                    }
                    for (int i = (model.gridext.n_rows - model.margin.left_margin); i < (model.gridext.n_rows - geo_const::pml_fdorder_half); i++)
                    {
                        float temp_x = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * left_pml_num * model.gridext.d_rows);

                        d_x.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = temp_x * pow((i - (model.gridext.n_rows - model.margin.left_margin - 1) - 0.5f) / left_pml_num, NUMB_POWER);
                        alpha_x.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = ALPHA_MAX_PML * (1.0f - (i - (model.gridext.n_rows - model.margin.left_margin - 1) - 0.5f) / left_pml_num);
                        d_x_half.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = temp_x * pow((i - (model.gridext.n_rows - model.margin.left_margin - 1)) / left_pml_num, NUMB_POWER);
                        alpha_x_half.bottom(i, j, k - (model.gridext.n_slices - model.margin.bottom_margin)) = ALPHA_MAX_PML * (1.0f - (i - (model.gridext.n_rows - model.margin.left_margin - 1)) / left_pml_num);
                    }
                }
            }

            for (int i = 0; i < grid_bottom.n_elem; i++)
            {
                if (alpha_x.bottom[i] < 0)
                    alpha_x.bottom[i] = 0;
                if (alpha_x_half.bottom[i] < 0)
                    alpha_x_half.bottom[i] = 0;
                if (d_x.bottom[i] >= 0)
                {
                    b_x.bottom[i] = exp(-(d_x.bottom[i] + alpha_x.bottom[i]) * seis_record_p->dt);
                    b_x_half.bottom[i] = exp(-(d_x_half.bottom[i] + alpha_x_half.bottom[i]) * seis_record_p->dt);
                }
                if (abs(d_x.bottom[i]) >= 1e-8)
                {
                    a_x.bottom[i] = d_x.bottom[i] * (b_x.bottom[i] - 1.0f) / ((d_x.bottom[i] + alpha_x.bottom[i]));
                    a_x_half.bottom[i] = d_x_half.bottom[i] * (b_x_half.bottom[i] - 1.0f) / ((d_x_half.bottom[i] + alpha_x_half.bottom[i]));
                }
                if (alpha_y.bottom[i] < 0)
                    alpha_y.bottom[i] = 0;
                if (alpha_y_half.bottom[i] < 0)
                    alpha_y_half.bottom[i] = 0;
                if (d_y.bottom[i] >= 0)
                {
                    b_y.bottom[i] = exp(-(d_y.bottom[i] + alpha_y.bottom[i]) * seis_record_p->dt);
                    b_y_half.bottom[i] = exp(-(d_y_half.bottom[i] + alpha_y_half.bottom[i]) * seis_record_p->dt);
                }
                if (abs(d_y.bottom[i]) >= 1e-8)
                {
                    a_y.bottom[i] = d_y.bottom[i] * (b_y.bottom[i] - 1.0f) / ((d_y.bottom[i] + alpha_y.bottom[i]));
                    a_y_half.bottom[i] = d_y_half.bottom[i] * (b_y_half.bottom[i] - 1.0f) / ((d_y_half.bottom[i] + alpha_y_half.bottom[i]));
                }
                if (alpha_z.bottom[i] < 0)
                    alpha_z.bottom[i] = 0;
                if (alpha_z_half.bottom[i] < 0)
                    alpha_z_half.bottom[i] = 0;
                if (d_z.bottom[i] >= 0)
                {
                    b_z.bottom[i] = exp(-(d_z.bottom[i] + alpha_z.bottom[i]) * seis_record_p->dt);
                    b_z_half.bottom[i] = exp(-(d_z_half.bottom[i] + alpha_z_half.bottom[i]) * seis_record_p->dt);
                }
                if (abs(d_z.bottom[i]) > 1e-8)
                {
                    a_z.bottom[i] = d_z.bottom[i] * (b_z.bottom[i] - 1.0f) / ((d_z.bottom[i] + alpha_z.bottom[i]));
                    a_z_half.bottom[i] = d_z_half.bottom[i] * (b_z_half.bottom[i] - 1.0f) / ((d_z_half.bottom[i] + alpha_z_half.bottom[i]));
                }
            }
        }
        // back
        if (is_pml.is_back)
        {
            int pad_right = 0;
            if (model.margin.right_margin >= geo_const::pml_fdorder_half)
                pad_right = geo_const::pml_fdorder_half;

            int pad_left = 0;
            if (model.margin.left_margin >= geo_const::pml_fdorder_half)
                pad_left = geo_const::pml_fdorder_half;

            for (int j = (model.gridext.n_cols - model.margin.back_margin); j < (model.gridext.n_cols - geo_const::pml_fdorder_half); j++)
            {
                for (int k = model.margin.top_margin; k < model.gridext.n_slices - model.margin.bottom_margin; k++)
                {
                    for (int i = pad_right; i < model.gridext.n_rows - pad_left; i++)
                    {
                        float temp_y = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * back_pml_num * model.gridext.d_cols);
                        d_y.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = temp_y * pow((j - (model.gridext.n_cols - model.margin.back_margin - 1) - 0.5f) / back_pml_num, NUMB_POWER);
                        alpha_y.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (j - (model.gridext.n_cols - model.margin.back_margin - 1) - 0.5f) / back_pml_num);
                        d_y_half.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = temp_y * pow((j - (model.gridext.n_cols - model.margin.back_margin - 1)) / back_pml_num, NUMB_POWER);
                        alpha_y_half.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (j - (model.gridext.n_cols - model.margin.back_margin - 1)) / back_pml_num);
                    }
                    for (int i = geo_const::pml_fdorder_half; i < model.margin.right_margin; i++)
                    {
                        float temp_x = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * right_pml_num * model.gridext.d_rows);
                        d_x.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = temp_x * pow((model.margin.right_margin - i) / right_pml_num, NUMB_POWER);
                        alpha_x.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (model.margin.right_margin - i) / right_pml_num);
                        d_x_half.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = temp_x * pow((model.margin.right_margin - i - 0.5f) / right_pml_num, NUMB_POWER);
                        alpha_x_half.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (model.margin.right_margin - i - 0.5f) / right_pml_num);
                    }
                    for (int i = (model.gridext.n_rows - model.margin.left_margin); i < (model.gridext.n_rows - geo_const::pml_fdorder_half); i++)
                    {
                        float temp_x = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * left_pml_num * model.gridext.d_rows);
                        d_x.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = temp_x * pow((i - (model.gridext.n_rows - model.margin.left_margin - 1) - 0.5f) / left_pml_num, NUMB_POWER);
                        alpha_x.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (i - (model.gridext.n_rows - model.margin.left_margin - 1) - 0.5f) / left_pml_num);
                        d_x_half.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = temp_x * pow((i - (model.gridext.n_rows - model.margin.left_margin - 1)) / left_pml_num, NUMB_POWER);
                        alpha_x_half.back(i, j - (model.gridext.n_cols - model.margin.back_margin), k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (i - (model.gridext.n_rows - model.margin.left_margin - 1)) / left_pml_num);
                    }
                }
            }

            for (int i = 0; i < grid_back.n_elem; i++)
            {
                if (alpha_x.back[i] < 0)
                    alpha_x.back[i] = 0;
                if (alpha_x_half.back[i] < 0)
                    alpha_x_half.back[i] = 0;
                if (d_x.back[i] >= 0)
                {
                    b_x.back[i] = exp(-(d_x.back[i] + alpha_x.back[i]) * seis_record_p->dt);
                    b_x_half.back[i] = exp(-(d_x_half.back[i] + alpha_x_half.back[i]) * seis_record_p->dt);
                }
                if (abs(d_x.back[i]) >= 1e-8)
                {
                    a_x.back[i] = d_x.back[i] * (b_x.back[i] - 1.0f) / ((d_x.back[i] + alpha_x.back[i]));
                    a_x_half.back[i] = d_x_half.back[i] * (b_x_half.back[i] - 1.0f) / ((d_x_half.back[i] + alpha_x_half.back[i]));
                }
                if (alpha_y.back[i] < 0)
                    alpha_y.back[i] = 0;
                if (alpha_y_half.back[i] < 0)
                    alpha_y_half.back[i] = 0;
                if (d_y.back[i] >= 0)
                {
                    b_y.back[i] = exp(-(d_y.back[i] + alpha_y.back[i]) * seis_record_p->dt);
                    b_y_half.back[i] = exp(-(d_y_half.back[i] + alpha_y_half.back[i]) * seis_record_p->dt);
                }
                if (abs(d_y.back[i]) >= 1e-8)
                {
                    a_y.back[i] = d_y.back[i] * (b_y.back[i] - 1.0f) / ((d_y.back[i] + alpha_y.back[i]));
                    a_y_half.back[i] = d_y_half.back[i] * (b_y_half.back[i] - 1.0f) / ((d_y_half.back[i] + alpha_y_half.back[i]));
                }
            }
        }
        // front
        if (is_pml.is_front)
        {
            int pad_right = 0;
            if (model.margin.right_margin >= geo_const::pml_fdorder_half)
                pad_right = geo_const::pml_fdorder_half;

            int pad_left = 0;
            if (model.margin.left_margin >= geo_const::pml_fdorder_half)
                pad_left = geo_const::pml_fdorder_half;

            for (int j = geo_const::pml_fdorder_half; j < model.margin.front_margin; j++)
            {
                for (int k = model.margin.top_margin; k < model.gridext.n_slices - model.margin.bottom_margin; k++)
                {
                    for (int i = pad_right; i < model.gridext.n_rows - pad_left; i++)
                    {
                        float temp_y = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * front_pml_num * model.gridext.d_cols);

                        d_y.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = temp_y * pow((model.margin.front_margin - j) / front_pml_num, NUMB_POWER);
                        alpha_y.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (model.margin.front_margin - j) / front_pml_num);
                        d_y_half.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = temp_y * pow((model.margin.front_margin - j - 0.5f) / front_pml_num, NUMB_POWER);
                        alpha_y_half.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (model.margin.front_margin - j - 0.5f) / front_pml_num);
                    }
                    for (int i = geo_const::pml_fdorder_half; i < model.margin.right_margin; i++)
                    {
                        float temp_x = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * right_pml_num * model.gridext.d_rows);

                        d_x.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = temp_x * pow((model.margin.right_margin - i) / right_pml_num, NUMB_POWER);
                        alpha_x.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (model.margin.right_margin - i) / right_pml_num);
                        d_x_half.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = temp_x * pow((model.margin.right_margin - i - 0.5f) / right_pml_num, NUMB_POWER);
                        alpha_x_half.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (model.margin.right_margin - i - 0.5f) / right_pml_num);
                    }
                    for (int i = (model.gridext.n_rows - model.margin.left_margin); i < (model.gridext.n_rows - geo_const::pml_fdorder_half); i++)
                    {
                        float temp_x = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * left_pml_num * model.gridext.d_rows);

                        d_x.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = temp_x * pow((i - (model.gridext.n_rows - model.margin.left_margin - 1) - 0.5f) / left_pml_num, NUMB_POWER);
                        alpha_x.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (i - (model.gridext.n_rows - model.margin.left_margin - 1) - 0.5f) / left_pml_num);
                        d_x_half.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = temp_x * pow((i - (model.gridext.n_rows - model.margin.left_margin - 1)) / left_pml_num, NUMB_POWER);
                        alpha_x_half.front(i, j - geo_const::pml_fdorder_half, k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (i - (model.gridext.n_rows - model.margin.left_margin - 1)) / left_pml_num);
                    }
                }
            }

            for (int i = 0; i < grid_front.n_elem; i++)
            {
                if (alpha_x.front[i] < 0)
                    alpha_x.front[i] = 0;
                if (alpha_x_half.front[i] < 0)
                    alpha_x_half.front[i] = 0;
                if (d_x.front[i] >= 0)
                {
                    b_x.front[i] = exp(-(d_x.front[i] + alpha_x.front[i]) * seis_record_p->dt);
                    b_x_half.front[i] = exp(-(d_x_half.front[i] + alpha_x_half.front[i]) * seis_record_p->dt);
                }
                if (abs(d_x.front[i]) >= 1e-8)
                {
                    a_x.front[i] = d_x.front[i] * (b_x.front[i] - 1.0f) / ((d_x.front[i] + alpha_x.front[i]));
                    a_x_half.front[i] = d_x_half.front[i] * (b_x_half.front[i] - 1.0f) / ((d_x_half.front[i] + alpha_x_half.front[i]));
                }
                if (alpha_y.front[i] < 0)
                    alpha_y.front[i] = 0;
                if (alpha_y_half.front[i] < 0)
                    alpha_y_half.front[i] = 0;
                if (d_y.front[i] >= 0)
                {
                    b_y.front[i] = exp(-(d_y.front[i] + alpha_y.front[i]) * seis_record_p->dt);
                    b_y_half.front[i] = exp(-(d_y_half.front[i] + alpha_y_half.front[i]) * seis_record_p->dt);
                }
                if (abs(d_y.front[i]) >= 1e-8)
                {
                    a_y.front[i] = d_y.front[i] * (b_y.front[i] - 1.0f) / ((d_y.front[i] + alpha_y.front[i]));
                    a_y_half.front[i] = d_y_half.front[i] * (b_y_half.front[i] - 1.0f) / ((d_y_half.front[i] + alpha_y_half.front[i]));
                }
            }
        }

        // right
        if (is_pml.is_right)
            for (int i = geo_const::pml_fdorder_half; i < model.margin.right_margin; i++)
            {
                for (int j = model.margin.front_margin; j < model.gridext.n_cols - model.margin.back_margin; j++)
                {
                    for (int k = model.margin.top_margin; k < model.gridext.n_slices - model.margin.bottom_margin; k++)
                    {
                        float temp_x = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * right_pml_num * model.gridext.d_rows);

                        d_x.right(i - geo_const::pml_fdorder_half, j - model.margin.front_margin, k - model.margin.top_margin) = temp_x * pow((model.margin.right_margin - i) / right_pml_num, NUMB_POWER);
                        alpha_x.right(i - geo_const::pml_fdorder_half, j - model.margin.front_margin, k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (model.margin.right_margin - i) / right_pml_num);
                        d_x_half.right(i - geo_const::pml_fdorder_half, j - model.margin.front_margin, k - model.margin.top_margin) = temp_x * pow((model.margin.right_margin - i - 0.5f) / right_pml_num, NUMB_POWER);
                        alpha_x_half.right(i - geo_const::pml_fdorder_half, j - model.margin.front_margin, k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (model.margin.right_margin - i - 0.5f) / right_pml_num);
                    }
                }
                for (int i = 0; i < grid_right.n_elem; i++)
                {
                    if (alpha_x.right[i] < 0)
                        alpha_x.right[i] = 0;
                    if (alpha_x_half.right[i] < 0)
                        alpha_x_half.right[i] = 0;
                    if (d_x.right[i] >= 0)
                    {
                        b_x.right[i] = exp(-(d_x.right[i] + alpha_x.right[i]) * seis_record_p->dt);
                        b_x_half.right[i] = exp(-(d_x_half.right[i] + alpha_x_half.right[i]) * seis_record_p->dt);
                    }
                    if (abs(d_x.right[i]) >= 1e-8)
                    {
                        a_x.right[i] = d_x.right[i] * (b_x.right[i] - 1.0f) / ((d_x.right[i] + alpha_x.right[i]));
                        a_x_half.right[i] = d_x_half.right[i] * (b_x_half.right[i] - 1.0f) / ((d_x_half.right[i] + alpha_x_half.right[i]));
                    }
                }
            }

        // left
        if (is_pml.is_left)
            for (int i = (model.gridext.n_rows - model.margin.left_margin); i < (model.gridext.n_rows - geo_const::pml_fdorder_half); i++)
            {
                for (int j = model.margin.front_margin; j < model.gridext.n_cols - model.margin.back_margin; j++)
                {
                    for (int k = model.margin.top_margin; k < model.gridext.n_slices - model.margin.bottom_margin; k++)
                    {
                        float temp_x = log10(1.0f / R) * DAMP_COEF * model.vp(i, j, k) / (2.0f * left_pml_num * model.gridext.d_rows);

                        d_x.left(i - (model.gridext.n_rows - model.margin.left_margin), j - model.margin.front_margin, k - model.margin.top_margin) = temp_x * pow((i - (model.gridext.n_rows - model.margin.left_margin - 1) - 0.5f) / left_pml_num, NUMB_POWER);
                        alpha_x.left(i - (model.gridext.n_rows - model.margin.left_margin), j - model.margin.front_margin, k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (i - (model.gridext.n_rows - model.margin.left_margin - 1) - 0.5f) / left_pml_num);
                        d_x_half.left(i - (model.gridext.n_rows - model.margin.left_margin), j - model.margin.front_margin, k - model.margin.top_margin) = temp_x * pow((i - (model.gridext.n_rows - model.margin.left_margin - 1)) / left_pml_num, NUMB_POWER);
                        alpha_x_half.left(i - (model.gridext.n_rows - model.margin.left_margin), j - model.margin.front_margin, k - model.margin.top_margin) = ALPHA_MAX_PML * (1.0f - (i - (model.gridext.n_rows - model.margin.left_margin - 1)) / left_pml_num);
                    }
                }
                for (int i = 0; i < grid_left.n_elem; i++)
                {
                    if (alpha_x.left[i] < 0)
                        alpha_x.left[i] = 0;
                    if (alpha_x_half.left[i] < 0)
                        alpha_x_half.left[i] = 0;
                    if (d_x.left[i] >= 0)
                    {
                        b_x.left[i] = exp(-(d_x.left[i] + alpha_x.left[i]) * seis_record_p->dt);
                        b_x_half.left[i] = exp(-(d_x_half.left[i] + alpha_x_half.left[i]) * seis_record_p->dt);
                    }
                    if (abs(d_x.left[i]) >= 1e-8)
                    {
                        a_x.left[i] = d_x.left[i] * (b_x.left[i] - 1.0f) / ((d_x.left[i] + alpha_x.left[i]));
                        a_x_half.left[i] = d_x_half.left[i] * (b_x_half.left[i] - 1.0f) / ((d_x_half.left[i] + alpha_x_half.left[i]));
                    }
                }
            }
    }

    inline void mpicuCPML::cu_copy_damp_coeff_h2d() // 5th
    {
        damp_device_memory.copy_h2d();
        cudaDeviceSynchronize();
        damp_host_memory.clear();
        damp_device_memory.host::clear();
    }
    //
    inline void mpicuCPML::clear()
    {
        if (this->is_pml.is_top)
            this->top.clear();
        if (this->is_pml.is_bottom)
            this->bottom.clear();
        if (this->is_pml.is_front)
            this->front.clear();
        if (this->is_pml.is_back)
            this->back.clear();
        if (this->is_pml.is_right)
            this->right.clear();
        if (this->is_pml.is_left)
            this->left.clear();
        damp_device_memory.device::clear();
        mpi_multi_halo_vel.clear();
        mpi_multi_halo_stress.clear();
    }
}
#endif