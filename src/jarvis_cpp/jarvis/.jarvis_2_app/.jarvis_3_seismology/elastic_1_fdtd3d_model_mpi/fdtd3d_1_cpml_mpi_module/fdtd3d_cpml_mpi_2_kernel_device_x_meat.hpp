#pragma once
#ifndef _FDTD3D_CPML_MPI_KERNEL_DEVICE_X_MEAT_HPP
#define _FDTD3D_CPML_MPI_KERNEL_DEVICE_X_MEAT_HPP
#include "fdtd3d_cpml_mpi_1_host_meat.hpp"
namespace jarvis
{
	namespace cpml
	{
		inline __device__ void mpi_update_vel_pml_x_dir(Frame gridext, Frame griddir,
														int ext_idx, int reg_idx,
														int ext_i, int ext_j, int ext_k, int i, int j, int k,
														bool is_pml_top, bool is_pml_bottom,
														bool is_pml_front, bool is_pml_back,
														float dt, float *pml_dc,
														float *rho,
														float *vx, float *vy, float *vz,
														float *sxx, float *syy, float *szz,
														float *sxy, float *sxz, float *syz,
														//
														float *dir_vx_x, float *dir_vy_x, float *dir_vz_x,

														float *dir_a_x,
														float *dir_a_x_half,
														float *dir_b_x,
														float *dir_b_x_half,
														//
														float *top_sxz, float *top_syz, float *top_szz,
														float *bottom_sxz, float *bottom_syz, float *bottom_szz,
														float *front_sxy, float *front_syy, float *front_syz,
														float *back_sxy, float *back_syy, float *back_syz)
		{
			set_frame_nd(gridext);
			set_frame_nd_pre(dir, griddir);
			//
			float sxx_x_tmp = 0.f, sxy_y_tmp = 0.f, sxz_z_tmp = 0.f,
				  syx_x_tmp = 0.f, syy_y_tmp = 0.f, syz_z_tmp = 0.f,
				  szx_x_tmp = 0.f, szy_y_tmp = 0.f, szz_z_tmp = 0.f;
			//
			float dt_rho = dt / rho[ext_idx];
			//
			int ir;
			//*Y
			if (j >= geo_const::pml_fdorder_half && j <= dir_n_cols - geo_const::pml_fdorder_half - 1)
			{
#pragma unroll
				for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
				{
					syy_y_tmp += pml_dc[ir] * (syy[ext_idx + (ir + 1) * n_rows] - syy[ext_idx - ir * n_rows]);
					sxy_y_tmp += pml_dc[ir] * (sxy[ext_idx + ir * n_rows] - sxy[ext_idx - (ir + 1) * n_rows]);
					szy_y_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_rows] - syz[ext_idx - (ir + 1) * n_rows]);
				}
			}
			else if (j < geo_const::pml_fdorder_half)
			{
				if (is_pml_front)
				{
#pragma unroll
					for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
					{
						syy_y_tmp += pml_dc[ir] * (syy[ext_idx + (ir + 1) * n_rows] - syy[ext_idx - ir * n_rows]);
						sxy_y_tmp += pml_dc[ir] * (sxy[ext_idx + ir * n_rows] - sxy[ext_idx - (ir + 1) * n_rows]);
						szy_y_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_rows] - syz[ext_idx - (ir + 1) * n_rows]);
					}
				}
				else
				{
					for (ir = 0; ir <= j; ++ir)
						syy_y_tmp += pml_dc[ir] * (syy[ext_idx + (ir + 1) * n_rows] - syy[ext_idx - ir * n_rows]);

					for (ir = j + 1; ir < geo_const::pml_fdorder_half; ++ir)
					{
						syy_y_tmp += pml_dc[ir] * (syy[ext_idx + (ir + 1) * n_rows] - front_syy[i + (j - ir - 1 + geo_const::pml_fdorder_half) * dir_n_rows + k * dir_n_rows * (geo_const::pml_fdorder_half - 1)]);
					}

					for (ir = 0; ir < j; ++ir)
					{
						sxy_y_tmp += pml_dc[ir] * (sxy[ext_idx + ir * n_rows] - sxy[ext_idx - (ir + 1) * n_rows]);
						szy_y_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_rows] - syz[ext_idx - (ir + 1) * n_rows]);
					}

					for (ir = j; ir < geo_const::pml_fdorder_half; ++ir)
					{
						sxy_y_tmp += pml_dc[ir] * (sxy[ext_idx + ir * n_rows] - front_sxy[i + (j - ir - 1 + geo_const::pml_fdorder_half) * dir_n_rows + k * dir_n_rows * geo_const::pml_fdorder_half]);
						szy_y_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_rows] - front_syz[i + (j - ir - 1 + geo_const::pml_fdorder_half) * dir_n_rows + k * dir_n_rows * geo_const::pml_fdorder_half]);
					}
				}
			}
			else if (j > dir_n_cols - geo_const::pml_fdorder_half - 1 && j < dir_n_cols)
			{
				if (is_pml_back)
				{
#pragma unroll
					for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
					{
						syy_y_tmp += pml_dc[ir] * (syy[ext_idx + (ir + 1) * n_rows] - syy[ext_idx - ir * n_rows]);
						sxy_y_tmp += pml_dc[ir] * (sxy[ext_idx + ir * n_rows] - sxy[ext_idx - (ir + 1) * n_rows]);
						szy_y_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_rows] - syz[ext_idx - (ir + 1) * n_rows]);
					}
				}
				else
				{
					for (ir = 0; ir < dir_n_cols - 1 - j; ++ir)
						syy_y_tmp += pml_dc[ir] * (syy[ext_idx + (ir + 1) * n_rows] - syy[ext_idx - ir * n_rows]);

					for (ir = dir_n_cols - 1 - j; ir < geo_const::pml_fdorder_half; ++ir)
						syy_y_tmp += pml_dc[ir] * (back_syy[i + (j + ir + 1 - dir_n_cols) * dir_n_rows + k * dir_n_rows * geo_const::pml_fdorder_half] - syy[ext_idx - ir * n_rows]);

					for (ir = 0; ir <= dir_n_cols - 1 - j; ++ir)
					{
						sxy_y_tmp += pml_dc[ir] * (sxy[ext_idx + ir * n_rows] - sxy[ext_idx - (ir + 1) * n_rows]);
						szy_y_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_rows] - syz[ext_idx - (ir + 1) * n_rows]);
					}

					for (ir = dir_n_cols - j; ir < geo_const::pml_fdorder_half; ++ir)
					{
						sxy_y_tmp += pml_dc[ir] * (back_sxy[i + (j + ir - dir_n_cols) * dir_n_rows + k * dir_n_rows * (geo_const::pml_fdorder_half - 1)] - sxy[ext_idx - (ir + 1) * n_rows]);
						szy_y_tmp += pml_dc[ir] * (back_syz[i + (j + ir - dir_n_cols) * dir_n_rows + k * dir_n_rows * (geo_const::pml_fdorder_half - 1)] - syz[ext_idx - (ir + 1) * n_rows]);
					}
				}
			}
			//*X
			//
#pragma unroll
			for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
			{
				sxx_x_tmp += pml_dc[ir] * (sxx[ext_idx + ir + 1] - sxx[ext_idx - ir]);
				syx_x_tmp += pml_dc[ir] * (sxy[ext_idx + ir] - sxy[ext_idx - ir - 1]);
				szx_x_tmp += pml_dc[ir] * (sxz[ext_idx + ir] - sxz[ext_idx - ir - 1]);
			}

			//*Z
			if (k >= geo_const::pml_fdorder_half && k <= dir_n_slices - geo_const::pml_fdorder_half - 1)
			{
#pragma unroll
				for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
				{
					szz_z_tmp += pml_dc[ir] * (szz[ext_idx + (ir + 1) * n_elem_slice] - szz[ext_idx - ir * n_elem_slice]);
					sxz_z_tmp += pml_dc[ir] * (sxz[ext_idx + ir * n_elem_slice] - sxz[ext_idx - (ir + 1) * n_elem_slice]);
					syz_z_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_elem_slice] - syz[ext_idx - (ir + 1) * n_elem_slice]);
				}
			}
			else if (k < geo_const::pml_fdorder_half)
			{
				if (is_pml_top)
				{
#pragma unroll
					for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
					{
						szz_z_tmp += pml_dc[ir] * (szz[ext_idx + (ir + 1) * n_elem_slice] - szz[ext_idx - ir * n_elem_slice]);
						sxz_z_tmp += pml_dc[ir] * (sxz[ext_idx + ir * n_elem_slice] - sxz[ext_idx - (ir + 1) * n_elem_slice]);
						syz_z_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_elem_slice] - syz[ext_idx - (ir + 1) * n_elem_slice]);
					}
				}
				else
				{
					for (ir = 0; ir <= k; ++ir)
						szz_z_tmp += pml_dc[ir] * (szz[ext_idx + (ir + 1) * n_elem_slice] - szz[ext_idx - ir * n_elem_slice]);

					for (ir = k + 1; ir < geo_const::pml_fdorder_half; ++ir)
					{
						szz_z_tmp += pml_dc[ir] * (szz[ext_idx + (ir + 1) * n_elem_slice] - top_szz[i + j * dir_n_rows + (k - ir - 1 + geo_const::pml_fdorder_half) * dir_n_elem_slice]);
					}

					for (ir = 0; ir < k; ++ir)
					{
						sxz_z_tmp += pml_dc[ir] * (sxz[ext_idx + ir * n_elem_slice] - sxz[ext_idx - (ir + 1) * n_elem_slice]);
						syz_z_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_elem_slice] - syz[ext_idx - (ir + 1) * n_elem_slice]);
					}

					for (ir = k; ir < geo_const::pml_fdorder_half; ++ir)
					{
						sxz_z_tmp += pml_dc[ir] * (sxz[ext_idx + ir * n_elem_slice] - top_sxz[i + j * dir_n_rows + (k - ir - 1 + geo_const::pml_fdorder_half) * dir_n_elem_slice]);
						syz_z_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_elem_slice] - top_syz[i + j * dir_n_rows + (k - ir - 1 + geo_const::pml_fdorder_half) * dir_n_elem_slice]);
					}
				}
			}
			else if (k > dir_n_slices - geo_const::pml_fdorder_half - 1 && k < dir_n_slices)
			{
				if (is_pml_bottom)
				{
#pragma unroll
					for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
					{
						szz_z_tmp += pml_dc[ir] * (szz[ext_idx + (ir + 1) * n_elem_slice] - szz[ext_idx - ir * n_elem_slice]);
						sxz_z_tmp += pml_dc[ir] * (sxz[ext_idx + ir * n_elem_slice] - sxz[ext_idx - (ir + 1) * n_elem_slice]);
						syz_z_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_elem_slice] - syz[ext_idx - (ir + 1) * n_elem_slice]);
					}
				}
				else
				{
					for (ir = 0; ir < dir_n_slices - 1 - k; ++ir)
						szz_z_tmp += pml_dc[ir] * (szz[ext_idx + (ir + 1) * n_elem_slice] - szz[ext_idx - ir * n_elem_slice]);

					for (ir = dir_n_slices - 1 - k; ir < geo_const::pml_fdorder_half; ++ir)
						szz_z_tmp += pml_dc[ir] * (bottom_szz[i + j * dir_n_rows + (k + ir + 1 - dir_n_slices) * dir_n_elem_slice] - szz[ext_idx - ir * n_elem_slice]);

					for (ir = 0; ir <= dir_n_slices - 1 - k; ++ir)
					{
						sxz_z_tmp += pml_dc[ir] * (sxz[ext_idx + ir * n_elem_slice] - sxz[ext_idx - (ir + 1) * n_elem_slice]);
						syz_z_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_elem_slice] - syz[ext_idx - (ir + 1) * n_elem_slice]);
					}

					for (ir = dir_n_slices - k; ir < geo_const::pml_fdorder_half; ++ir)
					{
						sxz_z_tmp += pml_dc[ir] * (bottom_sxz[i + j * dir_n_rows + (k + ir - dir_n_slices) * dir_n_elem_slice] - sxz[ext_idx - (ir + 1) * n_elem_slice]);
						syz_z_tmp += pml_dc[ir] * (bottom_syz[i + j * dir_n_rows + (k + ir - dir_n_slices) * dir_n_elem_slice] - syz[ext_idx - (ir + 1) * n_elem_slice]);
					}
				}
			}
			//
			sxx_x_tmp = sxx_x_tmp / d_rows;
			syy_y_tmp = syy_y_tmp / d_cols;
			szz_z_tmp = szz_z_tmp / d_slices;
			//
			sxy_y_tmp = sxy_y_tmp / d_cols;
			sxz_z_tmp = sxz_z_tmp / d_slices;
			syx_x_tmp = syx_x_tmp / d_rows;
			//
			syz_z_tmp = syz_z_tmp / d_slices;
			szx_x_tmp = szx_x_tmp / d_rows;
			szy_y_tmp = szy_y_tmp / d_cols;
			//
			dir_vx_x[reg_idx] = dir_b_x_half[reg_idx] * dir_vx_x[reg_idx] + dir_a_x_half[reg_idx] * sxx_x_tmp;
			dir_vy_x[reg_idx] = dir_b_x[reg_idx] * dir_vy_x[reg_idx] + dir_a_x[reg_idx] * syx_x_tmp;
			dir_vz_x[reg_idx] = dir_b_x[reg_idx] * dir_vz_x[reg_idx] + dir_a_x[reg_idx] * szx_x_tmp;
			//
			sxx_x_tmp += dir_vx_x[reg_idx];
			syx_x_tmp += dir_vy_x[reg_idx];
			szx_x_tmp += dir_vz_x[reg_idx];
			//
			vx[ext_idx] += dt_rho * (sxx_x_tmp + sxy_y_tmp + sxz_z_tmp);
			vy[ext_idx] += dt_rho * (syx_x_tmp + syy_y_tmp + syz_z_tmp);
			vz[ext_idx] += dt_rho * (szx_x_tmp + szy_y_tmp + szz_z_tmp);
		}

		inline __device__ void mpi_update_stress_pml_x_dir(Frame gridext, Frame griddir,
														   int ext_idx, int reg_idx,
														   int ext_i, int ext_j, int ext_k, int i, int j, int k,
														   bool is_pml_top, bool is_pml_bottom,
														   bool is_pml_front, bool is_pml_back,
														   float dt, float *pml_dc,
														   float *lambda, float *mu,
														   float *vx, float *vy, float *vz,
														   float *sxx, float *syy, float *szz,
														   float *sxy, float *sxz, float *syz,
														   float *sau,
														   float *dir_sxx_x, float *dir_sxy_x, float *dir_sxz_x,

														   float *dir_a_x,
														   float *dir_a_x_half,
														   float *dir_b_x,
														   float *dir_b_x_half,

														   float *top_vx, float *top_vy, float *top_vz,
														   float *bottom_vx, float *bottom_vy, float *bottom_vz,
														   float *front_vx, float *front_vy, float *front_vz,
														   float *back_vx, float *back_vy, float *back_vz)
		{
			set_frame_nd(gridext);
			set_frame_nd_pre(dir, griddir);
			float lambda_tmp, mu_tmp;
			float vx_x_tmp = 0.f, vy_y_tmp = 0.f, vz_z_tmp = 0.f;
			float vx_y_tmp = 0.f, vy_x_tmp = 0.f;
			float vx_z_tmp = 0.f, vz_x_tmp = 0.f;
			float vy_z_tmp = 0.f, vz_y_tmp = 0.f;
			//
			int ir;
			//*Y
			if (j >= geo_const::pml_fdorder_half && j <= dir_n_cols - geo_const::pml_fdorder_half - 1)
			{
#pragma unroll
				for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
				{
					vy_y_tmp += pml_dc[ir] * (vy[ext_idx + ir * n_rows] - vy[ext_idx - (ir + 1) * n_rows]);
					vx_y_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_rows] - vx[ext_idx - ir * n_rows]);
					vz_y_tmp += pml_dc[ir] * (vz[ext_idx + (ir + 1) * n_rows] - vz[ext_idx - ir * n_rows]);
				}
			}
			else if (j < geo_const::pml_fdorder_half)
			{
				if (is_pml_front)
				{
#pragma unroll
					for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vy_y_tmp += pml_dc[ir] * (vy[ext_idx + ir * n_rows] - vy[ext_idx - (ir + 1) * n_rows]);
						vx_y_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_rows] - vx[ext_idx - ir * n_rows]);
						vz_y_tmp += pml_dc[ir] * (vz[ext_idx + (ir + 1) * n_rows] - vz[ext_idx - ir * n_rows]);
					}
				}
				else
				{
					for (ir = 0; ir < j; ++ir)
					{
						vy_y_tmp += pml_dc[ir] * (vy[ext_idx + ir * n_rows] - vy[ext_idx - (ir + 1) * n_rows]);
					}

					for (ir = j; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vy_y_tmp += pml_dc[ir] * (vy[ext_idx + ir * n_rows] - front_vy[i + (j - ir - 1 + geo_const::pml_fdorder_half) * dir_n_rows + k * dir_n_rows * geo_const::pml_fdorder_half]);
					}

					for (ir = 0; ir <= j; ++ir)
					{
						vx_y_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_rows] - vx[ext_idx - ir * n_rows]);
						vz_y_tmp += pml_dc[ir] * (vz[ext_idx + (ir + 1) * n_rows] - vz[ext_idx - ir * n_rows]);
					}

					for (ir = j + 1; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vx_y_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_rows] - front_vx[i + (j - ir - 1 + geo_const::pml_fdorder_half) * dir_n_rows + k * dir_n_rows * (geo_const::pml_fdorder_half - 1)]);
						vz_y_tmp += pml_dc[ir] * (vz[ext_idx + (ir + 1) * n_rows] - front_vz[i + (j - ir - 1 + geo_const::pml_fdorder_half) * dir_n_rows + k * dir_n_rows * (geo_const::pml_fdorder_half - 1)]);
					}
				}
			}
			else if (j > dir_n_cols - geo_const::pml_fdorder_half - 1 && j < dir_n_cols)
			{
				if (is_pml_back)
				{
#pragma unroll
					for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vy_y_tmp += pml_dc[ir] * (vy[ext_idx + ir * n_rows] - vy[ext_idx - (ir + 1) * n_rows]);
						vx_y_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_rows] - vx[ext_idx - ir * n_rows]);
						vz_y_tmp += pml_dc[ir] * (vz[ext_idx + (ir + 1) * n_rows] - vz[ext_idx - ir * n_rows]);
					}
				}
				else
				{
					for (ir = 0; ir <= dir_n_cols - 1 - j; ++ir)
					{
						vy_y_tmp += pml_dc[ir] * (vy[ext_idx + ir * n_rows] - vy[ext_idx - (ir + 1) * n_rows]);
					}

					for (ir = dir_n_cols - j; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vy_y_tmp += pml_dc[ir] * (back_vy[i + (j + ir - dir_n_cols) * dir_n_rows + k * dir_n_rows * (geo_const::pml_fdorder_half - 1)] - vy[ext_idx - (ir + 1) * n_rows]);
					}

					for (ir = 0; ir < dir_n_cols - 1 - j; ++ir)
					{
						vx_y_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_rows] - vx[ext_idx - ir * n_rows]);
						vz_y_tmp += pml_dc[ir] * (vz[ext_idx + (ir + 1) * n_rows] - vz[ext_idx - ir * n_rows]);
					}

					for (ir = dir_n_cols - 1 - j; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vx_y_tmp += pml_dc[ir] * (back_vx[i + (j + ir + 1 - dir_n_cols) * dir_n_rows + k * dir_n_rows * geo_const::pml_fdorder_half] - vx[ext_idx - ir * n_rows]);
						vz_y_tmp += pml_dc[ir] * (back_vz[i + (j + ir + 1 - dir_n_cols) * dir_n_rows + k * dir_n_rows * geo_const::pml_fdorder_half] - vz[ext_idx - ir * n_rows]);
					}
				}
			}

			//*X
#pragma unroll
			for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
			{
				vx_x_tmp += pml_dc[ir] * (vx[ext_idx + ir] - vx[ext_idx - ir - 1]);
				vy_x_tmp += pml_dc[ir] * (vy[ext_idx + ir + 1] - vy[ext_idx - ir]);
				vz_x_tmp += pml_dc[ir] * (vz[ext_idx + ir + 1] - vz[ext_idx - ir]);
			}

			//*Z
			if (k >= geo_const::pml_fdorder_half && k <= dir_n_slices - geo_const::pml_fdorder_half - 1)
			{
#pragma unroll
				for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
				{
					vz_z_tmp += pml_dc[ir] * (vz[ext_idx + ir * n_elem_slice] - vz[ext_idx - (ir + 1) * n_elem_slice]);
					vx_z_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_elem_slice] - vx[ext_idx - ir * n_elem_slice]);
					vy_z_tmp += pml_dc[ir] * (vy[ext_idx + (ir + 1) * n_elem_slice] - vy[ext_idx - ir * n_elem_slice]);
				}
			}
			else if (k < geo_const::pml_fdorder_half)
			{

				if (is_pml_top)
				{
#pragma unroll
					for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vx_z_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_elem_slice] - vx[ext_idx - ir * n_elem_slice]);
						vy_z_tmp += pml_dc[ir] * (vy[ext_idx + (ir + 1) * n_elem_slice] - vy[ext_idx - ir * n_elem_slice]);
						vz_z_tmp += pml_dc[ir] * (vz[ext_idx + ir * n_elem_slice] - vz[ext_idx - (ir + 1) * n_elem_slice]);
					}
				}
				else
				{
					for (ir = 0; ir < k; ++ir)
					{
						vz_z_tmp += pml_dc[ir] * (vz[ext_idx + ir * n_elem_slice] - vz[ext_idx - (ir + 1) * n_elem_slice]);
					}

					for (ir = k; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vz_z_tmp += pml_dc[ir] * (vz[ext_idx + ir * n_elem_slice] - top_vz[i + j * dir_n_rows + (k - ir - 1 + geo_const::pml_fdorder_half) * dir_n_elem_slice]);
					}

					for (ir = 0; ir <= k; ++ir)
					{
						vx_z_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_elem_slice] - vx[ext_idx - ir * n_elem_slice]);
						vy_z_tmp += pml_dc[ir] * (vy[ext_idx + (ir + 1) * n_elem_slice] - vy[ext_idx - ir * n_elem_slice]);
					}

					for (ir = k + 1; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vx_z_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_elem_slice] - top_vx[i + j * dir_n_rows + (k - ir - 1 + geo_const::pml_fdorder_half) * dir_n_elem_slice]);
						vy_z_tmp += pml_dc[ir] * (vy[ext_idx + (ir + 1) * n_elem_slice] - top_vy[i + j * dir_n_rows + (k - ir - 1 + geo_const::pml_fdorder_half) * dir_n_elem_slice]);
					}
				}
			}
			else if (k > dir_n_slices - geo_const::pml_fdorder_half - 1 && k < dir_n_slices)
			{
				if (is_pml_bottom)
				{
#pragma unroll
					for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vz_z_tmp += pml_dc[ir] * (vz[ext_idx + ir * n_elem_slice] - vz[ext_idx - (ir + 1) * n_elem_slice]);
						vx_z_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_elem_slice] - vx[ext_idx - ir * n_elem_slice]);
						vy_z_tmp += pml_dc[ir] * (vy[ext_idx + (ir + 1) * n_elem_slice] - vy[ext_idx - ir * n_elem_slice]);
					}
				}
				else
				{
					for (ir = 0; ir <= dir_n_slices - 1 - k; ++ir)
					{
						vz_z_tmp += pml_dc[ir] * (vz[ext_idx + ir * n_elem_slice] - vz[ext_idx - (ir + 1) * n_elem_slice]);
					}

					for (ir = dir_n_slices - k; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vz_z_tmp += pml_dc[ir] * (bottom_vz[i + j * dir_n_rows + (k + ir - dir_n_slices) * dir_n_elem_slice] - vz[ext_idx - (ir + 1) * n_elem_slice]);
					}

					for (ir = 0; ir < dir_n_slices - 1 - k; ++ir)
					{
						vx_z_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_elem_slice] - vx[ext_idx - ir * n_elem_slice]);
						vy_z_tmp += pml_dc[ir] * (vy[ext_idx + (ir + 1) * n_elem_slice] - vy[ext_idx - ir * n_elem_slice]);
					}

					for (ir = dir_n_slices - 1 - k; ir < geo_const::pml_fdorder_half; ++ir)
					{
						vx_z_tmp += pml_dc[ir] * (bottom_vx[i + j * dir_n_rows + (k + ir + 1 - dir_n_slices) * dir_n_elem_slice] - vx[ext_idx - ir * n_elem_slice]);
						vy_z_tmp += pml_dc[ir] * (bottom_vy[i + j * dir_n_rows + (k + ir + 1 - dir_n_slices) * dir_n_elem_slice] - vy[ext_idx - ir * n_elem_slice]);
					}
				}
			}

			lambda_tmp = lambda[ext_idx];
			mu_tmp = mu[ext_idx];
			float lambda_2mu_tmp = lambda_tmp + 2.0f * mu_tmp;
			//
			vx_x_tmp = vx_x_tmp / d_rows;
			vy_y_tmp = vy_y_tmp / d_cols;
			vz_z_tmp = vz_z_tmp / d_slices;
			vy_x_tmp = vy_x_tmp / d_rows;
			vz_x_tmp = vz_x_tmp / d_rows;
			vx_y_tmp = vx_y_tmp / d_cols;
			vz_y_tmp = vz_y_tmp / d_cols;
			vx_z_tmp = vx_z_tmp / d_slices;
			vy_z_tmp = vy_z_tmp / d_slices;

			dir_sxx_x[reg_idx] = dir_b_x[reg_idx] * dir_sxx_x[reg_idx] + dir_a_x[reg_idx] * vx_x_tmp;
			dir_sxy_x[reg_idx] = dir_b_x_half[reg_idx] * dir_sxy_x[reg_idx] + dir_a_x_half[reg_idx] * vy_x_tmp;
			dir_sxz_x[reg_idx] = dir_b_x_half[reg_idx] * dir_sxz_x[reg_idx] + dir_a_x_half[reg_idx] * vz_x_tmp;
			//
			vx_x_tmp += dir_sxx_x[reg_idx];
			vy_x_tmp += dir_sxy_x[reg_idx];
			vz_x_tmp += dir_sxz_x[reg_idx];
			//
			sxx[ext_idx] += dt * (lambda_2mu_tmp * vx_x_tmp + lambda_tmp * vy_y_tmp + lambda_tmp * vz_z_tmp);
			syy[ext_idx] += dt * (lambda_tmp * vx_x_tmp + lambda_2mu_tmp * vy_y_tmp + lambda_tmp * vz_z_tmp);
			szz[ext_idx] += dt * (lambda_tmp * vx_x_tmp + lambda_tmp * vy_y_tmp + lambda_2mu_tmp * vz_z_tmp);
			sau[ext_idx] += dt * (lambda_2mu_tmp * vx_x_tmp + lambda_2mu_tmp * vy_y_tmp + lambda_2mu_tmp * vz_z_tmp);
			//
			sxy[ext_idx] += dt * mu_tmp * (vy_x_tmp + vx_y_tmp);
			sxz[ext_idx] += dt * mu_tmp * (vz_x_tmp + vx_z_tmp);
			syz[ext_idx] += dt * mu_tmp * (vz_y_tmp + vy_z_tmp);
		}
	}
}
#endif