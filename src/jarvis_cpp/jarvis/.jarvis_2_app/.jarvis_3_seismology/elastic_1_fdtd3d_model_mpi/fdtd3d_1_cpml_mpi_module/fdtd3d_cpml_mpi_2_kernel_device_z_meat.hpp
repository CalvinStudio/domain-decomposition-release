#pragma once
#ifndef _FDTD3D_CPML_MPI_KERNEL_DEVICE_Z_MEAT_HPP
#define _FDTD3D_CPML_MPI_KERNEL_DEVICE_Z_MEAT_HPP
#include "fdtd3d_cpml_mpi_2_kernel_device_y_meat.hpp"
namespace jarvis
{
	namespace cpml
	{
		inline __device__ void mpi_update_vel_pml_z_dir(Frame gridext,
														int ext_idx, int reg_idx,
														int ext_i, int ext_j, int ext_k, int i, int j, int k,
														float dt, float *pml_dc,
														float *rho,
														float *vx, float *vy, float *vz,
														float *sxx, float *syy, float *szz,
														float *sxy, float *sxz, float *syz,
														//
														float *dir_vx_x, float *dir_vx_y, float *dir_vx_z,
														float *dir_vy_x, float *dir_vy_y, float *dir_vy_z,
														float *dir_vz_x, float *dir_vz_y, float *dir_vz_z,
														//
														float *dir_a_x, float *dir_a_y, float *dir_a_z,
														float *dir_a_x_half, float *dir_a_y_half, float *dir_a_z_half,
														float *dir_b_x, float *dir_b_y, float *dir_b_z,
														float *dir_b_x_half, float *dir_b_y_half, float *dir_b_z_half,
														//
														float *left_sxx, float *left_sxy, float *left_sxz,
														float *right_sxx, float *right_sxy, float *right_sxz,
														float *front_sxy, float *front_syy, float *front_syz,
														float *back_sxy, float *back_syy, float *back_syz)
		{
			set_frame_nd(gridext);
			float sxx_x_tmp = 0.f, sxy_y_tmp = 0.f, sxz_z_tmp = 0.f,
				  syx_x_tmp = 0.f, syy_y_tmp = 0.f, syz_z_tmp = 0.f,
				  szx_x_tmp = 0.f, szy_y_tmp = 0.f, szz_z_tmp = 0.f;
			//
			float dt_rho = dt / rho[ext_idx];
			//
			int ir;
			//*X
			if (i >= geo_const::pml_fdorder_half && i <= n_rows - geo_const::pml_fdorder_half - 1)
			{
#pragma unroll
				for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
				{
					sxx_x_tmp += pml_dc[ir] * (sxx[ext_idx + ir + 1] - sxx[ext_idx - ir]);
					syx_x_tmp += pml_dc[ir] * (sxy[ext_idx + ir] - sxy[ext_idx - ir - 1]);
					szx_x_tmp += pml_dc[ir] * (sxz[ext_idx + ir] - sxz[ext_idx - ir - 1]);
				}
			}
			else if (i < geo_const::pml_fdorder_half)
			{
				for (ir = 0; ir <= i; ++ir)
					sxx_x_tmp += pml_dc[ir] * (sxx[ext_idx + ir + 1] - sxx[ext_idx - ir]);

				for (ir = i + 1; ir < geo_const::pml_fdorder_half; ++ir)
				{
					sxx_x_tmp += pml_dc[ir] * (sxx[ext_idx + ir + 1] - right_sxx[i - ir - 1 + geo_const::pml_fdorder_half + j * (geo_const::pml_fdorder_half - 1) + k * (geo_const::pml_fdorder_half - 1) * n_cols]);
				}

				for (ir = 0; ir < i; ++ir)
				{
					syx_x_tmp += pml_dc[ir] * (sxy[ext_idx + ir] - sxy[ext_idx - ir - 1]);
					szx_x_tmp += pml_dc[ir] * (sxz[ext_idx + ir] - sxz[ext_idx - ir - 1]);
				}

				for (ir = i; ir < geo_const::pml_fdorder_half; ++ir)
				{
					syx_x_tmp += pml_dc[ir] * (sxy[ext_idx + ir] - right_sxy[i - ir - 1 + geo_const::pml_fdorder_half + j * geo_const::pml_fdorder_half + k * geo_const::pml_fdorder_half * n_cols]);
					szx_x_tmp += pml_dc[ir] * (sxz[ext_idx + ir] - right_sxz[i - ir - 1 + geo_const::pml_fdorder_half + j * geo_const::pml_fdorder_half + k * geo_const::pml_fdorder_half * n_cols]);
				}
			}
			else if (i > n_rows - geo_const::pml_fdorder_half - 1 && i < n_rows)
			{
				for (ir = 0; ir < n_rows - 1 - i; ++ir)
					sxx_x_tmp += pml_dc[ir] * (sxx[ext_idx + ir + 1] - sxx[ext_idx - ir]);

				for (ir = n_rows - 1 - i; ir < geo_const::pml_fdorder_half; ++ir)
					sxx_x_tmp += pml_dc[ir] * (left_sxx[i + ir + 1 - n_rows + j * geo_const::pml_fdorder_half + k * geo_const::pml_fdorder_half * n_cols] - sxx[ext_idx - ir]);

				for (ir = 0; ir <= n_rows - 1 - i; ++ir)
				{
					syx_x_tmp += pml_dc[ir] * (sxy[ext_idx + ir] - sxy[ext_idx - ir - 1]);
					szx_x_tmp += pml_dc[ir] * (sxz[ext_idx + ir] - sxz[ext_idx - ir - 1]);
				}

				for (ir = n_rows - i; ir < geo_const::pml_fdorder_half; ++ir)
				{
					syx_x_tmp += pml_dc[ir] * (left_sxy[i + ir - n_rows + j * (geo_const::pml_fdorder_half - 1) + k * (geo_const::pml_fdorder_half - 1) * n_cols] - sxy[ext_idx - ir - 1]);
					szx_x_tmp += pml_dc[ir] * (left_sxz[i + ir - n_rows + j * (geo_const::pml_fdorder_half - 1) + k * (geo_const::pml_fdorder_half - 1) * n_cols] - sxz[ext_idx - ir - 1]);
				}
			}
			//*Y
			//
			if (j >= geo_const::pml_fdorder_half && j <= n_cols - geo_const::pml_fdorder_half - 1)
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
				for (ir = 0; ir <= j; ++ir)
					syy_y_tmp += pml_dc[ir] * (syy[ext_idx + (ir + 1) * n_rows] - syy[ext_idx - ir * n_rows]);

				for (ir = j + 1; ir < geo_const::pml_fdorder_half; ++ir)
				{
					syy_y_tmp += pml_dc[ir] * (syy[ext_idx + (ir + 1) * n_rows] - front_syy[i + (j - ir - 1 + geo_const::pml_fdorder_half) * n_rows + k * n_rows * (geo_const::pml_fdorder_half - 1)]);
				}

				for (ir = 0; ir < j; ++ir)
				{
					sxy_y_tmp += pml_dc[ir] * (sxy[ext_idx + ir * n_rows] - sxy[ext_idx - (ir + 1) * n_rows]);
					szy_y_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_rows] - syz[ext_idx - (ir + 1) * n_rows]);
				}

				for (ir = j; ir < geo_const::pml_fdorder_half; ++ir)
				{
					sxy_y_tmp += pml_dc[ir] * (sxy[ext_idx + ir * n_rows] - front_sxy[i + (j - ir - 1 + geo_const::pml_fdorder_half) * n_rows + k * n_rows * geo_const::pml_fdorder_half]);
					szy_y_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_rows] - front_syz[i + (j - ir - 1 + geo_const::pml_fdorder_half) * n_rows + k * n_rows * geo_const::pml_fdorder_half]);
				}
			}
			else if (j > n_cols - geo_const::pml_fdorder_half - 1 && j < n_cols)
			{
				for (ir = 0; ir < n_cols - 1 - j; ++ir)
					syy_y_tmp += pml_dc[ir] * (syy[ext_idx + (ir + 1) * n_rows] - syy[ext_idx - ir * n_rows]);
				for (ir = n_cols - 1 - j; ir < geo_const::pml_fdorder_half; ++ir)
					syy_y_tmp += pml_dc[ir] * (back_syy[i + (j + ir + 1 - n_cols) * n_rows + k * n_rows * geo_const::pml_fdorder_half] - syy[ext_idx - ir * n_rows]);
				for (ir = 0; ir <= n_cols - 1 - j; ++ir)
				{
					sxy_y_tmp += pml_dc[ir] * (sxy[ext_idx + ir * n_rows] - sxy[ext_idx - (ir + 1) * n_rows]);
					szy_y_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_rows] - syz[ext_idx - (ir + 1) * n_rows]);
				}
				for (ir = n_cols - j; ir < geo_const::pml_fdorder_half; ++ir)
				{
					sxy_y_tmp += pml_dc[ir] * (back_sxy[i + (j + ir - n_cols) * n_rows + k * n_rows * (geo_const::pml_fdorder_half - 1)] - sxy[ext_idx - (ir + 1) * n_rows]);
					szy_y_tmp += pml_dc[ir] * (back_syz[i + (j + ir - n_cols) * n_rows + k * n_rows * (geo_const::pml_fdorder_half - 1)] - syz[ext_idx - (ir + 1) * n_rows]);
				}
			}

			//*Z
#pragma unroll
			for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
			{
				szz_z_tmp += pml_dc[ir] * (szz[ext_idx + (ir + 1) * n_elem_slice] - szz[ext_idx - ir * n_elem_slice]);
				sxz_z_tmp += pml_dc[ir] * (sxz[ext_idx + ir * n_elem_slice] - sxz[ext_idx - (ir + 1) * n_elem_slice]);
				syz_z_tmp += pml_dc[ir] * (syz[ext_idx + ir * n_elem_slice] - syz[ext_idx - (ir + 1) * n_elem_slice]);
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
			dir_vy_y[reg_idx] = dir_b_y_half[reg_idx] * dir_vy_y[reg_idx] + dir_a_y_half[reg_idx] * syy_y_tmp;
			dir_vz_z[reg_idx] = dir_b_z_half[reg_idx] * dir_vz_z[reg_idx] + dir_a_z_half[reg_idx] * szz_z_tmp;

			dir_vx_y[reg_idx] = dir_b_y[reg_idx] * dir_vx_y[reg_idx] + dir_a_y[reg_idx] * sxy_y_tmp;
			dir_vx_z[reg_idx] = dir_b_z[reg_idx] * dir_vx_z[reg_idx] + dir_a_z[reg_idx] * sxz_z_tmp;

			dir_vy_x[reg_idx] = dir_b_x[reg_idx] * dir_vy_x[reg_idx] + dir_a_x[reg_idx] * syx_x_tmp;
			dir_vy_z[reg_idx] = dir_b_z[reg_idx] * dir_vy_z[reg_idx] + dir_a_z[reg_idx] * syz_z_tmp;

			dir_vz_x[reg_idx] = dir_b_x[reg_idx] * dir_vz_x[reg_idx] + dir_a_x[reg_idx] * szx_x_tmp;
			dir_vz_y[reg_idx] = dir_b_y[reg_idx] * dir_vz_y[reg_idx] + dir_a_y[reg_idx] * szy_y_tmp;

			sxx_x_tmp += dir_vx_x[reg_idx];
			syy_y_tmp += dir_vy_y[reg_idx];
			szz_z_tmp += dir_vz_z[reg_idx];

			sxy_y_tmp += dir_vx_y[reg_idx];
			sxz_z_tmp += dir_vx_z[reg_idx];

			syx_x_tmp += dir_vy_x[reg_idx];
			syz_z_tmp += dir_vy_z[reg_idx];

			szx_x_tmp += dir_vz_x[reg_idx];
			szy_y_tmp += dir_vz_y[reg_idx];

			vx[ext_idx] += dt_rho * (sxx_x_tmp + sxy_y_tmp + sxz_z_tmp);
			vy[ext_idx] += dt_rho * (syx_x_tmp + syy_y_tmp + syz_z_tmp);
			vz[ext_idx] += dt_rho * (szx_x_tmp + szy_y_tmp + szz_z_tmp);
		}

		inline __device__ void mpi_update_stress_pml_z_dir(Frame gridext,
														   int ext_idx, int reg_idx,
														   int ext_i, int ext_j, int ext_k, int i, int j, int k,
														   float dt, float *pml_dc,
														   float *lambda, float *mu,
														   float *vx, float *vy, float *vz,
														   float *sxx, float *syy, float *szz,
														   float *sxy, float *sxz, float *syz,
														   float *sau,
														   float *dir_sxx_x, float *dir_sxx_y, float *dir_sxx_z,
														   float *dir_syy_x, float *dir_syy_y, float *dir_syy_z,
														   float *dir_szz_x, float *dir_szz_y, float *dir_szz_z,
														   float *dir_sxy_x, float *dir_sxy_y,
														   float *dir_sxz_x, float *dir_sxz_z,
														   float *dir_syz_y, float *dir_syz_z,

														   float *dir_a_x, float *dir_a_y, float *dir_a_z,
														   float *dir_a_x_half, float *dir_a_y_half, float *dir_a_z_half,
														   float *dir_b_x, float *dir_b_y, float *dir_b_z,
														   float *dir_b_x_half, float *dir_b_y_half, float *dir_b_z_half,

														   float *left_vx, float *left_vy, float *left_vz,
														   float *right_vx, float *right_vy, float *right_vz,
														   float *front_vx, float *front_vy, float *front_vz,
														   float *back_vx, float *back_vy, float *back_vz)
		{
			set_frame_nd(gridext);

			float lambda_tmp, mu_tmp;
			float vx_x_tmp = 0.f, vy_y_tmp = 0.f, vz_z_tmp = 0.f;
			float vy_x_tmp = 0.f, vx_y_tmp = 0.f,
				  vz_x_tmp = 0.f, vx_z_tmp = 0.f,
				  vz_y_tmp = 0.f, vy_z_tmp = 0.f;
			//
			int ir;
			//*X
			if (i >= geo_const::pml_fdorder_half && i <= n_rows - geo_const::pml_fdorder_half - 1)
			{
#pragma unroll
				for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
				{
					vx_x_tmp += pml_dc[ir] * (vx[ext_idx + ir] - vx[ext_idx - ir - 1]);
					vy_x_tmp += pml_dc[ir] * (vy[ext_idx + ir + 1] - vy[ext_idx - ir]);
					vz_x_tmp += pml_dc[ir] * (vz[ext_idx + ir + 1] - vz[ext_idx - ir]);
				}
			}
			else if (i < geo_const::pml_fdorder_half)
			{
				for (ir = 0; ir < i; ++ir)
				{
					vx_x_tmp += pml_dc[ir] * (vx[ext_idx + ir] - vx[ext_idx - ir - 1]);
				}
				for (ir = i; ir < geo_const::pml_fdorder_half; ++ir)
				{
					vx_x_tmp += pml_dc[ir] * (vx[ext_idx + ir] - right_vx[(i - ir - 1 + geo_const::pml_fdorder_half) + j * geo_const::pml_fdorder_half + k * geo_const::pml_fdorder_half * n_cols]);
				}
				for (ir = 0; ir <= i; ++ir)
				{
					vy_x_tmp += pml_dc[ir] * (vy[ext_idx + ir + 1] - vy[ext_idx - ir]);
					vz_x_tmp += pml_dc[ir] * (vz[ext_idx + ir + 1] - vz[ext_idx - ir]);
				}
				for (ir = i + 1; ir < geo_const::pml_fdorder_half; ++ir)
				{
					vy_x_tmp += pml_dc[ir] * (vy[ext_idx + ir + 1] - right_vy[(i - ir - 1 + geo_const::pml_fdorder_half) + j * (geo_const::pml_fdorder_half - 1) + k * (geo_const::pml_fdorder_half - 1) * n_cols]);
					vz_x_tmp += pml_dc[ir] * (vz[ext_idx + ir + 1] - right_vz[(i - ir - 1 + geo_const::pml_fdorder_half) + j * (geo_const::pml_fdorder_half - 1) + k * (geo_const::pml_fdorder_half - 1) * n_cols]);
				}
			}
			else if (i > n_rows - geo_const::pml_fdorder_half - 1 && i < n_rows)
			{
				for (ir = 0; ir <= n_rows - 1 - i; ++ir)
				{
					vx_x_tmp += pml_dc[ir] * (vx[ext_idx + ir] - vx[ext_idx - ir - 1]);
				}
				for (ir = n_rows - i; ir < geo_const::pml_fdorder_half; ++ir)
				{
					vx_x_tmp += pml_dc[ir] * (left_vx[i + ir - n_rows + j * (geo_const::pml_fdorder_half - 1) + k * (geo_const::pml_fdorder_half - 1) * n_cols] - vx[ext_idx - ir - 1]);
				}
				for (ir = 0; ir < n_rows - 1 - i; ++ir)
				{
					vy_x_tmp += pml_dc[ir] * (vy[ext_idx + ir + 1] - vy[ext_idx - ir]);
					vz_x_tmp += pml_dc[ir] * (vz[ext_idx + ir + 1] - vz[ext_idx - ir]);
				}
				for (ir = n_rows - 1 - i; ir < geo_const::pml_fdorder_half; ++ir)
				{
					vy_x_tmp += pml_dc[ir] * (left_vy[i + ir + 1 - n_rows + j * geo_const::pml_fdorder_half + k * geo_const::pml_fdorder_half * n_cols] - vy[ext_idx - ir]);
					vz_x_tmp += pml_dc[ir] * (left_vz[i + ir + 1 - n_rows + j * geo_const::pml_fdorder_half + k * geo_const::pml_fdorder_half * n_cols] - vz[ext_idx - ir]);
				}
			}

			//*Y
			if (j >= geo_const::pml_fdorder_half && j <= n_cols - geo_const::pml_fdorder_half - 1)
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
				for (ir = 0; ir < j; ++ir)
				{
					vy_y_tmp += pml_dc[ir] * (vy[ext_idx + ir * n_rows] - vy[ext_idx - (ir + 1) * n_rows]);
				}
				for (ir = j; ir < geo_const::pml_fdorder_half; ++ir)
				{
					vy_y_tmp += pml_dc[ir] * (vy[ext_idx + ir * n_rows] - front_vy[i + (j - ir - 1 + geo_const::pml_fdorder_half) * n_rows + k * n_rows * geo_const::pml_fdorder_half]);
				}
				for (ir = 0; ir <= j; ++ir)
				{
					vx_y_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_rows] - vx[ext_idx - ir * n_rows]);
					vz_y_tmp += pml_dc[ir] * (vz[ext_idx + (ir + 1) * n_rows] - vz[ext_idx - ir * n_rows]);
				}
				for (ir = j + 1; ir < geo_const::pml_fdorder_half; ++ir)
				{
					vx_y_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_rows] - front_vx[i + (j - ir - 1 + geo_const::pml_fdorder_half) * n_rows + k * n_rows * (geo_const::pml_fdorder_half - 1)]);
					vz_y_tmp += pml_dc[ir] * (vz[ext_idx + (ir + 1) * n_rows] - front_vz[i + (j - ir - 1 + geo_const::pml_fdorder_half) * n_rows + k * n_rows * (geo_const::pml_fdorder_half - 1)]);
				}
			}
			else if (j > n_cols - geo_const::pml_fdorder_half - 1 && j < n_cols)
			{
				for (ir = 0; ir <= n_cols - 1 - j; ++ir)
				{
					vy_y_tmp += pml_dc[ir] * (vy[ext_idx + ir * n_rows] - vy[ext_idx - (ir + 1) * n_rows]);
				}
				for (ir = n_cols - j; ir < geo_const::pml_fdorder_half; ++ir)
				{
					vy_y_tmp += pml_dc[ir] * (back_vy[i + (j + ir - n_cols) * n_rows + k * n_rows * (geo_const::pml_fdorder_half - 1)] - vy[ext_idx - (ir + 1) * n_rows]);
				}
				for (ir = 0; ir < n_cols - 1 - j; ++ir)
				{
					vx_y_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_rows] - vx[ext_idx - ir * n_rows]);
					vz_y_tmp += pml_dc[ir] * (vz[ext_idx + (ir + 1) * n_rows] - vz[ext_idx - ir * n_rows]);
				}
				for (ir = n_cols - 1 - j; ir < geo_const::pml_fdorder_half; ++ir)
				{
					vx_y_tmp += pml_dc[ir] * (back_vx[i + (j + ir + 1 - n_cols) * n_rows + k * n_rows * geo_const::pml_fdorder_half] - vx[ext_idx - ir * n_rows]);
					vz_y_tmp += pml_dc[ir] * (back_vz[i + (j + ir + 1 - n_cols) * n_rows + k * n_rows * geo_const::pml_fdorder_half] - vz[ext_idx - ir * n_rows]);
				}
			}

			//*Z

#pragma unroll
			for (ir = 0; ir < geo_const::pml_fdorder_half; ++ir)
			{
				vz_z_tmp += pml_dc[ir] * (vz[ext_idx + ir * n_elem_slice] - vz[ext_idx - (ir + 1) * n_elem_slice]);
				vx_z_tmp += pml_dc[ir] * (vx[ext_idx + (ir + 1) * n_elem_slice] - vx[ext_idx - ir * n_elem_slice]);
				vy_z_tmp += pml_dc[ir] * (vy[ext_idx + (ir + 1) * n_elem_slice] - vy[ext_idx - ir * n_elem_slice]);
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
			dir_syy_y[reg_idx] = dir_b_y[reg_idx] * dir_syy_y[reg_idx] + dir_a_y[reg_idx] * vy_y_tmp;
			dir_szz_z[reg_idx] = dir_b_z[reg_idx] * dir_szz_z[reg_idx] + dir_a_z[reg_idx] * vz_z_tmp;
			dir_sxy_x[reg_idx] = dir_b_x_half[reg_idx] * dir_sxy_x[reg_idx] + dir_a_x_half[reg_idx] * vy_x_tmp;
			dir_sxy_y[reg_idx] = dir_b_y_half[reg_idx] * dir_sxy_y[reg_idx] + dir_a_y_half[reg_idx] * vx_y_tmp;
			dir_sxz_x[reg_idx] = dir_b_x_half[reg_idx] * dir_sxz_x[reg_idx] + dir_a_x_half[reg_idx] * vz_x_tmp;
			dir_sxz_z[reg_idx] = dir_b_z_half[reg_idx] * dir_sxz_z[reg_idx] + dir_a_z_half[reg_idx] * vx_z_tmp;
			dir_syz_y[reg_idx] = dir_b_y_half[reg_idx] * dir_syz_y[reg_idx] + dir_a_y_half[reg_idx] * vz_y_tmp;
			dir_syz_z[reg_idx] = dir_b_z_half[reg_idx] * dir_syz_z[reg_idx] + dir_a_z_half[reg_idx] * vy_z_tmp;
			//
			vx_x_tmp += dir_sxx_x[reg_idx];
			vy_y_tmp += dir_syy_y[reg_idx];
			vz_z_tmp += dir_szz_z[reg_idx];
			vy_x_tmp += dir_sxy_x[reg_idx];
			vx_y_tmp += dir_sxy_y[reg_idx];
			vz_x_tmp += dir_sxz_x[reg_idx];
			vx_z_tmp += dir_sxz_z[reg_idx];
			vz_y_tmp += dir_syz_y[reg_idx];
			vy_z_tmp += dir_syz_z[reg_idx];
			//
			sxx[ext_idx] += dt * (lambda_2mu_tmp * vx_x_tmp + lambda_tmp * vy_y_tmp + lambda_tmp * vz_z_tmp);
			syy[ext_idx] += dt * (lambda_tmp * vx_x_tmp + lambda_2mu_tmp * vy_y_tmp + lambda_tmp * vz_z_tmp);
			szz[ext_idx] += dt * (lambda_tmp * vx_x_tmp + lambda_tmp * vy_y_tmp + lambda_2mu_tmp * vz_z_tmp);
			sau[ext_idx] += dt * (lambda_2mu_tmp * vx_x_tmp + lambda_2mu_tmp * vy_y_tmp + lambda_2mu_tmp * vz_z_tmp);
			sxy[ext_idx] += dt * mu_tmp * (vy_x_tmp + vx_y_tmp);
			sxz[ext_idx] += dt * mu_tmp * (vz_x_tmp + vx_z_tmp);
			syz[ext_idx] += dt * mu_tmp * (vz_y_tmp + vy_z_tmp);
		}
	}
}
#endif