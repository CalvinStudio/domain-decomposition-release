#pragma once
#ifndef _FDTD3D_CPML_MPI_KERNEL_GLOBAL_Z_MEAT_HPP
#define _FDTD3D_CPML_MPI_KERNEL_GLOBAL_Z_MEAT_HPP
#include "fdtd3d_cpml_mpi_2_kernel_global_y_meat.hpp"
namespace jarvis
{
	namespace cpml
	{
		inline __global__ void mpi_update_vel_pml_top(Frame gridext, Frame griddir, Padding sub_pad,
													  float *rho, float dt,
													  float *pml_dc,
													  elasticWave::orgWave org_wave,

													  float *top_vx_x, float *top_vx_y, float *top_vx_z,
													  float *top_vy_x, float *top_vy_y, float *top_vy_z,
													  float *top_vz_x, float *top_vz_y, float *top_vz_z,

													  float *dir_a_x, float *dir_a_y, float *dir_a_z,
													  float *dir_a_x_half, float *dir_a_y_half, float *dir_a_z_half,
													  float *dir_b_x, float *dir_b_y, float *dir_b_z,
													  float *dir_b_x_half, float *dir_b_y_half, float *dir_b_z_half,

													  float *left_sxx, float *left_sxy, float *left_sxz,
													  float *right_sxx, float *right_sxy, float *right_sxz,

													  float *front_sxy, float *front_syy, float *front_syz,
													  float *back_sxy, float *back_syy, float *back_syz)
		{
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem &&
				i >= sub_pad.pad_right && i < ext_n_rows - sub_pad.pad_left &&
				j >= sub_pad.pad_front && j < gridext.n_cols - sub_pad.pad_back)
			{
				mpi_update_vel_pml_z_dir(gridext,
										 ext_idx, idx,
										 ext_i, ext_j, ext_k, i, j, k,
										 dt, pml_dc,
										 rho,
										 org_wave.vx, org_wave.vy, org_wave.vz,
										 org_wave.sxx, org_wave.syy, org_wave.szz,
										 org_wave.sxy, org_wave.sxz, org_wave.syz,
										 //
										 top_vx_x, top_vx_y, top_vx_z,
										 top_vy_x, top_vy_y, top_vy_z,
										 top_vz_x, top_vz_y, top_vz_z,

										 dir_a_x, dir_a_y, dir_a_z,
										 dir_a_x_half, dir_a_y_half, dir_a_z_half,
										 dir_b_x, dir_b_y, dir_b_z,
										 dir_b_x_half, dir_b_y_half, dir_b_z_half,

										 left_sxx, left_sxy, left_sxz,
										 right_sxx, right_sxy, right_sxz,

										 front_sxy, front_syy, front_syz,
										 back_sxy, back_syy, back_syz);
			}
		}
		inline __global__ void mpi_update_stress_pml_top(Frame gridext, Frame griddir, Padding sub_pad,
														 float *lambda, float *mu, float dt,
														 float *pml_dc,
														 elasticWave::orgWave org_wave,
														 float *sau,
														 float *top_sxx_x, float *top_sxx_y, float *top_sxx_z,
														 float *top_syy_x, float *top_syy_y, float *top_syy_z,
														 float *top_szz_x, float *top_szz_y, float *top_szz_z,
														 float *top_sxy_x, float *top_sxy_y,
														 float *top_sxz_x, float *top_sxz_z,
														 float *top_syz_y, float *top_syz_z,

														 float *dir_a_x, float *dir_a_y, float *dir_a_z,
														 float *dir_a_x_half, float *dir_a_y_half, float *dir_a_z_half,

														 float *dir_b_x, float *dir_b_y, float *dir_b_z,
														 float *dir_b_x_half, float *dir_b_y_half, float *dir_b_z_half,

														 float *left_vx, float *left_vy, float *left_vz,
														 float *right_vx, float *right_vy, float *right_vz,
														 float *front_vx, float *front_vy, float *front_vz,
														 float *back_vx, float *back_vy, float *back_vz)
		{
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem &&
				i >= sub_pad.pad_right && i < ext_n_rows - sub_pad.pad_left &&
				j >= sub_pad.pad_front && j < gridext.n_cols - sub_pad.pad_back)
			{
				mpi_update_stress_pml_z_dir(gridext,
											ext_idx, idx,
											ext_i, ext_j, ext_k, i, j, k,
											dt, pml_dc,
											lambda, mu,
											org_wave.vx, org_wave.vy, org_wave.vz,
											org_wave.sxx, org_wave.syy, org_wave.szz,
											org_wave.sxy, org_wave.sxz, org_wave.syz,
											sau,
											top_sxx_x, top_sxx_y, top_sxx_z,
											top_syy_x, top_syy_y, top_syy_z,
											top_szz_x, top_szz_y, top_szz_z,
											top_sxy_x, top_sxy_y,
											top_sxz_x, top_sxz_z,
											top_syz_y, top_syz_z,

											dir_a_x, dir_a_y, dir_a_z,
											dir_a_x_half, dir_a_y_half, dir_a_z_half,
											dir_b_x, dir_b_y, dir_b_z,
											dir_b_x_half, dir_b_y_half, dir_b_z_half,

											left_vx, left_vy, left_vz,
											right_vx, right_vy, right_vz,
											front_vx, front_vy, front_vz,
											back_vx, back_vy, back_vz);
			}
		}

		// bottom
		inline __global__ void mpi_update_vel_pml_bottom(Frame gridext, Frame griddir, Padding sub_pad,
														 float *rho, float dt,
														 float *pml_dc,
														 elasticWave::orgWave org_wave,

														 float *bottom_vx_x, float *bottom_vx_y, float *bottom_vx_z,
														 float *bottom_vy_x, float *bottom_vy_y, float *bottom_vy_z,
														 float *bottom_vz_x, float *bottom_vz_y, float *bottom_vz_z,

														 float *dir_a_x, float *dir_a_y, float *dir_a_z,
														 float *dir_a_x_half, float *dir_a_y_half, float *dir_a_z_half,
														 float *dir_b_x, float *dir_b_y, float *dir_b_z,
														 float *dir_b_x_half, float *dir_b_y_half, float *dir_b_z_half,

														 float *left_sxx, float *left_sxy, float *left_sxz,
														 float *right_sxx, float *right_sxy, float *right_sxz,

														 float *front_sxy, float *front_syy, float *front_syz,
														 float *back_sxy, float *back_syy, float *back_syz)
		{
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem &&
				i >= sub_pad.pad_right && i < ext_n_rows - sub_pad.pad_left &&
				j >= sub_pad.pad_front && j < gridext.n_cols - sub_pad.pad_back)
			{
				mpi_update_vel_pml_z_dir(gridext,
										 ext_idx, idx,
										 ext_i, ext_j, ext_k, i, j, k,
										 dt, pml_dc,
										 rho,
										 org_wave.vx, org_wave.vy, org_wave.vz,
										 org_wave.sxx, org_wave.syy, org_wave.szz,
										 org_wave.sxy, org_wave.sxz, org_wave.syz,

										 bottom_vx_x, bottom_vx_y, bottom_vx_z,
										 bottom_vy_x, bottom_vy_y, bottom_vy_z,
										 bottom_vz_x, bottom_vz_y, bottom_vz_z,

										 dir_a_x, dir_a_y, dir_a_z,
										 dir_a_x_half, dir_a_y_half, dir_a_z_half,
										 dir_b_x, dir_b_y, dir_b_z,
										 dir_b_x_half, dir_b_y_half, dir_b_z_half,

										 left_sxx, left_sxy, left_sxz,
										 right_sxx, right_sxy, right_sxz,

										 front_sxy, front_syy, front_syz,
										 back_sxy, back_syy, back_syz);
			}
		}

		inline __global__ void mpi_update_stress_pml_bottom(Frame gridext, Frame griddir, Padding sub_pad,
															float *lambda, float *mu, float dt,
															float *pml_dc,
															elasticWave::orgWave org_wave,
															float *sau,
															float *bottom_sxx_x, float *bottom_sxx_y, float *bottom_sxx_z,
															float *bottom_syy_x, float *bottom_syy_y, float *bottom_syy_z,
															float *bottom_szz_x, float *bottom_szz_y, float *bottom_szz_z,
															float *bottom_sxy_x, float *bottom_sxy_y,
															float *bottom_sxz_x, float *bottom_sxz_z,
															float *bottom_syz_y, float *bottom_syz_z,
															float *dir_a_x, float *dir_a_y, float *dir_a_z,
															float *dir_a_x_half, float *dir_a_y_half, float *dir_a_z_half,
															float *dir_b_x, float *dir_b_y, float *dir_b_z,
															float *dir_b_x_half, float *dir_b_y_half, float *dir_b_z_half,
															float *left_vx, float *left_vy, float *left_vz,
															float *right_vx, float *right_vy, float *right_vz,
															float *front_vx, float *front_vy, float *front_vz,
															float *back_vx, float *back_vy, float *back_vz)
		{
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem &&
				i >= sub_pad.pad_right && i < ext_n_rows - sub_pad.pad_left &&
				j >= sub_pad.pad_front && j < gridext.n_cols - sub_pad.pad_back)
			{
				mpi_update_stress_pml_z_dir(gridext,
											ext_idx, idx,
											ext_i, ext_j, ext_k, i, j, k,
											dt, pml_dc,
											lambda, mu,
											org_wave.vx, org_wave.vy, org_wave.vz,
											org_wave.sxx, org_wave.syy, org_wave.szz,
											org_wave.sxy, org_wave.sxz, org_wave.syz,
											sau,
											bottom_sxx_x, bottom_sxx_y, bottom_sxx_z,
											bottom_syy_x, bottom_syy_y, bottom_syy_z,
											bottom_szz_x, bottom_szz_y, bottom_szz_z,
											bottom_sxy_x, bottom_sxy_y,
											bottom_sxz_x, bottom_sxz_z,
											bottom_syz_y, bottom_syz_z,

											dir_a_x, dir_a_y, dir_a_z,
											dir_a_x_half, dir_a_y_half, dir_a_z_half,
											dir_b_x, dir_b_y, dir_b_z,
											dir_b_x_half, dir_b_y_half, dir_b_z_half,

											left_vx, left_vy, left_vz,
											right_vx, right_vy, right_vz,
											front_vx, front_vy, front_vz,
											back_vx, back_vy, back_vz);
			}
		}
	}
}
#endif