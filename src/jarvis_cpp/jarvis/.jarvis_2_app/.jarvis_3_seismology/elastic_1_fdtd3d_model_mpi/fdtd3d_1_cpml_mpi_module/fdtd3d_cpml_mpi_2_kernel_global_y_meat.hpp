#pragma once
#ifndef _FDTD3D_CPML_MPI_KERNEL_GLOBAL_Y_MEAT_HPP
#define _FDTD3D_CPML_MPI_KERNEL_GLOBAL_Y_MEAT_HPP
#include "fdtd3d_cpml_mpi_2_kernel_global_x_meat.hpp"
namespace jarvis
{
	namespace cpml
	{
		inline __global__ void mpi_update_vel_pml_front(Frame gridext, Frame griddir, Padding sub_pad, isPosition is_pml,
														float *rho, float dt,
														float *pml_dc,
														elasticWave::orgWave org_wave,

														float *front_vx_x, float *front_vx_y,
														float *front_vy_x, float *front_vy_y,
														float *front_vz_x, float *front_vz_y,

														float *dir_a_x, float *dir_a_y,
														float *dir_a_x_half, float *dir_a_y_half,
														float *dir_b_x, float *dir_b_y,
														float *dir_b_x_half, float *dir_b_y_half,

														float *top_sxz, float *top_syz, float *top_szz,
														float *bottom_sxz, float *bottom_syz, float *bottom_szz,

														float *right_sxx, float *right_sxy, float *right_sxz,
														float *left_sxx, float *left_sxy, float *left_sxz)
		{
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem &&
				i >= sub_pad.pad_right && i < ext_n_rows - sub_pad.pad_left)
			{
				mpi_updata_vel_pml_y_dir(gridext, griddir,
										 ext_idx, idx,
										 ext_i, ext_j, ext_k, i, j, k,
										 is_pml.is_top, is_pml.is_bottom,
										 dt, pml_dc, rho,
										 org_wave.vx, org_wave.vy, org_wave.vz,
										 org_wave.sxx, org_wave.syy, org_wave.szz,
										 org_wave.sxy, org_wave.sxz, org_wave.syz,
										 front_vx_x, front_vx_y,
										 front_vy_x, front_vy_y,
										 front_vz_x, front_vz_y,

										 dir_a_x, dir_a_y,
										 dir_a_x_half, dir_a_y_half,
										 dir_b_x, dir_b_y,
										 dir_b_x_half, dir_b_y_half,

										 top_sxz, top_syz, top_szz,
										 bottom_sxz, bottom_syz, bottom_szz,

										 right_sxx, right_sxy, right_sxz,
										 left_sxx, left_sxy, left_sxz);
			}
		}
		//
		inline __global__ void mpi_update_stress_pml_front(Frame gridext, Frame griddir, Padding sub_pad, isPosition is_pml,
														   float *lambda, float *mu,
														   float dt, float *pml_dc,
														   elasticWave::orgWave org_wave,
														   float *sau,
														   float *front_sxx_x, float *front_sxy_x, float *front_sxz_x,
														   float *front_syy_y, float *front_sxy_y, float *front_syz_y,

														   float *dir_a_x, float *dir_a_y,
														   float *dir_a_x_half, float *dir_a_y_half,
														   float *dir_b_x, float *dir_b_y,
														   float *dir_b_x_half, float *dir_b_y_half,

														   float *top_vx, float *top_vy, float *top_vz,
														   float *bottom_vx, float *bottom_vy, float *bottom_vz,
														   float *right_vx, float *right_vy, float *right_vz,
														   float *left_vx, float *left_vy, float *left_vz)
		{
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem &&
				i >= sub_pad.pad_right && i < ext_n_rows - sub_pad.pad_left)
			{
				mpi_update_stress_pml_y_dir(gridext, griddir,
											ext_idx, idx,
											ext_i, ext_j, ext_k, i, j, k,
											is_pml.is_top, is_pml.is_bottom,
											dt, pml_dc,
											lambda, mu,
											org_wave.vx, org_wave.vy, org_wave.vz,
											org_wave.sxx, org_wave.syy, org_wave.szz,
											org_wave.sxy, org_wave.sxz, org_wave.syz,
											sau,
											front_sxx_x, front_sxy_x, front_sxz_x,
											front_syy_y, front_sxy_y, front_syz_y,

											dir_a_x, dir_a_y,
											dir_a_x_half, dir_a_y_half,
											dir_b_x, dir_b_y,
											dir_b_x_half, dir_b_y_half,

											top_vx, top_vy, top_vz,
											bottom_vx, bottom_vy, bottom_vz,
											right_vx, right_vy, right_vz,
											left_vx, left_vy, left_vz);
			}
		}

		inline __global__ void mpi_update_vel_pml_back(Frame gridext, Frame griddir, Padding sub_pad, isPosition is_pml,
													   float *rho, float dt,
													   float *pml_dc,
													   elasticWave::orgWave org_wave,

													   float *back_vx_x, float *back_vx_y,
													   float *back_vy_x, float *back_vy_y,
													   float *back_vz_x, float *back_vz_y,

													   float *dir_a_x, float *dir_a_y,
													   float *dir_a_x_half, float *dir_a_y_half,
													   float *dir_b_x, float *dir_b_y,
													   float *dir_b_x_half, float *dir_b_y_half,
													   //
													   float *top_sxz, float *top_syz, float *top_szz,
													   float *bottom_sxz, float *bottom_syz, float *bottom_szz,

													   float *right_sxx, float *right_sxy, float *right_sxz,
													   float *left_sxx, float *left_sxy, float *left_sxz)
		{
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem && i >= sub_pad.pad_right && i < ext_n_rows - sub_pad.pad_left)
			{
				mpi_updata_vel_pml_y_dir(gridext, griddir,
										 ext_idx, idx,
										 ext_i, ext_j, ext_k, i, j, k,
										 is_pml.is_top, is_pml.is_bottom,
										 dt, pml_dc, rho,
										 org_wave.vx, org_wave.vy, org_wave.vz,
										 org_wave.sxx, org_wave.syy, org_wave.szz,
										 org_wave.sxy, org_wave.sxz, org_wave.syz,
										 //
										 back_vx_x, back_vx_y,
										 back_vy_x, back_vy_y,
										 back_vz_x, back_vz_y,
										 //
										 dir_a_x, dir_a_y,
										 dir_a_x_half, dir_a_y_half,
										 dir_b_x, dir_b_y,
										 dir_b_x_half, dir_b_y_half,
										 //
										 top_sxz, top_syz, top_szz,
										 bottom_sxz, bottom_syz, bottom_szz,
										 //
										 right_sxx, right_sxy, right_sxz,
										 left_sxx, left_sxy, left_sxz);
			}
		}
		inline __global__ void mpi_update_stress_pml_back(Frame gridext, Frame griddir, Padding sub_pad, isPosition is_pml,
														  float *lambda, float *miu, float dt,
														  float *pml_dc,
														  elasticWave::orgWave org_wave,
														  float *sau,
														  float *back_sxx_x, float *back_sxy_x, float *back_sxz_x,
														  float *back_syy_y, float *back_sxy_y, float *back_syz_y,

														  float *dir_a_x, float *dir_a_y,
														  float *dir_a_x_half, float *dir_a_y_half,
														  float *dir_b_x, float *dir_b_y,
														  float *dir_b_x_half, float *dir_b_y_half,

														  float *top_vx, float *top_vy, float *top_vz,
														  float *bottom_vx, float *bottom_vy, float *bottom_vz,
														  float *right_vx, float *right_vy, float *right_vz,
														  float *left_vx, float *left_vy, float *left_vz)
		{
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem && i >= sub_pad.pad_right && i < ext_n_rows - sub_pad.pad_left)
			{
				mpi_update_stress_pml_y_dir(gridext, griddir,
											ext_idx, idx,
											ext_i, ext_j, ext_k, i, j, k,
											is_pml.is_top, is_pml.is_bottom,
											dt, pml_dc,
											lambda, miu,
											org_wave.vx, org_wave.vy, org_wave.vz,
											org_wave.sxx, org_wave.syy, org_wave.szz,
											org_wave.sxy, org_wave.sxz, org_wave.syz,
											sau,
											back_sxx_x, back_sxy_x, back_sxz_x,
											back_syy_y, back_sxy_y, back_syz_y,

											dir_a_x, dir_a_y,
											dir_a_x_half, dir_a_y_half,
											dir_b_x, dir_b_y,
											dir_b_x_half, dir_b_y_half,

											top_vx, top_vy, top_vz,
											bottom_vx, bottom_vy, bottom_vz,
											right_vx, right_vy, right_vz,
											left_vx, left_vy, left_vz);
			}
		}
	}
}
#endif