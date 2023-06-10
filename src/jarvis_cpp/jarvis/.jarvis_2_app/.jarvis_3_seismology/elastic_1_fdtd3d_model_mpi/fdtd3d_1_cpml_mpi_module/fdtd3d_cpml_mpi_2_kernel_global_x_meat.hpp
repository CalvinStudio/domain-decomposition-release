#pragma once
#ifndef _FDTD3D_CPML_MPI_KERNEL_GLOBAL_X_MEAT_HPP
#define _FDTD3D_CPML_MPI_KERNEL_GLOBAL_X_MEAT_HPP
#include "fdtd3d_cpml_mpi_2_kernel_device_z_meat.hpp"
namespace jarvis
{
	namespace cpml
	{
		inline __global__ void mpi_update_vel_pml_right(Frame gridext, Frame griddir, Padding sub_pad, isPosition is_pml,
														float *rho, float dt,
														float *pml_dc,
														elasticWave::orgWave org_wave,

														float *right_vx_x,
														float *right_vy_x,
														float *right_vz_x,

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
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem)
			{
				mpi_update_vel_pml_x_dir(gridext, griddir,
										 ext_idx, idx,
										 ext_i, ext_j, ext_k, i, j, k,
										 is_pml.is_top, is_pml.is_bottom,
										 is_pml.is_front, is_pml.is_back,
										 dt, pml_dc, rho,
										 org_wave.vx, org_wave.vy, org_wave.vz,
										 org_wave.sxx, org_wave.syy, org_wave.szz,
										 org_wave.sxy, org_wave.sxz, org_wave.syz,
										 right_vx_x,
										 right_vy_x,
										 right_vz_x,

										 dir_a_x,
										 dir_a_x_half,
										 dir_b_x,
										 dir_b_x_half,

										 top_sxz, top_syz, top_szz,
										 bottom_sxz, bottom_syz, bottom_szz,

										 front_sxy, front_syy, front_syz,
										 back_sxy, back_syy, back_syz);
			}
		}
		inline __global__ void mpi_update_stress_pml_right(Frame gridext, Frame griddir,
														   isPosition is_pml,
														   float *lambda, float *mu, float dt,
														   float *pml_dc,
														   elasticWave::orgWave org_wave,
														   float *sau,
														   float *right_sxx_x,
														   float *right_sxy_x,
														   float *right_sxz_x,
														   float *dir_a_x,
														   float *dir_a_x_half,
														   float *dir_b_x,
														   float *dir_b_x_half,
														   float *top_vx, float *top_vy, float *top_vz,
														   float *bottom_vx, float *bottom_vy, float *bottom_vz,
														   float *front_vx, float *front_vy, float *front_vz,
														   float *back_vx, float *back_vy, float *back_vz)
		{
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem)
			{
				mpi_update_stress_pml_x_dir(gridext, griddir,
											ext_idx, idx,
											ext_i, ext_j, ext_k, i, j, k,
											is_pml.is_top, is_pml.is_bottom,
											is_pml.is_front, is_pml.is_back,
											dt, pml_dc,
											lambda, mu,
											org_wave.vx, org_wave.vy, org_wave.vz,
											org_wave.sxx, org_wave.syy, org_wave.szz,
											org_wave.sxy, org_wave.sxz, org_wave.syz,
											sau,
											right_sxx_x,
											right_sxy_x,
											right_sxz_x,

											dir_a_x,
											dir_a_x_half,
											dir_b_x,
											dir_b_x_half,

											top_vx, top_vy, top_vz,
											bottom_vx, bottom_vy, bottom_vz,
											front_vx, front_vy, front_vz,
											back_vx, back_vy, back_vz);
			}
		}

		//* left
		inline __global__ void mpi_update_vel_pml_left(Frame gridext, Frame griddir, Padding sub_pad, isPosition is_pml,
													   float *rho, float dt,
													   float *pml_dc,
													   elasticWave::orgWave org_wave,

													   float *left_vx_x,
													   float *left_vy_x,
													   float *left_vz_x,

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
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem)
			{
				mpi_update_vel_pml_x_dir(gridext, griddir,
										 ext_idx, idx,
										 ext_i, ext_j, ext_k, i, j, k,
										 is_pml.is_top, is_pml.is_bottom,
										 is_pml.is_front, is_pml.is_back,
										 dt, pml_dc, rho,
										 org_wave.vx, org_wave.vy, org_wave.vz,
										 org_wave.sxx, org_wave.syy, org_wave.szz,
										 org_wave.sxy, org_wave.sxz, org_wave.syz,
										 left_vx_x,
										 left_vy_x,
										 left_vz_x,

										 dir_a_x,
										 dir_a_x_half,
										 dir_b_x,
										 dir_b_x_half,

										 top_sxz, top_syz, top_szz,
										 bottom_sxz, bottom_syz, bottom_szz,

										 front_sxy, front_syy, front_syz,
										 back_sxy, back_syy, back_syz);
			}
		}
		//
		inline __global__ void mpi_update_stress_pml_left(Frame gridext, Frame griddir, isPosition is_pml,
														  float *lambda, float *mu, float dt,
														  float *pml_dc,
														  elasticWave::orgWave org_wave,
														  float *sau,
														  float *left_sxx_x,
														  float *left_sxy_x,
														  float *left_sxz_x,

														  float *dir_a_x,
														  float *dir_a_x_half,
														  float *dir_b_x,
														  float *dir_b_x_half,

														  float *top_vx, float *top_vy, float *top_vz,
														  float *bottom_vx, float *bottom_vy, float *bottom_vz,
														  float *front_vx, float *front_vy, float *front_vz,
														  float *back_vx, float *back_vy, float *back_vz)
		{
			set_cufield_grid_3d_idx(gridext, griddir);
			if (idx < griddir.n_elem)
			{
				mpi_update_stress_pml_x_dir(gridext, griddir,
											ext_idx, idx,
											ext_i, ext_j, ext_k, i, j, k,
											is_pml.is_top, is_pml.is_bottom,
											is_pml.is_front, is_pml.is_back,
											dt, pml_dc,
											lambda, mu,
											org_wave.vx, org_wave.vy, org_wave.vz,
											org_wave.sxx, org_wave.syy, org_wave.szz,
											org_wave.sxy, org_wave.sxz, org_wave.syz,
											sau,
											left_sxx_x,
											left_sxy_x,
											left_sxz_x,

											dir_a_x,
											dir_a_x_half,
											dir_b_x,
											dir_b_x_half,

											top_vx, top_vy, top_vz,
											bottom_vx, bottom_vy, bottom_vz,
											front_vx, front_vy, front_vz,
											back_vx, back_vy, back_vz);
			}
		}
	}
}
#endif