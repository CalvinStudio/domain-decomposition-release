#ifndef JARVIS_MPICU_HALO_STORE_OPERATOR_MEAT
#define JARVIS_MPICU_HALO_STORE_OPERATOR_MEAT
//**********************************Developer*************************************
// 2022.06.13 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_1_mpicu_halo_store_operator_bones.hpp"
namespace jarvis
{
	template <typename elem_type>
	inline void cu_halo_store_Operator<elem_type>::extract_into_halo()
	{
		if (have.is_top)
		{
			top.extract_data_from(*box, jarvis_mpi_cuda_stream_p);
		}
		if (have.is_bottom)
		{
			bottom.extract_data_from(*box, jarvis_mpi_cuda_stream_p);
		}
		if (have.is_front)
		{
			front.extract_data_from(*box, jarvis_mpi_cuda_stream_p);
		}
		if (have.is_back)
		{
			back.extract_data_from(*box, jarvis_mpi_cuda_stream_p);
		}
		if (have.is_right)
		{
			right.extract_data_from(*box, jarvis_mpi_cuda_stream_p);
		}
		if (have.is_left)
		{
			left.extract_data_from(*box, jarvis_mpi_cuda_stream_p);
		}
	}
	template <typename elem_type>
	inline void cu_halo_store_Operator<elem_type>::inject_from_halo()
	{
		if (have.is_top)
		{
			top.inject_data_into(*box, jarvis_mpi_cuda_stream_p);
		}
		if (have.is_bottom)
		{
			bottom.inject_data_into(*box, jarvis_mpi_cuda_stream_p);
		}
		if (have.is_front)
		{
			front.inject_data_into(*box, jarvis_mpi_cuda_stream_p);
		}
		if (have.is_back)
		{
			back.inject_data_into(*box, jarvis_mpi_cuda_stream_p);
		}
		if (have.is_right)
		{
			right.inject_data_into(*box, jarvis_mpi_cuda_stream_p);
		}
		if (have.is_left)
		{
			left.inject_data_into(*box, jarvis_mpi_cuda_stream_p);
		}
	}
	template <typename elem_type>
	inline void cu_halo_store_Operator<elem_type>::copy_halo_into_host(int it)
	{
		if (have.is_top)
		{
			top.copy_into_host(top_store_slice[it], jarvis_mpi_cuda_stream_p->d2h_stream());
		}
		if (have.is_bottom)
		{
			bottom.copy_into_host(bottom_store_slice[it], jarvis_mpi_cuda_stream_p->d2h_stream());
		}
		if (have.is_front)
		{
			front.copy_into_host(front_store_slice[it], jarvis_mpi_cuda_stream_p->d2h_stream());
		}
		if (have.is_back)
		{
			back.copy_into_host(back_store_slice[it], jarvis_mpi_cuda_stream_p->d2h_stream());
		}
		if (have.is_right)
		{
			right.copy_into_host(right_store_slice[it], jarvis_mpi_cuda_stream_p->d2h_stream());
		}
		if (have.is_left)
		{
			left.copy_into_host(left_store_slice[it], jarvis_mpi_cuda_stream_p->d2h_stream());
		}
	}
	template <typename elem_type>
	inline void cu_halo_store_Operator<elem_type>::copy_halo_from_host(int it)
	{
		if (have.is_top)
		{
			top.copy_from_host(top_store_slice[it], jarvis_mpi_cuda_stream_p->h2d_stream());
		}
		if (have.is_bottom)
		{
			bottom.copy_from_host(bottom_store_slice[it], jarvis_mpi_cuda_stream_p->h2d_stream());
		}
		if (have.is_front)
		{
			front.copy_from_host(front_store_slice[it], jarvis_mpi_cuda_stream_p->h2d_stream());
		}
		if (have.is_back)
		{
			back.copy_from_host(back_store_slice[it], jarvis_mpi_cuda_stream_p->h2d_stream());
		}
		if (have.is_right)
		{
			right.copy_from_host(right_store_slice[it], jarvis_mpi_cuda_stream_p->h2d_stream());
		}
		if (have.is_left)
		{
			left.copy_from_host(left_store_slice[it], jarvis_mpi_cuda_stream_p->h2d_stream());
		}
	}
	template <typename elem_type>
	inline void cu_halo_store_Operator<elem_type>::set_margin_halo_grid(Field<elem_type, MemType::device> *_box, Frame &_grid_base, int _slice_num, int top_width, int bottom_width, int front_width, int back_width, int right_width, int left_width)
	{
		jarvis_error_if_ptr_is_null(_box, "box");
		jarvis_error_if_ptr_is_null(jarvis_mpi_cuda_stream_p, "jarvis_mpi_cuda_stream_p");
		slice_num = _slice_num;
		box = _box;
		grid_base.copy(_grid_base);
		top.set_ndl(grid_base.n_rows, grid_base.n_cols, top_width, grid_base.d_rows, grid_base.d_cols, grid_base.d_slices, grid_base.l_rows, grid_base.l_cols, grid_base.l_slices - top_width * grid_base.d_slices);
		bottom.set_ndl(grid_base.n_rows, grid_base.n_cols, bottom_width, grid_base.d_rows, grid_base.d_cols, grid_base.d_slices, grid_base.l_rows, grid_base.l_cols, grid_base.r_slices + grid_base.d_slices);
		front.set_ndl(grid_base.n_rows, front_width, grid_base.n_slices, grid_base.d_rows, grid_base.d_cols, grid_base.d_slices, grid_base.l_rows, grid_base.l_cols - front_width * grid_base.d_cols, grid_base.l_slices);
		back.set_ndl(grid_base.n_rows, back_width, grid_base.n_slices, grid_base.d_rows, grid_base.d_cols, grid_base.d_slices, grid_base.l_rows, grid_base.r_cols + grid_base.d_cols, grid_base.l_slices);
		right.set_ndl(right_width, grid_base.n_cols, grid_base.n_slices, grid_base.d_rows, grid_base.d_cols, grid_base.d_slices, grid_base.l_rows - right_width * grid_base.d_rows, grid_base.l_cols, grid_base.l_slices);
		left.set_ndl(left_width, grid_base.n_cols, grid_base.n_slices, grid_base.d_rows, grid_base.d_cols, grid_base.d_slices, grid_base.r_rows + grid_base.d_rows, grid_base.l_cols, grid_base.l_slices);
	}
	template <typename elem_type>
	inline void cu_halo_store_Operator<elem_type>::alloc_halo()
	{
		if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_top || top.Frame::is_empty())
		{
			have.is_top = false;
		}
		if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_bottom || bottom.Frame::is_empty())
		{
			have.is_bottom = false;
		}
		if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_front || front.Frame::is_empty())
		{
			have.is_front = false;
		}
		if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_back || back.Frame::is_empty())
		{
			have.is_back = false;
		}
		if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_right || right.Frame::is_empty())
		{
			have.is_right = false;
		}
		if (jarvis_mpi_cuda_stream_p->mpi_frame.near.is_left || left.Frame::is_empty())
		{
			have.is_left = false;
		}
		if (have.is_top)
		{
			top.alloc((Frame)top);
			top_store_slice.alloc(slice_num);
			for (int i_slice = 0; i_slice < slice_num; ++i_slice)
			{
				top_store_slice[i_slice].alloc(top.n_elem);
			}
		}
		if (have.is_bottom)
		{
			bottom.alloc((Frame)bottom);
			bottom_store_slice.alloc(slice_num);
			for (int i_slice = 0; i_slice < slice_num; ++i_slice)
			{
				bottom_store_slice[i_slice].alloc(bottom.n_elem);
			}
		}
		if (have.is_front)
		{
			front.alloc((Frame)front);
			front_store_slice.alloc(slice_num);
			for (int i_slice = 0; i_slice < slice_num; ++i_slice)
			{
				front_store_slice[i_slice].alloc(front.n_elem);
			}
		}
		if (have.is_back)
		{
			back.alloc((Frame)back);
			back_store_slice.alloc(slice_num);
			for (int i_slice = 0; i_slice < slice_num; ++i_slice)
			{
				back_store_slice[i_slice].alloc(back.n_elem);
			}
		}
		if (have.is_right)
		{
			right.alloc((Frame)right);
			right_store_slice.alloc(slice_num);
			for (int i_slice = 0; i_slice < slice_num; ++i_slice)
			{
				right_store_slice[i_slice].alloc(right.n_elem);
			}
		}
		if (have.is_left)
		{
			left.alloc((Frame)left);
			left_store_slice.alloc(slice_num);
			for (int i_slice = 0; i_slice < slice_num; ++i_slice)
			{
				left_store_slice[i_slice].alloc(left.n_elem);
			}
		}
	}
}
#endif