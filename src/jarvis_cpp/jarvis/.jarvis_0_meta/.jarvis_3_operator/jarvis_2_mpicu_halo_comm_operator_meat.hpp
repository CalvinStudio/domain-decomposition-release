#ifndef JARVIS_MPICU_HALO_COMM_OPERATOR_MEAT
#define JARVIS_MPICU_HALO_COMM_OPERATOR_MEAT
//**********************************Developer*************************************
// 2022.06.13 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_2_mpicu_halo_comm_operator_bones.hpp"
namespace jarvis
{
	inline void check_is_push_back_single_mpi_stream(int _mpi_rank, int _send_rank, int _recv_rank, bool is_send, bool is_recv, int _err_flag)
	{
		if (_mpi_rank == _recv_rank)
		{
			MPI_Send(&is_recv, 1, MPI_BYTE, _send_rank, 0, MPI_COMM_WORLD);
		}
		if (_mpi_rank == _send_rank)
		{
			bool is_recv_tmp;
			MPI_Recv(&is_recv_tmp, 1, MPI_BYTE, _recv_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (is_send != is_recv_tmp)
			{
				printf("\033[41;37m[error]:\033[0m check_is_push_back_single_mpi_stream ERROR [%d]!!!\n", _err_flag);
				std::abort();
			}
		}
	}
	inline void mpiFrame::set_mpi_frame(int _mpi_rank, int _mpi_size)
	{
		mpi_rank = _mpi_rank;
		mpi_size = _mpi_size;
		if (n_rows > 0 && n_cols > 0 && n_slices > 0)
		{
			n_elem_slice = n_rows * n_cols;
			n_elem = n_elem_slice * n_slices;
			if (mpi_rank > n_elem - 1)
			{
				printf("\033[41;37m set mpi error!\033[0m");
				std::abort();
			}
			i_rows = mpi_rank % (n_elem_slice) % n_rows;
			i_cols = mpi_rank % (n_elem_slice) / n_rows;
			i_slices = mpi_rank / (n_elem_slice);
		}
		else
		{
			printf("\033[41;37m modelr rank n_elem error!\033[0m\n");
			std::abort();
		}
		if (n_elem != mpi_size)
		{
			printf("\033[41;37m mpi num config error!\033[0m\n");
			std::abort();
		}
		if (n_rows == 1)
		{
			near.is_left = false;
			near.is_right = false;
		}
		else
		{
			if (i_rows == 0)
			{
				near.is_right = false;
			}
			if (i_rows == n_rows - 1)
			{
				near.is_left = false;
			}
		}
		if (n_cols == 1)
		{
			near.is_front = false;
			near.is_back = false;
		}
		else
		{
			if (i_cols == 0)
			{
				near.is_front = false;
			}
			if (i_cols == n_cols - 1)
			{
				near.is_back = false;
			}
		}
		if (n_slices == 1)
		{
			near.is_top = false;
			near.is_bottom = false;
		}
		else
		{
			if (i_slices == 0)
			{
				near.is_top = false;
			}
			if (i_slices == n_slices - 1)
			{
				near.is_bottom = false;
			}
		}
	}
	template <typename elem_type, MemBlock mem_block>
	inline void mpicu_halo_comm_Operator_memory_Base<elem_type, mem_block>::mpi_stream_init()
	{
		jarvis_error_if_ptr_is_null(this->jarvis_mpi_cuda_stream_p, "jarvis_mpi_cuda_stream_p");
		for (int i_rank = 0; i_rank < jarvis_mpi_cuda_stream_p->mpi_frame.n_elem; ++i_rank)
		{
			int model_x = i_rank % (jarvis_mpi_cuda_stream_p->mpi_frame.n_elem_slice) % jarvis_mpi_cuda_stream_p->mpi_frame.n_rows;
			int model_y = i_rank % (jarvis_mpi_cuda_stream_p->mpi_frame.n_elem_slice) / jarvis_mpi_cuda_stream_p->mpi_frame.n_rows;
			int model_z = i_rank / (jarvis_mpi_cuda_stream_p->mpi_frame.n_elem_slice);
			if (model_z - 1 >= 0)
			{
				check_is_push_back_single_mpi_stream(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank - jarvis_mpi_cuda_stream_p->mpi_frame.n_elem_slice, send.is_top, recv.is_bottom, 5);
				if (send.is_top || recv.is_bottom)
				{
					jarvis_mpi_cuda_stream_p->single_mpi_stream_push_back(top_send, bottom_recv, jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank - jarvis_mpi_cuda_stream_p->mpi_frame.n_elem_slice, TransType::d2d, operator_id);
				}
			}
			if (model_z + 1 < jarvis_mpi_cuda_stream_p->mpi_frame.n_slices)
			{
				check_is_push_back_single_mpi_stream(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank + jarvis_mpi_cuda_stream_p->mpi_frame.n_elem_slice, send.is_bottom, recv.is_top, 6);
				if (send.is_bottom || recv.is_top)
				{
					jarvis_mpi_cuda_stream_p->single_mpi_stream_push_back(bottom_send, top_recv, jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank + jarvis_mpi_cuda_stream_p->mpi_frame.n_elem_slice, TransType::d2d, operator_id);
				}
			}
			if (model_y - 1 >= 0)
			{
				check_is_push_back_single_mpi_stream(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank - jarvis_mpi_cuda_stream_p->mpi_frame.n_rows, send.is_front, recv.is_back, 3);
				if (send.is_front || recv.is_back)
				{
					jarvis_mpi_cuda_stream_p->single_mpi_stream_push_back(front_send, back_recv, jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank - jarvis_mpi_cuda_stream_p->mpi_frame.n_rows, TransType::d2d, operator_id);
				}
			}
			if (model_y + 1 < jarvis_mpi_cuda_stream_p->mpi_frame.n_cols)
			{
				check_is_push_back_single_mpi_stream(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank + jarvis_mpi_cuda_stream_p->mpi_frame.n_rows, send.is_back, recv.is_front, 4);
				if (send.is_back || recv.is_front)
				{
					jarvis_mpi_cuda_stream_p->single_mpi_stream_push_back(back_send, front_recv, jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank + jarvis_mpi_cuda_stream_p->mpi_frame.n_rows, TransType::d2d, operator_id);
				}
			}
			if (model_x - 1 >= 0)
			{
				check_is_push_back_single_mpi_stream(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank - 1, send.is_right, recv.is_left, 1);
				if (send.is_right || recv.is_left)
				{
					jarvis_mpi_cuda_stream_p->single_mpi_stream_push_back(right_send, left_recv, jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank - 1, TransType::d2d, operator_id);
				}
			}
			if (model_x + 1 < jarvis_mpi_cuda_stream_p->mpi_frame.n_rows)
			{
				check_is_push_back_single_mpi_stream(jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank + 1, send.is_left, recv.is_right, 2);
				if (send.is_left || recv.is_right)
				{
					jarvis_mpi_cuda_stream_p->single_mpi_stream_push_back(left_send, right_recv, jarvis_mpi_cuda_stream_p->mpi_frame.mpi_rank, i_rank, i_rank + 1, TransType::d2d, operator_id);
				}
			}
		}
	}
	template <typename elem_type, MemBlock mem_block>
	inline void mpicu_halo_comm_Operator_memory_Base<elem_type, mem_block>::set_zero()
	{
		if (recv.is_top)
		{
			top_recv.device::set_zero();
		}
		if (recv.is_bottom)
		{
			bottom_recv.device::set_zero();
		}
		if (recv.is_front)
		{
			front_recv.device::set_zero();
		}
		if (recv.is_back)
		{
			back_recv.device::set_zero();
		}
		if (recv.is_right)
		{
			right_recv.device::set_zero();
		}
		if (recv.is_left)
		{
			left_recv.device::set_zero();
		}
	}
	template <typename elem_type, MemBlock mem_block>
	inline void mpicu_halo_comm_Operator_memory_Base<elem_type, mem_block>::clear()
	{
		if (recv.is_top)
		{
			top_recv.clear();
		}
		if (recv.is_bottom)
		{
			bottom_recv.clear();
		}
		if (recv.is_front)
		{
			front_recv.clear();
		}
		if (recv.is_back)
		{
			back_recv.clear();
		}
		if (recv.is_right)
		{
			right_recv.clear();
		}
		if (recv.is_left)
		{
			left_recv.clear();
		}
	}
	template <typename elem_type, MemBlock mem_block>
	inline void mpicu_halo_comm_Operator_memory_Base<elem_type, mem_block>::set_is_position()
	{
		if (this->top_recv.device::is_zero_size())
		{
			this->recv.is_top = false;
		}
		if (this->top_send.device::is_zero_size())
		{
			this->send.is_top = false;
		}
		if (this->bottom_recv.device::is_zero_size())
		{
			this->recv.is_bottom = false;
		}
		if (this->bottom_send.device::is_zero_size())
		{
			this->send.is_bottom = false;
		}
		if (this->front_recv.device::is_zero_size())
		{
			this->recv.is_front = false;
		}
		if (this->front_send.device::is_zero_size())
		{
			this->send.is_front = false;
		}
		if (this->back_recv.device::is_zero_size())
		{
			this->recv.is_back = false;
		}
		if (this->back_send.device::is_zero_size())
		{
			this->send.is_back = false;
		}
		if (this->right_recv.device::is_zero_size())
		{
			this->recv.is_right = false;
		}
		if (this->right_send.device::is_zero_size())
		{
			this->send.is_right = false;
		}
		if (this->left_recv.device::is_zero_size())
		{
			this->recv.is_left = false;
		}
		if (this->left_send.device::is_zero_size())
		{
			this->send.is_left = false;
		}
		if (!this->jarvis_mpi_cuda_stream_p->mpi_frame.near.is_top)
		{
			this->send.is_top = false;
			this->recv.is_top = false;
		}
		if (!this->jarvis_mpi_cuda_stream_p->mpi_frame.near.is_bottom)
		{
			this->send.is_bottom = false;
			this->recv.is_bottom = false;
		}
		if (!this->jarvis_mpi_cuda_stream_p->mpi_frame.near.is_front)
		{
			this->send.is_front = false;
			this->recv.is_front = false;
		}
		if (!this->jarvis_mpi_cuda_stream_p->mpi_frame.near.is_back)
		{
			this->send.is_back = false;
			this->recv.is_back = false;
		}
		if (!this->jarvis_mpi_cuda_stream_p->mpi_frame.near.is_right)
		{
			this->send.is_right = false;
			this->recv.is_right = false;
		}
		if (!this->jarvis_mpi_cuda_stream_p->mpi_frame.near.is_left)
		{
			this->send.is_left = false;
			this->recv.is_left = false;
		}
	}
	template <typename elem_type, MemBlock mem_block>
	inline void mpicu_halo_comm_Operator_Base<elem_type, mem_block>::extract_into_halo()
	{
		if (this->send.is_top)
		{
			this->top_send.extract_data_from(*box, this->jarvis_mpi_cuda_stream_p);
		}
		if (this->send.is_bottom)
		{
			this->bottom_send.extract_data_from(*box, this->jarvis_mpi_cuda_stream_p);
		}
		if (this->send.is_front)
		{
			this->front_send.extract_data_from(*box, this->jarvis_mpi_cuda_stream_p);
		}
		if (this->send.is_back)
		{
			this->back_send.extract_data_from(*box, this->jarvis_mpi_cuda_stream_p);
		}
		if (this->send.is_right)
		{
			this->right_send.extract_data_from(*box, this->jarvis_mpi_cuda_stream_p);
		}
		if (this->send.is_left)
		{
			this->left_send.extract_data_from(*box, this->jarvis_mpi_cuda_stream_p);
		}
	}
	template <typename elem_type, MemBlock mem_block>
	inline void mpicu_halo_comm_Operator_Base<elem_type, mem_block>::inject_from_halo()
	{
		if (this->send.is_top)
		{
			this->top_send.inject_data_into(*box, this->jarvis_mpi_cuda_stream_p);
		}
		if (this->send.is_bottom)
		{
			this->bottom_send.inject_data_into(*box, this->jarvis_mpi_cuda_stream_p);
		}
		if (this->send.is_front)
		{
			this->front_send.inject_data_into(*box, this->jarvis_mpi_cuda_stream_p);
		}
		if (this->send.is_back)
		{
			this->back_send.inject_data_into(*box, this->jarvis_mpi_cuda_stream_p);
		}
		if (this->send.is_right)
		{
			this->right_send.inject_data_into(*box, this->jarvis_mpi_cuda_stream_p);
		}
		if (this->send.is_left)
		{
			this->left_send.inject_data_into(*box, this->jarvis_mpi_cuda_stream_p);
		}
	}
	template <typename elem_type>
	inline void mpicu_halo_comm_Operator<elem_type, MemBlock::puppet>::set_padding_halo_grid(
		mpicu_halo_comm_Operator<elem_type, MemBlock::multiple> &_mpicu_multiple_halo_operator,
		Field<elem_type, MemType::device> *_box, Frame &_grid_base, int top_width, int bottom_width, int front_width, int back_width, int right_width, int left_width)
	{
		jarvis_error_if_ptr_is_null(this->jarvis_mpi_cuda_stream_p, "jarvis_mpi_cuda_stream_p");
		jarvis_error_if_ptr_is_null(_box, "box_p");
		this->box = _box;
		this->grid_base.copy(_grid_base);
		if (top_width > 0)
		{
			if (this->recv.is_top)
				this->top_recv.alloc(_mpicu_multiple_halo_operator.top_recv, Frame(_grid_base.n_rows, _grid_base.n_cols, top_width, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.l_rows, _grid_base.l_cols, _grid_base.r_slices - (top_width - 1) * _grid_base.d_slices));
			if (this->send.is_bottom)
				this->bottom_send.alloc(_mpicu_multiple_halo_operator.bottom_send, Frame(_grid_base.n_rows, _grid_base.n_cols, top_width, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.l_rows, _grid_base.l_cols, _grid_base.r_slices - (top_width - 1) * _grid_base.d_slices));
		}
		if (bottom_width > 0)
		{
			if (this->recv.is_bottom)
				this->bottom_recv.alloc(_mpicu_multiple_halo_operator.bottom_recv, Frame(_grid_base.n_rows, _grid_base.n_cols, bottom_width, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.l_rows, _grid_base.l_cols, _grid_base.l_slices));
			if (this->send.is_top)
				this->top_send.alloc(_mpicu_multiple_halo_operator.top_send, Frame(_grid_base.n_rows, _grid_base.n_cols, bottom_width, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.l_rows, _grid_base.l_cols, _grid_base.l_slices));
		}
		if (front_width > 0)
		{
			if (this->recv.is_front)
				this->front_recv.alloc(_mpicu_multiple_halo_operator.front_recv, Frame(_grid_base.n_rows, front_width, _grid_base.n_slices, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.l_rows, _grid_base.r_cols - (front_width - 1) * _grid_base.d_cols, _grid_base.l_slices));
			if (this->send.is_back)
				this->back_send.alloc(_mpicu_multiple_halo_operator.back_send, Frame(_grid_base.n_rows, front_width, _grid_base.n_slices, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.l_rows, _grid_base.r_cols - (front_width - 1) * _grid_base.d_cols, _grid_base.l_slices));
		}
		if (back_width > 0)
		{
			if (this->recv.is_back)
				this->back_recv.alloc(_mpicu_multiple_halo_operator.back_recv, Frame(_grid_base.n_rows, back_width, _grid_base.n_slices, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.l_rows, _grid_base.l_cols, _grid_base.l_slices));
			if (this->send.is_front)
				this->front_send.alloc(_mpicu_multiple_halo_operator.front_send, Frame(_grid_base.n_rows, back_width, _grid_base.n_slices, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.l_rows, _grid_base.l_cols, _grid_base.l_slices));
		}
		if (right_width > 0)
		{
			if (this->recv.is_right)
				this->right_recv.alloc(_mpicu_multiple_halo_operator.right_recv, Frame(right_width, _grid_base.n_cols, _grid_base.n_slices, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.r_rows - (right_width - 1) * _grid_base.d_rows, _grid_base.l_cols, _grid_base.l_slices));
			if (this->send.is_left)
				this->left_send.alloc(_mpicu_multiple_halo_operator.left_send, Frame(right_width, _grid_base.n_cols, _grid_base.n_slices, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.r_rows - (right_width - 1) * _grid_base.d_rows, _grid_base.l_cols, _grid_base.l_slices));
		}
		if (left_width > 0)
		{
			if (this->recv.is_left)
				this->left_recv.alloc(_mpicu_multiple_halo_operator.left_recv, Frame(left_width, _grid_base.n_cols, _grid_base.n_slices, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.l_rows, _grid_base.l_cols, _grid_base.l_slices));
			if (this->send.is_right)
				this->right_send.alloc(_mpicu_multiple_halo_operator.right_send, Frame(left_width, _grid_base.n_cols, _grid_base.n_slices, _grid_base.d_rows, _grid_base.d_cols, _grid_base.d_slices, _grid_base.l_rows, _grid_base.l_cols, _grid_base.l_slices));
		}
		this->set_is_position();
	}
	template <typename elem_type>
	inline void mpicu_halo_comm_Operator<elem_type, MemBlock::multiple>::alloc_halo()
	{
		jarvis_error_if_ptr_is_null(this->jarvis_mpi_cuda_stream_p, "jarvis_mpi_cuda_stream_p");
		this->set_is_position();
		if (this->send.is_top)
		{
			this->top_send.alloc();
		}
		if (this->recv.is_top)
		{
			this->top_recv.alloc();
		}
		if (this->send.is_bottom)
		{
			this->bottom_send.alloc();
		}
		if (this->recv.is_bottom)
		{
			this->bottom_recv.alloc();
		}
		if (this->send.is_front)
		{
			this->front_send.alloc();
		}
		if (this->recv.is_front)
		{
			this->front_recv.alloc();
		}
		if (this->send.is_back)
		{
			this->back_send.alloc();
		}
		if (this->recv.is_back)
		{
			this->back_recv.alloc();
		}
		if (this->send.is_right)
		{
			this->right_send.alloc();
		}
		if (this->recv.is_right)
		{
			this->right_recv.alloc();
		}
		if (this->send.is_left)
		{
			this->left_send.alloc();
		}
		if (this->recv.is_left)
		{
			this->left_recv.alloc();
		}
	}
}
#endif