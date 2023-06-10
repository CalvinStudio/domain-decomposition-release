#ifndef JARVIS_MPICU_STREAM_MEAT
#define JARVIS_MPICU_STREAM_MEAT
//**********************************Developer*************************************
// 2022.06.13 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_0_mpicu_stream_bones.hpp"
namespace jarvis
{
	template <typename field_type>
	inline void jarvis_mpi_cuda_Stream<field_type>::initialize(int _x_divide_num, int _y_divide_num, int _z_divide_num, int _mpi_rank, int _mpi_size, int _slots_per_node)
	{
		mpi_frame.n_rows = _x_divide_num;
		mpi_frame.n_cols = _y_divide_num;
		mpi_frame.n_slices = _z_divide_num;
		if (_mpi_size / _slots_per_node > 1)
			mpi_frame.is_cross_node = true;
		else
			mpi_frame.is_cross_node = false;
		mpi_frame.set_mpi_frame(_mpi_rank, _mpi_size);
		if (mpi_frame.mpi_rank == 0)
			if (mpi_frame.is_cross_node)
				printf("cross node calculation is currently in progress!\n");
	}
	template <typename field_type>
	inline void jarvis_mpi_cuda_Stream<field_type>::single_mpi_stream_push_back(field_type &_send, field_type &_recv, int _current_mpi_rank, int _send_rank, int _recv_rank, TransType _trans_type, int _batch_rank)
	{
		if (_current_mpi_rank < 0)
		{
			printf("\n\033[41;37m[error]:\033[0mmpi_rank error!\n\n");
			std::abort();
		}
		if (_send_rank < 0 || _recv_rank < 0)
		{
			printf("\n\033[41;37m[error]:\033[0mrank error!\n\n");
			std::abort();
		}
		mpi_single_Stream<field_type> mpi_single_stream_tmp;
		int mpi_graph_size = mpi_graph.size();
		//* set mpi single stream
		mpi_single_stream_tmp.send_field_p = &_send;
		mpi_single_stream_tmp.recv_field_p = &_recv;
		mpi_single_stream_tmp.trans_type = _trans_type;
		mpi_single_stream_tmp.mpi_rank = _current_mpi_rank;
		mpi_single_stream_tmp.send_rank = _send_rank;
		mpi_single_stream_tmp.recv_rank = _recv_rank;
		mpi_single_stream_tmp.batch_rank = _batch_rank;
		//* check is equal size
		if (_current_mpi_rank == _recv_rank)
		{
			MPI_Send(&(_recv.device::size()), 1, MPI_INT, _send_rank, 0, MPI_COMM_WORLD);
		}
		if (_current_mpi_rank == _send_rank)
		{
			int _recv_n_elem = -1;
			MPI_Recv(&_recv_n_elem, 1, MPI_INT, _recv_rank, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			if (_send.device::size() != _recv_n_elem)
			{
				printf("\033[41;37m[error]:\033[0m mpiField data size mismatch!\n");
				std::abort();
			}
			else if (_send.device::is_zero_size())
			{
				printf("\033[41;37m[error]:\033[0m mpiField data size is zero!!! mpi_rank==%d\n", _send_rank);
				std::abort();
			}
		}
		//*mpi init
		switch (mpi_single_stream_tmp.trans_type)
		{
		case TransType::h2h:
			if (_current_mpi_rank == _send_rank)
			{
				MPI_Send_init(_send.host::ptr(), _send.host::size() * _send.host::type_size(), MPI_BYTE, _recv_rank, mpi_graph_size, MPI_COMM_WORLD, &mpi_single_stream_tmp.send_request);
			}
			if (_current_mpi_rank == _recv_rank)
			{
				MPI_Recv_init(_recv.host::ptr(), _recv.host::size() * _recv.host::type_size(), MPI_BYTE, _send_rank, mpi_graph_size, MPI_COMM_WORLD, &mpi_single_stream_tmp.recv_request);
			}
			break;
		case TransType::h2d:
			if (_current_mpi_rank == _send_rank)
			{
				MPI_Send_init(_send.host::ptr(), _send.host::size() * _send.host::type_size(), MPI_BYTE, _recv_rank, mpi_graph_size, MPI_COMM_WORLD, &mpi_single_stream_tmp.send_request);
			}
			if (_current_mpi_rank == _recv_rank)
			{
				MPI_Recv_init(_recv.device::ptr(), _recv.device::size() * _recv.device::type_size(), MPI_BYTE, _send_rank, mpi_graph_size, MPI_COMM_WORLD, &mpi_single_stream_tmp.recv_request);
			}
			break;
		case TransType::d2h:
			if (_current_mpi_rank == _send_rank)
			{
				MPI_Send_init(_send.device::ptr(), _send.device::size() * _send.device::type_size(), MPI_BYTE, _recv_rank, mpi_graph_size, MPI_COMM_WORLD, &mpi_single_stream_tmp.send_request);
			}
			if (_current_mpi_rank == _recv_rank)
			{
				MPI_Recv_init(_recv.host::ptr(), _recv.host::size() * _recv.host::type_size(), MPI_BYTE, _send_rank, mpi_graph_size, MPI_COMM_WORLD, &mpi_single_stream_tmp.recv_request);
			}
			break;
		case TransType::d2d:
			if (!mpi_frame.is_cross_node)
			{
				if (_current_mpi_rank == _send_rank)
				{
					MPI_Send_init(_send.device::ptr(), _send.device::size() * _send.device::type_size(), MPI_BYTE, _recv_rank, mpi_graph_size, MPI_COMM_WORLD, &mpi_single_stream_tmp.send_request);
				}
				if (_current_mpi_rank == _recv_rank)
				{
					MPI_Recv_init(_recv.device::ptr(), _recv.device::size() * _recv.device::type_size(), MPI_BYTE, _send_rank, mpi_graph_size, MPI_COMM_WORLD, &mpi_single_stream_tmp.recv_request);
				}
			}
			else
			{
				if (_current_mpi_rank == _send_rank)
				{
					MPI_Send_init(_send.host::ptr(), _send.host::size() * _send.host::type_size(), MPI_BYTE, _recv_rank, mpi_graph_size, MPI_COMM_WORLD, &mpi_single_stream_tmp.send_request);
				}
				if (_current_mpi_rank == _recv_rank)
				{
					MPI_Recv_init(_recv.device::ptr(), _recv.device::size() * _recv.device::type_size(), MPI_BYTE, _send_rank, mpi_graph_size, MPI_COMM_WORLD, &mpi_single_stream_tmp.recv_request);
				}
			}
			break;
		}
		mpi_graph.push_back(mpi_single_stream_tmp);
	}
	template <typename field_type>
	inline void jarvis_mpi_cuda_Stream<field_type>::mpi_start_batch_id(int _batch_rank)
	{
		if (mpi_frame.is_cross_node)
		{
			for (auto it = mpi_graph.begin(); it != mpi_graph.end(); ++it)
			{
				if ((*it).batch_rank == _batch_rank)
				{
					if ((*it).mpi_rank == (*it).send_rank)
					{
						(*it).send_field_p->copy_d2h(this->d2h_stream());
					}
				}
			}
			cudaStreamSynchronize(this->d2h_stream());
		}
		for (auto it = mpi_graph.begin(); it != mpi_graph.end(); ++it)
		{
			if ((*it).batch_rank == _batch_rank)
			{
				if ((*it).mpi_rank == (*it).send_rank)
				{
					MPI_Start(&((*it).send_request));
				}
				if ((*it).mpi_rank == (*it).recv_rank)
				{
					MPI_Start(&((*it).recv_request));
				}
			}
		}
	}
	template <typename field_type>
	inline void jarvis_mpi_cuda_Stream<field_type>::mpi_wait_batch_id(int _batch_rank)
	{
		for (auto it = mpi_graph.begin(); it != mpi_graph.end(); ++it)
		{
			if ((*it).batch_rank == _batch_rank)
			{
				if ((*it).mpi_rank == (*it).send_rank)
				{
					MPI_Wait(&((*it).send_request), &((*it).send_status));
				}
				if ((*it).mpi_rank == (*it).recv_rank)
				{
					MPI_Wait(&((*it).recv_request), &((*it).recv_status));
				}
			}
		}
	}
	template <typename field_type>
	inline void jarvis_mpi_cuda_Stream<field_type>::clear()
	{
		mpi_graph.clear();
	}
}
#endif