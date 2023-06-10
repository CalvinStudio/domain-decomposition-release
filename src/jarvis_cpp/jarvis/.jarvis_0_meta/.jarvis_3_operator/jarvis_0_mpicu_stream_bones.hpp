#ifndef JARVIS_MPICU_STREAM_BONES
#define JARVIS_MPICU_STREAM_BONES
//**********************************Developer*************************************
// 2022.06.13 BY CAIWEI CALVIN CAI
//********************************************************************************
#include ".jarvis_0_operator_in.h"
#include <mpi.h>
namespace jarvis
{
	struct isPosition
	{
		bool is_top = true, is_bottom = true;
		bool is_front = true, is_back = true;
		bool is_right = true, is_left = true;
		void print_info(int _mpi_rank)
		{
			printf("---------------------\n");
			printf("mpi_rank==\t%d\n", _mpi_rank);
			printf("is_top:\t%d\n", is_top);
			printf("is_bottom:\t%d\n", is_bottom);
			printf("is_front:\t%d\n", is_front);
			printf("is_back:\t%d\n", is_back);
			printf("is_right:\t%d\n", is_right);
			printf("is_left:\t%d\n", is_left);
			printf("---------------------\n");
		}
	};
	enum class TransType
	{
		h2h = 0,
		h2d,
		d2h,
		d2d
	};
	struct mpiFrame
	{
		bool is_cross_node;
		isPosition near;
		int mpi_rank;
		int mpi_size;
		int n_rows = 0;
		int n_cols = 0;
		int n_slices = 0;
		int n_elem;
		int n_elem_slice;
		int i_rows;
		int i_cols;
		int i_slices;
		void set_mpi_frame(int _mpi_rank, int _mpi_size);
	};
	template <typename field_type>
	struct mpi_single_Stream
	{
		int mpi_rank;
		int send_rank;
		int recv_rank;
		int batch_rank;
		TransType trans_type;
		field_type *send_field_p = nullptr;
		field_type *recv_field_p = nullptr;
		MPI_Request send_request;
		MPI_Request recv_request;
		MPI_Status send_status;
		MPI_Status recv_status;
	};
	template <typename field_type>
	class jarvis_mpi_cuda_Stream : public jarvis_cuda_Stream
	{
	private:
		using host = vector<float, MemType::pinned, MemBlock::multiple>;
		using device = vector<float, MemType::device, MemBlock::multiple>;

	public:
		mpiFrame mpi_frame;
		std::vector<mpi_single_Stream<field_type>> mpi_graph;

		void initialize(int _x_divide_num, int _y_divide_num, int _z_divide_num, int _mpi_rank, int _mpi_size, int _slots_per_node);
		void single_mpi_stream_push_back(field_type &_send, field_type &_recv, int _current_mpi_rank, int _send_rank, int _recv_rank, TransType _trans_type, int _batch_rank);
		void mpi_start_batch_id(int _batch_rank);
		void mpi_wait_batch_id(int _batch_rank);
		void clear();
	};
	using jarvis_global_mpicu_Stream_t = jarvis_mpi_cuda_Stream<Field<float, MemType::pinned_device, MemBlock::multiple>>;
	//******************************************************************
#define jarvis_mpi_init(argc, argv, mpi_rank, mpi_size) \
	MPI_Init(&argc, &argv);                             \
	int mpi_rank, mpi_size;                             \
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);           \
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);           \
	if (mpi_rank == 0)                                  \
		printf("[mpi process num]: %d\n", mpi_size);    \
	MPI_Barrier(MPI_COMM_WORLD)

#define jarvis_mpi_final MPI_Finalize()

#define mpi_tic(s) \
	double mpi_start##s = MPI_Wtime();

#define mpi_toc(s, mpi_rank, _string) \
	double mpi_end##s = MPI_Wtime();  \
	if (mpi_rank == 0)                \
		cout << _string << (double)(mpi_end##s - mpi_start##s) << "s" << endl;
}
#endif
