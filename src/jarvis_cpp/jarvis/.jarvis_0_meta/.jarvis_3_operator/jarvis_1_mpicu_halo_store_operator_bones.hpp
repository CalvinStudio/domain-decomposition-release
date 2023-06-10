#ifndef JARVIS_MPICU_HALO_STORE_OPERATOR_BONES
#define JARVIS_MPICU_HALO_STORE_OPERATOR_BONES
//**********************************Developer*************************************
// 2022.08.11 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_0_mpicu_stream_meat.hpp"
namespace jarvis
{
	template <typename elem_type>
	class cu_halo_store_Operator
	{
	public:
		static jarvis_mpi_cuda_Stream<Field<elem_type, MemType::pinned_device, MemBlock::multiple>> *jarvis_mpi_cuda_stream_p;

	protected:
		Field<elem_type, MemType::device> *box = nullptr;
		isPosition have;
		Frame grid_base;
		int slice_num;
		//
		Field<elem_type, MemType::device> top;
		Field<elem_type, MemType::device> bottom;
		Field<elem_type, MemType::device> front;
		Field<elem_type, MemType::device> back;
		Field<elem_type, MemType::device> right;
		Field<elem_type, MemType::device> left;
		//
		vector<vector<elem_type, MemType::pinned>> top_store_slice;
		vector<vector<elem_type, MemType::pinned>> bottom_store_slice;
		vector<vector<elem_type, MemType::pinned>> front_store_slice;
		vector<vector<elem_type, MemType::pinned>> back_store_slice;
		vector<vector<elem_type, MemType::pinned>> right_store_slice;
		vector<vector<elem_type, MemType::pinned>> left_store_slice;

	public:
		void set_margin_halo_grid(Field<elem_type, MemType::device> *_box, Frame &_grid_base, int _slice_num, int top_width, int bottom_width, int front_width, int back_width, int right_width, int left_width);
		void alloc_halo();
		void extract_into_halo();
		void inject_from_halo();
		void copy_halo_into_host(int i_slice);
		void copy_halo_from_host(int i_slice);
	};
	template <typename elem_type>
	jarvis_mpi_cuda_Stream<Field<elem_type, MemType::pinned_device, MemBlock::multiple>> *cu_halo_store_Operator<elem_type>::jarvis_mpi_cuda_stream_p = nullptr;
}
#endif