#ifndef JARVIS_MPICU_HALO_COMM_OPERATOR_BONES
#define JARVIS_MPICU_HALO_COMM_OPERATOR_BONES
//**********************************Developer*************************************
// 2022.08.11 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_1_mpicu_halo_store_operator_meat.hpp"
namespace jarvis
{
	template <typename elem_type, MemBlock mem_block>
	class mpicu_halo_comm_Operator_memory_Base
	{
	private:
		using device = vector<elem_type, MemType::device, mem_block>;

	protected:
		static int operator_count;

	public:
		static jarvis_mpi_cuda_Stream<Field<elem_type, MemType::pinned_device, MemBlock::multiple>> *jarvis_mpi_cuda_stream_p;
		int operator_id;

		isPosition send;
		isPosition recv;
		//
		Field<elem_type, MemType::pinned_device, mem_block> top_send;
		Field<elem_type, MemType::pinned_device, mem_block> bottom_recv;
		//
		Field<elem_type, MemType::pinned_device, mem_block> bottom_send;
		Field<elem_type, MemType::pinned_device, mem_block> top_recv;
		//
		Field<elem_type, MemType::pinned_device, mem_block> front_send;
		Field<elem_type, MemType::pinned_device, mem_block> back_recv;
		//
		Field<elem_type, MemType::pinned_device, mem_block> back_send;
		Field<elem_type, MemType::pinned_device, mem_block> front_recv;
		//
		Field<elem_type, MemType::pinned_device, mem_block> right_send;
		Field<elem_type, MemType::pinned_device, mem_block> left_recv;
		//
		Field<elem_type, MemType::pinned_device, mem_block> left_send;
		Field<elem_type, MemType::pinned_device, mem_block> right_recv;
		//
		void mpi_stream_init();
		void set_zero();
		void clear();

	protected:
		void set_is_position();
	};
	template <typename elem_type, MemBlock mem_block>
	int mpicu_halo_comm_Operator_memory_Base<elem_type, mem_block>::operator_count = 0;
	template <typename elem_type, MemBlock mem_block>
	jarvis_mpi_cuda_Stream<Field<elem_type, MemType::pinned_device, MemBlock::multiple>> *mpicu_halo_comm_Operator_memory_Base<elem_type, mem_block>::jarvis_mpi_cuda_stream_p = nullptr;

	template <typename elem_type, MemBlock mem_block>
	class mpicu_halo_comm_Operator_Base : public mpicu_halo_comm_Operator_memory_Base<elem_type, mem_block>
	{
	public:
		Field<elem_type, MemType::device> *box;
		Frame grid_base;
		mpicu_halo_comm_Operator_Base() { ++this->operator_count, this->operator_id = -this->operator_count; }
		void extract_into_halo();
		void inject_from_halo();
	};

	template <typename elem_type, MemBlock mem_block = MemBlock::single>
	class mpicu_halo_comm_Operator
	{
	};
	template <typename elem_type>
	class mpicu_halo_comm_Operator<elem_type, MemBlock::puppet> : public mpicu_halo_comm_Operator_Base<elem_type, MemBlock::puppet>
	{
	public:
		void set_padding_halo_grid(mpicu_halo_comm_Operator<elem_type, MemBlock::multiple> &_mpicu_multiple_halo_operator,
								   Field<elem_type, MemType::device> *_box, Frame &_grid_base, int top_width, int bottom_width, int front_width, int back_width, int right_width, int left_width);
	};
	template <typename elem_type>
	class mpicu_halo_comm_Operator<elem_type, MemBlock::multiple> : public mpicu_halo_comm_Operator_memory_Base<elem_type, MemBlock::multiple>
	{
	public:
		mpicu_halo_comm_Operator() { ++this->operator_count, this->operator_id = this->operator_count; }
		void alloc_halo();
	};
}
#endif