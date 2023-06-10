#ifndef JARVIS_DEVICE_VECTOR_BONES
#define JARVIS_DEVICE_VECTOR_BONES
//**********************************Developer*************************************
// 2021.04.9 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "../.jarvis_0_math/./.jarvis_2_device_math_header_out.h"
#include "jarvis_1_host_vector_meat.hpp"
#define __jarvis_device__
namespace jarvis
{
	template <typename elem_type, MemBlock mem_block>
	class device_vector : public jarvis_memory_block<elem_type, mem_block>
	{
	public:
		void copy_from_host(vector<elem_type, MemType::paged, mem_block> &_src_host_vector);
		void copy_into_host(vector<elem_type, MemType::paged, mem_block> &_dst_host_vector);
		void copy_from_host(vector<elem_type, MemType::pinned, mem_block> &_src_host_vector, const cudaStream_t &_cuda_copy_stream);
		void copy_into_host(vector<elem_type, MemType::pinned, mem_block> &_dst_host_vector, const cudaStream_t &_cuda_copy_stream);
		void fill(elem_type _elem_value, uint64_t _first_index, uint64_t _last_index);
		void set_value(elem_type _elem_value, uint64_t _index);
		elem_type get_value(uint64_t _index);
		void set_zero();
	};
	template <typename elem_type>
	class vector<elem_type, MemType::pinned, MemBlock::single> : public host_vector<elem_type, MemBlock::single>
	{
	public:
		void alloc(uint64_t _elem_num);
		void clear();
	};
	template <typename elem_type>
	class vector<elem_type, MemType::pinned, MemBlock::puppet> : public host_vector<elem_type, MemBlock::puppet>
	{
	public:
		void alloc(vector<elem_type, MemType::pinned, MemBlock::multiple> &_multiple_vector, uint64_t _elem_num);
	};
	template <typename elem_type>
	class vector<elem_type, MemType::pinned, MemBlock::multiple> : public host_vector<elem_type, MemBlock::multiple>
	{
	public:
		void alloc();
		void load(string _file_path);
		void clear();
	};
	template <typename elem_type>
	class vector<elem_type, MemType::device, MemBlock::single> : public device_vector<elem_type, MemBlock::single>
	{
	public:
		void alloc(uint64_t _elem_num);
		void clear();
	};
	template <typename elem_type>
	class vector<elem_type, MemType::device, MemBlock::puppet> : public device_vector<elem_type, MemBlock::puppet>
	{
	public:
		void alloc(vector<elem_type, MemType::device, MemBlock::multiple> &_multiple_vector, uint64_t _elem_num);
	};
	template <typename elem_type>
	class vector<elem_type, MemType::device, MemBlock::multiple> : public device_vector<elem_type, MemBlock::multiple>
	{
	public:
		void alloc();
		void clear();
	};
}
#endif