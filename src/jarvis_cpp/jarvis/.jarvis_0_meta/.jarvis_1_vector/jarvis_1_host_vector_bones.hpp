#ifndef JARVIS_HOST_VECTOR_BONES
#define JARVIS_HOST_VECTOR_BONES
//**********************************Developer*************************************
// 2021.04.9 BY CAIWEI CALVIN CAI
//********************************************************************************
#include ".jarvis_0_vector_in.h"
namespace jarvis
{
	enum class MemType
	{
		paged,
		pinned,
		device,
		paged_device,
		pinned_device
	};
	enum class MemBlock
	{
		single,
		puppet,
		multiple
	};
	template <typename elem_type>
	class jarvis_memory
	{
	public:
		elem_type *&ptr();
		uint64_t &size();
		elem_type *ptr() const;
		uint64_t size() const;
		uint8_t type_size() const;
		bool is_zero_size() const;
		bool is_null_ptr() const;

	protected:
		void debug_error_if_size_is_zero(string _func_name) const;
		void debug_error_if_ptr_is_null(string _func_name) const;
		void debug_error_if_size_is_unequal(const jarvis_memory &_jarvis_memory_o, string _func_name) const;
		void debug_error_if_index_out_of_range(uint64_t _first_index, uint64_t _last_index, string _func_name) const;

	private:
		elem_type *elem_ptr = nullptr;
		uint64_t elem_size = 0;
	};
	template <typename elem_type, MemBlock mem_block>
	class jarvis_memory_block : public jarvis_memory<elem_type>
	{
	};
	template <typename elem_type>
	class jarvis_memory_block<elem_type, MemBlock::multiple> : public jarvis_memory<elem_type>
	{
	public:
		void push_back(jarvis_memory_block<elem_type, MemBlock::puppet> *_puppet_vector_p);

	protected:
		std::vector<jarvis_memory_block<elem_type, MemBlock::puppet> *> puppet_memory_list;
	};
	template <typename elem_type, MemBlock mem_block>
	class host_vector : public jarvis_memory_block<elem_type, mem_block>
	{
	public:
		elem_type &operator[](uint64_t _i_elem);
		elem_type operator[](uint64_t _i_elem) const;
		elem_type max() const;
		elem_type min() const;
		elem_type mean() const;
		elem_type abs_sum() const;
		elem_type abs_max() const;
		elem_type abs_mean() const;
		void set_zero();
		void fill(elem_type _elem_value, uint64_t _first_index, uint64_t _last_index);
		void save(string _file_path);
	};
	//
	template <typename elem_type, MemType mem_type = MemType::paged, MemBlock mem_block = MemBlock::single>
	class vector
	{
	};
	template <typename elem_type>
	class vector<elem_type, MemType::paged, MemBlock::single> : public host_vector<elem_type, MemBlock::single>
	{
	public:
		void alloc(uint64_t _elem_num);
		void clear();
	};
	template <typename elem_type>
	class vector<elem_type, MemType::paged, MemBlock::puppet> : public host_vector<elem_type, MemBlock::puppet>
	{
	public:
		void alloc(vector<elem_type, MemType::paged, MemBlock::multiple> &_multiple_vector, uint64_t _elem_num);
	};
	template <typename elem_type>
	class vector<elem_type, MemType::paged, MemBlock::multiple> : public host_vector<elem_type, MemBlock::multiple>
	{
	public:
		void alloc();
		void load(string _file_path);
		void clear();
	};
}
#endif