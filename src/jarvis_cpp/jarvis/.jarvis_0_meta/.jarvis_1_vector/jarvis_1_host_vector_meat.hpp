#ifndef JARVIS_HOST_VECTOR_MEAT
#define JARVIS_HOST_VECTOR_MEAT
//**********************************Developer*************************************
// 2021.04.9 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_1_host_vector_bones.hpp"
namespace jarvis
{
	template <typename elem_type>
	inline uint64_t &jarvis_memory<elem_type>::size()
	{
		return elem_size;
	}
	template <typename elem_type>
	inline elem_type *&jarvis_memory<elem_type>::ptr()
	{
		return elem_ptr;
	}
	template <typename elem_type>
	inline uint64_t jarvis_memory<elem_type>::size() const
	{
		return elem_size;
	}
	template <typename elem_type>
	inline elem_type *jarvis_memory<elem_type>::ptr() const
	{
		return elem_ptr;
	}
	template <typename elem_type>
	inline bool jarvis_memory<elem_type>::is_zero_size() const
	{
		if (elem_size == 0)
			return true;
		else
			return false;
	}
	template <typename elem_type>
	inline bool jarvis_memory<elem_type>::is_null_ptr() const
	{
		if (elem_ptr == nullptr)
			return true;
		else
			return false;
	}
	template <typename elem_type>
	inline uint8_t jarvis_memory<elem_type>::type_size() const
	{
		return sizeof(elem_type);
	}
	template <typename elem_type>
	inline void jarvis_memory<elem_type>::debug_error_if_size_is_zero(string _func_name) const
	{
		if (is_zero_size())
		{
			std::cout << _func_name + "\033[41;37m[error]:\033[0m:size is zero!" << std::endl;
			std::abort();
		}
	}
	template <typename elem_type>
	inline void jarvis_memory<elem_type>::debug_error_if_ptr_is_null(string _func_name) const
	{
		if (is_null_ptr())
		{
			std::cout << _func_name + "\033[41;37m[error]:\033[0m:prt is null!" << std::endl;
			std::abort();
		}
	}
	template <typename elem_type>
	inline void jarvis_memory<elem_type>::debug_error_if_index_out_of_range(uint64_t _first_index, uint64_t _last_index, string _func_name) const
	{
		if (_first_index < 0 || _last_index >= size() || _first_index > _last_index)
		{
			std::cout << _func_name + "\033[41;37m[error]:\033[0m:index_out_of_range!" << std::endl;
			std::abort();
		}
	}
	template <typename elem_type>
	inline void jarvis_memory<elem_type>::debug_error_if_size_is_unequal(const jarvis_memory &_jarvis_memory_o, string _func_name) const
	{
		if (size() != _jarvis_memory_o.size())
		{
			std::cout << _func_name + "\033[41;37m[error]:\033[0m:size_is_unequal!" << std::endl;
			std::abort();
		}
	}
	template <typename elem_type>
	inline void jarvis_memory_block<elem_type, MemBlock::multiple>::push_back(jarvis_memory_block<elem_type, MemBlock::puppet> *_puppet_vector_p)
	{
		puppet_memory_list.push_back(_puppet_vector_p);
		jarvis_memory<elem_type>::size() += _puppet_vector_p->size();
	}
	template <typename elem_type, MemBlock mem_block>
	inline elem_type &host_vector<elem_type, mem_block>::operator[](uint64_t _i_elem)
	{
		return *(jarvis_memory<elem_type>::ptr() + _i_elem);
	}
	template <typename elem_type, MemBlock mem_block>
	inline elem_type host_vector<elem_type, mem_block>::operator[](uint64_t _i_elem) const
	{
		return *(jarvis_memory<elem_type>::ptr() + _i_elem);
	}
	template <typename elem_type, MemBlock mem_block>
	inline elem_type host_vector<elem_type, mem_block>::max() const
	{
		return jarvis::max(jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size());
	}
	template <typename elem_type, MemBlock mem_block>
	inline elem_type host_vector<elem_type, mem_block>::min() const
	{
		return jarvis::min(jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size());
	}
	template <typename elem_type, MemBlock mem_block>
	inline elem_type host_vector<elem_type, mem_block>::abs_sum() const
	{
		return jarvis::abs_sum(jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size());
	}
	template <typename elem_type, MemBlock mem_block>
	inline elem_type host_vector<elem_type, mem_block>::abs_max() const
	{
		return jarvis::abs_max(jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size());
	}
	template <typename elem_type, MemBlock mem_block>
	inline elem_type host_vector<elem_type, mem_block>::abs_mean() const
	{
		return jarvis::abs_mean(jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size());
	}
	template <typename elem_type, MemBlock mem_block>
	inline elem_type host_vector<elem_type, mem_block>::mean() const
	{
		return jarvis::mean(jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size());
	}
	template <typename elem_type, MemBlock mem_block>
	inline void host_vector<elem_type, mem_block>::set_zero()
	{
		jarvis_memory<elem_type>::debug_error_if_ptr_is_null("jarvis::host_vector::set_zero()");
		std::memset(jarvis_memory<elem_type>::ptr(), 0, jarvis_memory<elem_type>::size() * sizeof(elem_type));
	}
	template <typename elem_type, MemBlock mem_block>
	inline void host_vector<elem_type, mem_block>::fill(elem_type _elem_value, uint64_t _first_index, uint64_t _last_index)
	{
		jarvis_memory<elem_type>::debug_error_if_index_out_of_range(_first_index, _last_index, "fill()");
		for (uint64_t _index = _first_index; _index <= _last_index; ++_index)
		{
			*(jarvis_memory<elem_type>::ptr() + _index) = _elem_value;
		}
	}
	template <typename elem_type, MemBlock mem_block>
	inline void host_vector<elem_type, mem_block>::save(string _file_path)
	{
		jarvis_memory<elem_type>::debug_error_if_ptr_is_null("jarvis::vector::save()");
		_file_path += ".jvs";
		FILE *fp = fopen(_file_path.c_str(), "wb");
		fwrite(&jarvis_memory<elem_type>::size(), 8, 1, fp);
		fwrite(jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size() * sizeof(elem_type), 1, fp);
		fclose(fp);
		printf("\033[45;30m[save jvs-format file]:\033[0m%s\n", _file_path.c_str());
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::paged, MemBlock::single>::alloc(uint64_t _elem_num)
	{
		jarvis_memory<elem_type>::size() = _elem_num;
		jarvis_memory<elem_type>::ptr() = (elem_type *)calloc(jarvis_memory<elem_type>::size(), sizeof(elem_type));
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::paged, MemBlock::single>::clear()
	{
		jarvis_memory<elem_type>::debug_error_if_ptr_is_null("jarvis::vector::clear()");
		free(jarvis_memory<elem_type>::ptr());
		jarvis_memory<elem_type>::ptr() = nullptr;
		jarvis_memory<elem_type>::size() = 0;
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::paged, MemBlock::puppet>::alloc(vector<elem_type, MemType::paged, MemBlock::multiple> &_multiple_vector, uint64_t _elem_num)
	{
		jarvis_memory<elem_type>::size() = _elem_num;
		_multiple_vector.push_back(this);
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::paged, MemBlock::multiple>::alloc()
	{
		jarvis_memory<elem_type>::debug_error_if_size_is_zero("jarvis::vector::alloc()");
		jarvis_memory<elem_type>::ptr() = (elem_type *)calloc(jarvis_memory<elem_type>::size(), sizeof(elem_type));
		uint64_t offset = 0;
		for (auto it = jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.begin();
			 it != jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.end(); ++it)
		{
			(*it)->ptr() = jarvis_memory<elem_type>::ptr() + offset;
			offset += (*it)->size();
		}
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::paged, MemBlock::multiple>::load(string _file_path)
	{
		FILE *fp;
		fp = fopen(_file_path.c_str(), "rb");
		if (!fp)
		{
			std::cout << "load():\033[41;37m[error]:\033[0mfile open error!" << std::endl;
			std::abort();
		}
		fread(&jarvis_memory<elem_type>::size(), 8, 1, fp);
		alloc();
		fread(jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size() * sizeof(elem_type), 1, fp);
		fclose(fp);
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::paged, MemBlock::multiple>::clear()
	{
		jarvis_memory<elem_type>::debug_error_if_ptr_is_null("jarvis::host_vector::clear()");
		free(jarvis_memory<elem_type>::ptr());
		jarvis_memory<elem_type>::ptr() = nullptr;
		jarvis_memory<elem_type>::size() = 0;
		for (auto it = jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.begin();
			 it != jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.end(); ++it)
		{
			(*it)->ptr() = nullptr;
			(*it)->size() = 0;
		}
	}
}
#endif