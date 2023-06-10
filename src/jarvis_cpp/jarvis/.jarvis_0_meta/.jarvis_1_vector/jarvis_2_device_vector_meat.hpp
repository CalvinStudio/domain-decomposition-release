#ifndef JARVIS_DEVICE_VECTOR_MEAT
#define JARVIS_DEVICE_VECTOR_MEAT
//**********************************Developer*************************************
// 2021.04.9 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_2_device_vector_bones.hpp"
namespace jarvis
{
	template <typename elem_type>
	inline __global__ void _jarvis_fill_cuda_kernel(elem_type *elem_ptr, uint64_t _elem_num, elem_type _elem_value)
	{
		int idx = threadIdx.x + blockDim.x * blockIdx.x;
		if (idx < _elem_num)
			elem_ptr[idx] = _elem_value;
	}
	template <typename elem_type>
	inline __global__ void _jarvis_set_value_cuda_kernel(elem_type *elem_ptr, elem_type _elem_value) { *elem_ptr = _elem_value; }
	//
	template <typename elem_type, MemBlock mem_block>
	inline void device_vector<elem_type, mem_block>::set_zero()
	{
		jarvis_memory<elem_type>::debug_error_if_ptr_is_null("jarvis::vector::set_zero()");
		jarvis_handle_cuda_error(cudaMemset(jarvis_memory<elem_type>::ptr(), 0, jarvis_memory<elem_type>::size() * sizeof(elem_type)));
	}
	template <typename elem_type, MemBlock mem_block>
	inline void device_vector<elem_type, mem_block>::copy_from_host(vector<elem_type, MemType::paged, mem_block> &_src_host_vector)
	{
		jarvis_memory<elem_type>::debug_error_if_size_is_unequal(_src_host_vector, "copy_from_host()");
		jarvis_handle_cuda_error(cudaMemcpy(jarvis_memory<elem_type>::ptr(), _src_host_vector.ptr(), jarvis_memory<elem_type>::size() * sizeof(elem_type), cudaMemcpyHostToDevice));
	}
	template <typename elem_type, MemBlock mem_block>
	inline void device_vector<elem_type, mem_block>::copy_from_host(vector<elem_type, MemType::pinned, mem_block> &_src_host_vector, const cudaStream_t &_cuda_copy_stream)
	{
		jarvis_memory<elem_type>::debug_error_if_size_is_unequal(_src_host_vector, "copy_from_host()");
		jarvis_handle_cuda_error(cudaMemcpyAsync(jarvis_memory<elem_type>::ptr(), _src_host_vector.ptr(), jarvis_memory<elem_type>::size() * sizeof(elem_type), cudaMemcpyHostToDevice, _cuda_copy_stream));
	}
	template <typename elem_type, MemBlock mem_block>
	inline void device_vector<elem_type, mem_block>::copy_into_host(vector<elem_type, MemType::paged, mem_block> &_dst_host_vector)
	{
		jarvis_memory<elem_type>::debug_error_if_size_is_unequal(_dst_host_vector, "copy_into_host()");
		jarvis_handle_cuda_error(cudaMemcpy(_dst_host_vector.ptr(), jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size() * sizeof(elem_type), cudaMemcpyDeviceToHost));
	}
	template <typename elem_type, MemBlock mem_block>
	inline void device_vector<elem_type, mem_block>::copy_into_host(vector<elem_type, MemType::pinned, mem_block> &_dst_host_vector, const cudaStream_t &_cuda_copy_stream)
	{
		jarvis_memory<elem_type>::debug_error_if_size_is_unequal(_dst_host_vector, "copy_into_host()");
		jarvis_handle_cuda_error(cudaMemcpyAsync(_dst_host_vector.ptr(), jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size() * sizeof(elem_type), cudaMemcpyDeviceToHost, _cuda_copy_stream));
	}
	template <typename elem_type, MemBlock mem_block>
	inline void device_vector<elem_type, mem_block>::fill(elem_type _elem_value, uint64_t _first_index, uint64_t _last_index)
	{
		jarvis_memory<elem_type>::debug_error_if_index_out_of_range(_first_index, _last_index, "fill()");
		elem_type *ptr_tmp = jarvis_memory<elem_type>::ptr() + _first_index;
		uint64_t size_tmp = _last_index - _first_index + 1;
		void *args_list[] = {&ptr_tmp, &size_tmp, &_elem_value};
		cudaLaunchKernel((void *)_jarvis_fill_cuda_kernel<elem_type>, jarvis_cuda_kernel_size(size_tmp), args_list, 0, jarvis_default_cuda_stream);
	};
	template <typename elem_type, MemBlock mem_block>
	inline void device_vector<elem_type, mem_block>::set_value(elem_type _elem_value, uint64_t _index)
	{
		elem_type *ptr_tmp = jarvis_memory<elem_type>::ptr() + _index;
		void *args_list[] = {&ptr_tmp, &_elem_value};
		cudaLaunchKernel((void *)_jarvis_set_value_cuda_kernel<elem_type>, 1, 1, args_list, 0, jarvis_default_cuda_stream);
	}
	template <typename elem_type, MemBlock mem_block>
	inline elem_type device_vector<elem_type, mem_block>::get_value(uint64_t _index)
	{
		elem_type tmp;
		cudaMemcpy(&tmp, jarvis_memory<elem_type>::ptr() + _index, sizeof(elem_type), cudaMemcpyDeviceToHost);
		cudaStreamSynchronize(jarvis_default_cuda_stream);
		return tmp;
	}
	//
	template <typename elem_type>
	inline void vector<elem_type, MemType::pinned, MemBlock::single>::alloc(uint64_t _elem_num)
	{
		jarvis_memory<elem_type>::size() = _elem_num;
		jarvis_memory<elem_type>::debug_error_if_size_is_zero("jarvis::vector::alloc(uint64_t _elem_num)");
		jarvis_handle_cuda_error(cudaHostAlloc((void **)&jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size() * sizeof(elem_type), cudaHostAllocDefault));
		std::memset(jarvis_memory<elem_type>::ptr(), 0, jarvis_memory<elem_type>::size() * sizeof(elem_type));
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::pinned, MemBlock::single>::clear()
	{
		jarvis_memory<elem_type>::debug_error_if_ptr_is_null("jarvis::vector::clear()");
		jarvis_handle_cuda_error(cudaFreeHost(jarvis_memory<elem_type>::ptr()));
		jarvis_memory<elem_type>::ptr() = nullptr;
		jarvis_memory<elem_type>::size() = 0;
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::pinned, MemBlock::puppet>::alloc(vector<elem_type, MemType::pinned, MemBlock::multiple> &_multiple_vector, uint64_t _elem_num)
	{
		jarvis_memory<elem_type>::size() = _elem_num;
		_multiple_vector.push_back(this);
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::pinned, MemBlock::multiple>::alloc()
	{
		jarvis_memory<elem_type>::debug_error_if_size_is_zero("jarvis::vector::alloc()");
		jarvis_handle_cuda_error(cudaHostAlloc((void **)&(jarvis_memory<elem_type>::ptr()), jarvis_memory<elem_type>::size() * sizeof(elem_type), cudaHostAllocDefault));
		std::memset(jarvis_memory<elem_type>::ptr(), 0, jarvis_memory<elem_type>::size() * sizeof(elem_type));
		uint64_t offset = 0;
		for (auto it = jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.begin();
			 it != jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.end(); ++it)
		{
			(*it)->jarvis_memory<elem_type>::ptr() = jarvis_memory<elem_type>::ptr() + offset;
			offset += (*it)->jarvis_memory<elem_type>::size();
		}
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::pinned, MemBlock::multiple>::load(string _file_path)
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
	inline void vector<elem_type, MemType::pinned, MemBlock::multiple>::clear()
	{
		jarvis_memory<elem_type>::debug_error_if_ptr_is_null("jarvis::vector::clear()");
		jarvis_handle_cuda_error(cudaFreeHost(jarvis_memory<elem_type>::ptr()));
		jarvis_memory<elem_type>::ptr() = nullptr;
		jarvis_memory<elem_type>::size() = 0;
		for (auto it = jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.begin();
			 it != jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.end(); ++it)
		{
			(*it)->jarvis_memory<elem_type>::ptr() = nullptr;
			(*it)->jarvis_memory<elem_type>::size() = 0;
		}
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::device, MemBlock::single>::alloc(uint64_t _elem_num)
	{
		jarvis_memory<elem_type>::size() = _elem_num;
		jarvis_memory<elem_type>::debug_error_if_size_is_zero("jarvis::vector::alloc(uint64_t _elem_num)");
		jarvis_handle_cuda_error(cudaMalloc((void **)&jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size() * sizeof(elem_type)));
		jarvis_handle_cuda_error(cudaMemset(jarvis_memory<elem_type>::ptr(), 0, jarvis_memory<elem_type>::size() * sizeof(elem_type)));
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::device, MemBlock::single>::clear()
	{
		jarvis_memory<elem_type>::debug_error_if_ptr_is_null("jarvis::vector::clear()");
		jarvis_handle_cuda_error(cudaFree(jarvis_memory<elem_type>::ptr()));
		jarvis_memory<elem_type>::ptr() = nullptr;
		jarvis_memory<elem_type>::size() = 0;
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::device, MemBlock::puppet>::alloc(vector<elem_type, MemType::device, MemBlock::multiple> &_multiple_vector, uint64_t _elem_num)
	{
		jarvis_memory<elem_type>::size() = _elem_num;
		_multiple_vector.push_back(this);
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::device, MemBlock::multiple>::alloc()
	{
		jarvis_memory<elem_type>::debug_error_if_size_is_zero("jarvis::vector::alloc()");
		jarvis_handle_cuda_error(cudaMalloc((void **)&jarvis_memory<elem_type>::ptr(), jarvis_memory<elem_type>::size() * sizeof(elem_type)));
		jarvis_handle_cuda_error(cudaMemset(jarvis_memory<elem_type>::ptr(), 0, jarvis_memory<elem_type>::size() * sizeof(elem_type)));
		uint64_t offset = 0;
		for (auto it = jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.begin();
			 it != jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.end(); ++it)
		{
			(*it)->jarvis_memory<elem_type>::ptr() = jarvis_memory<elem_type>::ptr() + offset;
			offset += (*it)->jarvis_memory<elem_type>::size();
		}
	}
	template <typename elem_type>
	inline void vector<elem_type, MemType::device, MemBlock::multiple>::clear()
	{
		jarvis_memory<elem_type>::debug_error_if_ptr_is_null("jarvis::vector::clear()");
		jarvis_handle_cuda_error(cudaFree(jarvis_memory<elem_type>::ptr()));
		jarvis_memory<elem_type>::ptr() = nullptr;
		jarvis_memory<elem_type>::size() = 0;
		for (auto it = jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.begin();
			 it != jarvis_memory_block<elem_type, MemBlock::multiple>::puppet_memory_list.end(); ++it)
		{
			(*it)->jarvis_memory<elem_type>::ptr() = nullptr;
			(*it)->jarvis_memory<elem_type>::size() = 0;
		}
	}
}
#endif