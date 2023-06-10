#ifndef JARVIS_DEVICE_FIELD_MEAT
#define JARVIS_DEVICE_FIELD_MEAT
//**********************************Developer*************************************
// 2021.05.14 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_3_device_field_bones.hpp"
namespace jarvis
{
	inline __host__ __device__ int get_gap_num(float l, float r, float d)
	{
		int n_l, n_r;
		float v = 1.0f / d;
		if (l > 0)
			n_l = l * v + 0.125f;
		else
			n_l = l * v - 0.125f;
		if (r > 0)
			n_r = r * v + 0.125f;
		else
			n_r = r * v - 0.125f;
		return n_r - n_l;
	}
	template <typename elem_type>
	inline __global__ void _jarvis_extract_cuda_kernel(Frame _bigger_grid, Frame _grid, elem_type *_bigger_ptr, elem_type *_ptr)
	{
		int idx = threadIdx.x + blockDim.x * blockIdx.x;
		if (idx < _grid.n_elem)
		{
			int n_rows = _grid.n_rows;
			int n_elem_slice = _grid.n_elem_slice;
			int a = idx % n_elem_slice;
			int i = a % n_rows;
			int j = a / n_rows;
			int k = idx / n_elem_slice;
			int bigger_i = i + get_gap_num(_bigger_grid.l_rows, _grid.l_rows, _bigger_grid.d_rows);
			int bigger_j = j + get_gap_num(_bigger_grid.l_cols, _grid.l_cols, _bigger_grid.d_cols);
			int bigger_k = k + get_gap_num(_bigger_grid.l_slices, _grid.l_slices, _bigger_grid.d_slices);
			int bigger_idx = bigger_i + bigger_j * _bigger_grid.n_rows + bigger_k * _bigger_grid.n_elem_slice;
			_ptr[idx] = _bigger_ptr[bigger_idx];
		}
	}
	template <typename elem_type>
	inline __global__ void _jarvis_inject_cuda_kernel(Frame _bigger_grid, Frame _grid, elem_type *_bigger_ptr, elem_type *_ptr)
	{
		int idx = threadIdx.x + blockDim.x * blockIdx.x;
		if (idx < _grid.n_elem)
		{
			int n_rows = _grid.n_rows;
			int n_elem_slice = _grid.n_elem_slice;
			int a = idx % n_elem_slice;
			int i = a % n_rows;
			int j = a / n_rows;
			int k = idx / n_elem_slice;
			int bigger_i = i + get_gap_num(_bigger_grid.l_rows, _grid.l_rows, _bigger_grid.d_rows);
			int bigger_j = j + get_gap_num(_bigger_grid.l_cols, _grid.l_cols, _bigger_grid.d_cols);
			int bigger_k = k + get_gap_num(_bigger_grid.l_slices, _grid.l_slices, _bigger_grid.d_slices);
			int bigger_idx = bigger_i + bigger_j * _bigger_grid.n_rows + bigger_k * _bigger_grid.n_elem_slice;
			_bigger_ptr[bigger_idx] = _ptr[idx];
		}
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	inline void device_single_field<elem_type, mem_type, mem_block>::extract_data_from(device_single_field<elem_type, MemType::device, MemBlock::single> &_bigger_cufield, jarvis_cuda_Stream *_jarvis_cuda_stream_p)
	{
		void *args_list[] = {(Frame *)&_bigger_cufield, (Frame *)this, &_bigger_cufield.ptr(), &device::ptr()};
		cudaLaunchKernel((void *)_jarvis_extract_cuda_kernel<elem_type>, jarvis_cuda_kernel_size(Frame::n_elem), args_list, 0, _jarvis_cuda_stream_p->cal_stream());
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	inline void device_single_field<elem_type, mem_type, mem_block>::inject_data_into(device_single_field<elem_type, MemType::device, MemBlock::single> &_bigger_cufield, jarvis_cuda_Stream *_jarvis_cuda_stream_p)
	{
		void *args_list[] = {(Frame *)&_bigger_cufield, (Frame *)this, &_bigger_cufield.ptr(), &device::ptr()};
		cudaLaunchKernel((void *)_jarvis_inject_cuda_kernel<elem_type>, jarvis_cuda_kernel_size(Frame::n_elem), args_list, 0, _jarvis_cuda_stream_p->cal_stream());
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	inline void device_single_field<elem_type, mem_type, mem_block>::set_value(elem_type _elem_value, uint64_t _i_rows, uint32_t _i_cols, uint32_t _i_slices)
	{
		device::set_value(_elem_value, _i_rows + _i_cols * Frame::n_rows + _i_slices * Frame::n_elem_slice);
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	inline elem_type device_single_field<elem_type, mem_type, mem_block>::get_value(uint64_t _i_rows, uint32_t _i_cols, uint32_t _i_slices)
	{
		return device::get_value(_i_rows + _i_cols * Frame::n_rows + _i_slices * Frame::n_elem_slice);
	}
	//
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned, MemBlock::single>::alloc(const Frame &_frame)
	{
		_frame.debug_error_if_frame_is_empty("alloc()");
		Frame::copy(_frame);
		host::alloc(Frame::n_elem);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned, MemBlock::puppet>::alloc(Field<elem_type, MemType::pinned, MemBlock::multiple> &_multiple_field, const Frame &_frame)
	{
		Frame::copy(_frame);
		host::size() = _frame.n_elem;
		_multiple_field.push_back(this);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned, MemBlock::multiple>::alloc()
	{
		host::alloc();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::device, MemBlock::single>::alloc(const Frame &_frame)
	{
		_frame.debug_error_if_frame_is_empty("alloc()");
		Frame::copy(_frame);
		device::alloc(Frame::n_elem);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::device, MemBlock::single>::cu_save(string _file_path, SaveFormat _format)
	{
		Frame::debug_error_if_frame_is_empty("cu_save()");
		Field<elem_type, MemType::pinned, MemBlock::single> _host_tmp;
		_host_tmp.alloc((Frame)(*this));
		device::copy_into_host(_host_tmp, jarvis_default_cuda_stream);
		cudaDeviceSynchronize();
		_host_tmp.save(_file_path, _format);
		_host_tmp.clear();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::device, MemBlock::single>::clear()
	{
		device::clear();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::device, MemBlock::puppet>::alloc(Field<elem_type, MemType::device, MemBlock::multiple> &_multiple_field, const Frame &_frame)
	{
		Frame::copy(_frame);
		device::size() = Frame::n_elem;
		_multiple_field.push_back((device *)this);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::device, MemBlock::multiple>::alloc()
	{
		device::alloc();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::device, MemBlock::multiple>::clear()
	{
		device::clear();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged_device, MemBlock::single>::alloc(const Frame &_frame)
	{
		_frame.debug_error_if_frame_is_empty("alloc()");
		Frame::copy(_frame);
		host::alloc(Frame::n_elem);
		device::alloc(Frame::n_elem);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged_device, MemBlock::single>::copy_h2d()
	{
		jarvis_handle_cuda_error(cudaMemcpy(device::ptr(), host::ptr(), Frame::n_elem * sizeof(elem_type), cudaMemcpyHostToDevice));
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged_device, MemBlock::single>::copy_d2h()
	{
		jarvis_handle_cuda_error(cudaMemcpy(host::ptr(), device::ptr(), Frame::n_elem * sizeof(elem_type), cudaMemcpyDeviceToHost));
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged_device, MemBlock::single>::cu_save(string _file_path, SaveFormat _format)
	{
		Frame::debug_error_if_frame_is_empty("cu_save()");
		device::copy_into_host((host)(*this));
		cudaDeviceSynchronize();
		host::save(_file_path, _format);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged_device, MemBlock::single>::read_raw(string _file_path, const Frame &_frame)
	{
		FILE *fp;
		fp = fopen(_file_path.c_str(), "rb");
		if (!fp)
		{
			std::cout << "load():\033[41;37m[error]:\033[0m File open error!" << std::endl;
			std::abort();
		}
		(*this).copy(_frame);
		host::size() = Frame::n_elem;
		device::size() = Frame::n_elem;
		host::alloc(host::size());
		device::alloc(device::size());
		fread(host::ptr(), host::size() * sizeof(elem_type), 1, fp);
		fclose(fp);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged_device, MemBlock::single>::clear()
	{
		host::clear();
		device::clear();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged_device, MemBlock::puppet>::alloc(Field<elem_type, MemType::paged_device, MemBlock::multiple> &_multiple_field, const Frame &_frame)
	{
		Frame::copy(_frame);
		host::size() = Frame::n_elem;
		device::size() = Frame::n_elem;
		_multiple_field.vector<elem_type, MemType::paged, MemBlock::multiple>::push_back((host *)this);
		_multiple_field.vector<elem_type, MemType::device, MemBlock::multiple>::push_back((device *)this);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged_device, MemBlock::multiple>::alloc()
	{
		host::alloc();
		device::alloc();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged_device, MemBlock::multiple>::copy_h2d()
	{
		jarvis_handle_cuda_error(cudaMemcpy(device::ptr(), host::ptr(), host::size() * sizeof(elem_type), cudaMemcpyHostToDevice));
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged_device, MemBlock::multiple>::copy_d2h()
	{
		jarvis_handle_cuda_error(cudaMemcpy(host::ptr(), device::ptr(), device::size() * sizeof(elem_type), cudaMemcpyDeviceToHost));
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged_device, MemBlock::multiple>::clear()
	{
		host::clear();
		device::clear();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned_device, MemBlock::single>::alloc(const Frame &_frame)
	{
		_frame.debug_error_if_frame_is_empty("alloc()");
		Frame::copy(_frame);
		host::alloc(Frame::n_elem);
		device::alloc(Frame::n_elem);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned_device, MemBlock::single>::copy_h2d(const cudaStream_t &_h2d_stream)
	{
		jarvis_handle_cuda_error(cudaMemcpyAsync(device::ptr(), host::ptr(), Frame::n_elem * sizeof(elem_type), cudaMemcpyHostToDevice, _h2d_stream));
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned_device, MemBlock::single>::copy_d2h(const cudaStream_t &_d2h_stream)
	{
		jarvis_handle_cuda_error(cudaMemcpyAsync(host::ptr(), device::ptr(), Frame::n_elem * sizeof(elem_type), cudaMemcpyDeviceToHost, _d2h_stream));
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned_device, MemBlock::single>::cu_save(string _file_path, SaveFormat _format, const cudaStream_t &_d2h_stream)
	{
		Frame::debug_error_if_frame_is_empty("cu_save()");
		device::copy_into_host((host)(*this), _d2h_stream);
		cudaDeviceSynchronize();
		host::save(_file_path, _format);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned_device, MemBlock::single>::clear()
	{
		host::clear();
		device::clear();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned_device, MemBlock::puppet>::alloc(Field<elem_type, MemType::pinned_device, MemBlock::multiple> &_multiple_field, const Frame &_frame)
	{
		Frame::copy(_frame);
		host::size() = Frame::n_elem;
		device::size() = Frame::n_elem;
		_multiple_field.vector<elem_type, MemType::pinned, MemBlock::multiple>::push_back((host *)this);
		_multiple_field.vector<elem_type, MemType::device, MemBlock::multiple>::push_back((device *)this);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned_device, MemBlock::multiple>::copy_h2d(const cudaStream_t &_h2d_stream)
	{
		jarvis_handle_cuda_error(cudaMemcpyAsync(device::ptr(), host::ptr(), host::size() * sizeof(elem_type), cudaMemcpyHostToDevice, _h2d_stream));
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned_device, MemBlock::multiple>::copy_d2h(const cudaStream_t &_d2h_stream)
	{
		jarvis_handle_cuda_error(cudaMemcpyAsync(host::ptr(), device::ptr(), device::size() * sizeof(elem_type), cudaMemcpyDeviceToHost, _d2h_stream));
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned_device, MemBlock::multiple>::alloc()
	{
		host::alloc();
		device::alloc();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::pinned_device, MemBlock::multiple>::clear()
	{
		host::clear();
		device::clear();
	}
}
#endif