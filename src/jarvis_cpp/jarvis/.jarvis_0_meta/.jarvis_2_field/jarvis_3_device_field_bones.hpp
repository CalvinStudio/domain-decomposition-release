#ifndef JARVIS_DEVICE_FIELD_BONES
#define JARVIS_DEVICE_FIELD_BONES
//**********************************Developer*************************************
// 2021.05.14 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_2_device_stream_meat.hpp"
namespace jarvis
{
	template <typename elem_type, MemType host_mem_type, MemBlock mem_block>
	class device_single_field : public host_single_field<elem_type, host_mem_type, mem_block>, public vector<elem_type, MemType::device, mem_block>
	{
	private:
		using host = vector<elem_type, host_mem_type, mem_block>;
		using device = vector<elem_type, MemType::device, mem_block>;

	public:
		void extract_data_from(device_single_field<elem_type, MemType::device, MemBlock::single> &_bigger_cufield, jarvis_cuda_Stream *_jarvis_cuda_stream_p);
		void inject_data_into(device_single_field<elem_type, MemType::device, MemBlock::single> &_bigger_cufield, jarvis_cuda_Stream *_jarvis_cuda_stream_p);
		void set_value(elem_type _elem_value, uint64_t _i_rows, uint32_t _i_cols = 0, uint32_t _i_slices = 0);
		elem_type get_value(uint64_t _i_rows, uint32_t _i_cols = 0, uint32_t _i_slices = 0);
	};
	template <typename elem_type>
	class Field<elem_type, MemType::pinned, MemBlock::single> : public host_single_field<elem_type, MemType::pinned, MemBlock::single>
	{
	private:
		using host = vector<elem_type, MemType::pinned, MemBlock::single>;

	public:
		void alloc(const Frame &_frame);
	};
	template <typename elem_type>
	class Field<elem_type, MemType::pinned, MemBlock::puppet> : public host_single_field<elem_type, MemType::pinned, MemBlock::puppet>
	{
	private:
		using host = vector<elem_type, MemType::pinned, MemBlock::puppet>;

	public:
		void alloc(Field<elem_type, MemType::pinned, MemBlock::multiple> &_multiple_field, const Frame &_frame);
	};
	template <typename elem_type>
	class Field<elem_type, MemType::pinned, MemBlock::multiple> : public host_single_field<elem_type, MemType::pinned, MemBlock::multiple>
	{
	private:
		using host = vector<elem_type, MemType::pinned, MemBlock::multiple>;

	public:
		void alloc();
	};
	template <typename elem_type>
	class Field<elem_type, MemType::device, MemBlock::single> : public device_single_field<elem_type, MemType::device, MemBlock::single>
	{
	private:
		using device = vector<elem_type, MemType::device, MemBlock::single>;

	public:
		void alloc(const Frame &_frame);
		void cu_save(string _file_path, SaveFormat _format);
		void clear();
	};
	template <typename elem_type>
	class Field<elem_type, MemType::device, MemBlock::puppet> : public device_single_field<elem_type, MemType::device, MemBlock::puppet>
	{
	private:
		using device = vector<elem_type, MemType::device, MemBlock::puppet>;

	public:
		void alloc(Field<elem_type, MemType::device, MemBlock::multiple> &_multiple_field, const Frame &_frame);
	};
	template <typename elem_type>
	class Field<elem_type, MemType::device, MemBlock::multiple> : public device_single_field<elem_type, MemType::device, MemBlock::multiple>
	{
	private:
		using device = vector<elem_type, MemType::device, MemBlock::multiple>;

	public:
		void alloc();
		void clear();
	};
	template <typename elem_type>
	class Field<elem_type, MemType::paged_device, MemBlock::single> : public device_single_field<elem_type, MemType::paged, MemBlock::single>
	{
	public:
		using host = vector<elem_type, MemType::paged, MemBlock::single>;
		using device = vector<elem_type, MemType::device, MemBlock::single>;

	public:
		void alloc(const Frame &_frame);
		void copy_h2d();
		void copy_d2h();
		void read_raw(string _file_path, const Frame &_frame);
		void cu_save(string _file_path, SaveFormat _format);
		void clear();
	};
	template <typename elem_type>
	class Field<elem_type, MemType::paged_device, MemBlock::puppet> : public device_single_field<elem_type, MemType::paged, MemBlock::puppet>
	{
	public:
		using host = vector<elem_type, MemType::paged, MemBlock::puppet>;
		using device = vector<elem_type, MemType::device, MemBlock::puppet>;

	public:
		void alloc(Field<elem_type, MemType::paged_device, MemBlock::multiple> &_multiple_field, const Frame &_frame);
	};
	template <typename elem_type>
	class Field<elem_type, MemType::paged_device, MemBlock::multiple> : public vector<elem_type, MemType::paged, MemBlock::multiple>, public vector<elem_type, MemType::device, MemBlock::multiple>
	{
	public:
		using host = vector<elem_type, MemType::paged, MemBlock::multiple>;
		using device = vector<elem_type, MemType::device, MemBlock::multiple>;

	public:
		void alloc();
		void copy_h2d();
		void copy_d2h();
		void clear();
	};
	template <typename elem_type>
	class Field<elem_type, MemType::pinned_device, MemBlock::single> : public device_single_field<elem_type, MemType::pinned, MemBlock::single>
	{
	public:
		using host = vector<elem_type, MemType::pinned, MemBlock::single>;
		using device = vector<elem_type, MemType::device, MemBlock::single>;

	public:
		void alloc(const Frame &_frame);
		void copy_h2d(const cudaStream_t &_h2d_stream);
		void copy_d2h(const cudaStream_t &_d2h_stream);
		void cu_save(string _file_path, SaveFormat _format, const cudaStream_t &_d2h_stream);
		void clear();
	};
	template <typename elem_type>
	class Field<elem_type, MemType::pinned_device, MemBlock::puppet> : public device_single_field<elem_type, MemType::pinned, MemBlock::puppet>
	{
	public:
		using host = vector<elem_type, MemType::pinned, MemBlock::puppet>;
		using device = vector<elem_type, MemType::device, MemBlock::puppet>;

	public:
		void alloc(Field<elem_type, MemType::pinned_device, MemBlock::multiple> &_multiple_field, const Frame &_frame);
	};
	template <typename elem_type>
	class Field<elem_type, MemType::pinned_device, MemBlock::multiple> : public vector<elem_type, MemType::pinned, MemBlock::multiple>, public vector<elem_type, MemType::device, MemBlock::multiple>
	{
	public:
		using host = vector<elem_type, MemType::pinned, MemBlock::multiple>;
		using device = vector<elem_type, MemType::device, MemBlock::multiple>;

	public:
		void alloc();
		void copy_h2d(const cudaStream_t &_h2d_stream);
		void copy_d2h(const cudaStream_t &_d2h_stream);
		void clear();
	};
	//
#define set_cufield_2d_idx(_mesh, _idx, _i_rows, _i_cols) \
	int _idx = threadIdx.x + blockDim.x * blockIdx.x;     \
	int n_rows = (mesh).n_rows;                           \
	int n_cols = (mesh).n_cols;                           \
	int _i_rows = _idx % n_rows;                          \
	int _i_cols = _idx / n_rows

#define set_cufield_3d_idx(_grid, _idx, _i_rows, _i_cols, _i_slices) \
	int _idx = threadIdx.x + blockDim.x * blockIdx.x;                \
	int n_rows = (_grid).n_rows;                                     \
	int n_elem_slice = (_grid).n_elem_slice;                         \
	int _i_rows = _idx % n_elem_slice % n_rows;                      \
	int _i_cols = _idx % n_elem_slice / n_rows;                      \
	int _i_slices = _idx / n_elem_slice

#define set_cufield_grid_3d_idx(ext_grid, grid)                                         \
	int idx = threadIdx.x + blockDim.x * blockIdx.x;                                    \
	int n_rows = (grid).n_rows;                                                         \
	int n_elem_slice = (grid).n_elem_slice;                                             \
	int i = idx % n_elem_slice % n_rows;                                                \
	int j = idx % n_elem_slice / n_rows;                                                \
	int k = idx / n_elem_slice;                                                         \
	int ext_n_rows = (ext_grid).n_rows;                                                 \
	int ext_n_elem_slice = (ext_grid).n_elem_slice;                                     \
	int ext_i = i + get_gap_num((ext_grid).l_rows, (grid).l_rows, (grid).d_rows);       \
	int ext_j = j + get_gap_num((ext_grid).l_cols, (grid).l_cols, (grid).d_cols);       \
	int ext_k = k + get_gap_num((ext_grid).l_slices, (grid).l_slices, (grid).d_slices); \
	int ext_idx = ext_i + ext_j * ext_n_rows + ext_k * ext_n_elem_slice
}
#endif