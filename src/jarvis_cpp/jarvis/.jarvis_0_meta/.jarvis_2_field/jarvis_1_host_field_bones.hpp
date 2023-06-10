#ifndef JARVIS_HOST_FIELD_BONES
#define JARVIS_HOST_FIELD_BONES
//**********************************Developer*************************************
// 2021.04.9 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_0_frame_meat.hpp"
namespace jarvis
{
	enum class SaveFormat
	{
		binary_fld,
		binary_raw,
		ascii_xyz,
		ascii_grd,
		ascii_txt,
	};
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	class host_single_field : public Frame, public vector<elem_type, mem_type, mem_block>
	{
	private:
		using host = vector<elem_type, mem_type, mem_block>;

	public:
		elem_type &operator()(uint64_t _i_rows, uint32_t _i_cols = 0, uint32_t _i_slices = 0);
		elem_type operator()(uint64_t _i_rows, uint32_t _i_cols = 0, uint32_t _i_slices = 0) const;
		void load(string _file_path);
		void read_raw(string _file_path, const Frame &_frame);
		void save(string _file_path, SaveFormat _format);

	protected:
		void save_as_binary_fld(string _file_path);
		void save_as_binary_raw(string _file_path);
		void save_as_ascii_xyz(string _file_path);
		void save_as_ascii_grd(string _file_path);
		void save_as_ascii_txt(string _file_path);
	};
	template <typename elem_type>
	class host_single_field<elem_type, MemType::device, MemBlock::single> : public Frame
	{
	};
	template <typename elem_type>
	class host_single_field<elem_type, MemType::device, MemBlock::puppet> : public Frame
	{
	};
	template <typename elem_type>
	class host_single_field<elem_type, MemType::device, MemBlock::multiple>
	{
	};
	template <typename elem_type, MemType mem_type = MemType::paged, MemBlock mem_block = MemBlock::single>
	class Field
	{
	};
	template <typename elem_type>
	class Field<elem_type, MemType::paged, MemBlock::single> : public host_single_field<elem_type, MemType::paged, MemBlock::single>
	{
	private:
		using host = vector<elem_type, MemType::paged, MemBlock::single>;

	public:
		void alloc(const Frame &_frame);
		void clear();
	};
	template <typename elem_type>
	class Field<elem_type, MemType::paged, MemBlock::puppet> : public host_single_field<elem_type, MemType::paged, MemBlock::puppet>
	{
	private:
		using host = vector<elem_type, MemType::paged, MemBlock::puppet>;

	public:
		void alloc(Field<elem_type, MemType::paged, MemBlock::multiple> &_multiple_field, const Frame &_frame);
	};
	template <typename elem_type>
	class Field<elem_type, MemType::paged, MemBlock::multiple> : public vector<elem_type, MemType::paged, MemBlock::multiple>
	{
	private:
		using host = vector<elem_type, MemType::paged, MemBlock::multiple>;

	public:
		void alloc();
		void clear();
	};
}
#endif