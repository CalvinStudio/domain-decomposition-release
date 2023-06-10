#ifndef JARVIS_HOST_FIELD_MEAT
#define JARVIS_HOST_FIELD_MEAT
//**********************************Developer*************************************
// 2021.04.9 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_1_host_field_bones.hpp"
namespace jarvis
{
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	elem_type &host_single_field<elem_type, mem_type, mem_block>::operator()(uint64_t _i_rows, uint32_t _i_cols, uint32_t _i_slices)
	{
		return *(host::ptr() + _i_rows + _i_cols * n_rows + _i_slices * n_elem_slice);
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	elem_type host_single_field<elem_type, mem_type, mem_block>::operator()(uint64_t _i_rows, uint32_t _i_cols, uint32_t _i_slices) const
	{
		return *(host::ptr() + _i_rows + _i_cols * n_rows + _i_slices * n_elem_slice);
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	void host_single_field<elem_type, mem_type, mem_block>::save(string _file_path, SaveFormat _format)
	{
		switch (_format)
		{
		case SaveFormat::binary_raw:
			save_as_binary_raw(_file_path);
			break;
		case SaveFormat::binary_fld:
			save_as_binary_fld(_file_path);
			break;
		case SaveFormat::ascii_xyz:
			save_as_ascii_xyz(_file_path);
			break;
		case SaveFormat::ascii_grd:
			save_as_ascii_grd(_file_path);
			break;
		case SaveFormat::ascii_txt:
			save_as_ascii_txt(_file_path);
			break;
		}
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	void host_single_field<elem_type, mem_type, mem_block>::load(string _file_path)
	{
		string file_end_str = _file_path.substr(_file_path.size() - 4, _file_path.size());
		string end_str = ".fld";
		if (file_end_str != end_str)
		{
			std::cout << "load:\033[41;37m[error]:\033[0m " + _file_path + " Can't read field-raw file %s ,exit" << std::endl;
			std::abort();
		}
		FILE *fp;
		fp = fopen(_file_path.c_str(), "rb");
		if (!fp)
		{
			std::cout << "load():\033[41;37m[error]:\033[0m File open error!" << std::endl;
			std::abort();
		}
		fread((Frame *)this, sizeof(Frame), 1, fp);
		host::size() = Frame::n_elem;
		host::alloc(host::size());
		fread(host::ptr(), host::size() * sizeof(elem_type), 1, fp);
		fclose(fp);
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	void host_single_field<elem_type, mem_type, mem_block>::read_raw(string _file_path, const Frame &_frame)
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
		host::alloc(host::size());
		fread(host::ptr(), host::size() * sizeof(elem_type), 1, fp);
		fclose(fp);
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	void host_single_field<elem_type, mem_type, mem_block>::save_as_binary_fld(string _file_path)
	{
		_file_path += ".fld";
		FILE *fp = fopen(_file_path.c_str(), "wb");
		fwrite((Frame *)this, sizeof(Frame), 1, fp);
		fwrite(host::ptr(), host::size() * sizeof(elem_type), 1, fp);
		fclose(fp);
		printf("\033[45;30m[save fld-format file]:\033[0m%s\n", _file_path.c_str());
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	void host_single_field<elem_type, mem_type, mem_block>::save_as_binary_raw(string _file_path)
	{
		_file_path += size_info() + ".raw";
		FILE *fp = fopen(_file_path.c_str(), "wb");
		fwrite(host::ptr(), host::size() * sizeof(elem_type), 1, fp);
		fclose(fp);
		printf("\033[45;30m[save raw-format file]:\033[0m%s\n", _file_path.c_str());
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	void host_single_field<elem_type, mem_type, mem_block>::save_as_ascii_xyz(string _file_path)
	{
		if (type != Frame::Type::grid)
		{
			std::cout << "save():\033[41;37m[error]:\033[0m Format Error! check type!" << std::endl;
			std::abort();
		}
		FILE *fp;
		_file_path += ".txt";
		if ((fp = fopen(_file_path.c_str(), "w")) == NULL)
		{
			printf("grid_output:\033[41;37m[error]:\033[0mCan't create output file %s ,exit\n", _file_path.c_str());
			std::abort();
		}
		for (uint64_t idx = 0; idx < n_elem; ++idx)
		{
			uint64_t _i_rows = idx % n_elem_slice % n_rows;
			uint32_t _i_cols = idx % n_elem_slice / n_rows;
			uint32_t _i_slices = idx / n_elem_slice;
			double x = l_rows + _i_rows * d_rows;
			double y = l_cols + _i_cols * d_cols;
			double z = l_slices + _i_slices * d_slices;
			fprintf(fp, "%7.3f,%7.3f,%7.3f,%8.6e\n", x, y, z, *(host::ptr() + idx));
		}
		fclose(fp);
		printf("save X-Y-Z-V File:%s\n", _file_path.c_str());
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	void host_single_field<elem_type, mem_type, mem_block>::save_as_ascii_grd(string _file_path)
	{
		_file_path += ".grd";
		std::cout << "save .grd File:" << _file_path << std::endl;
		ofstream fout(_file_path);
		if (type != Frame::Type::mesh)
		{
			std::cout << "save():\033[41;37m[error]:\033[0m Format Error! check type!" << std::endl;
			std::abort();
		}
		if (!fout.is_open())
		{
			printf("save_as_ascii_grd():\033[41;37m[error]:\033[0mCan't create GRD output file %s ,exit\n", _file_path.c_str());
			std::abort();
		}
		fout << "DSAA" << std::endl;
		fout << n_rows << " " << n_cols << std::endl;
		fout << l_rows << " " << r_rows << std::endl;
		fout << l_cols << " " << r_cols << std::endl;
		fout << host::min() << " " << host::max() << std::endl;
		for (uint32_t _i_cols = 0; _i_cols < n_cols; ++_i_cols)
		{
			for (uint64_t _i_rows = 0; _i_rows < n_rows; ++_i_rows)
			{
				fout << *(host::ptr() + _i_rows + _i_cols * n_rows) << " ";
			}
			fout << std::endl;
		}
		fout.close();
	}
	template <typename elem_type, MemType mem_type, MemBlock mem_block>
	void host_single_field<elem_type, mem_type, mem_block>::save_as_ascii_txt(string _file_path)
	{
		FILE *fp;
		_file_path += ".txt";
		fp = fopen(_file_path.c_str(), "wt");
		if (!fp)
		{
			printf("save_as_ascii_txt():\033[41;37m[error]:\033[0mcan't create txt output file %s ,exit\n", _file_path.c_str());
			std::abort();
		}
		if (type == Frame::Type::null)
		{
			printf("save_as_ascii_txt():\033[41;37m[error]:\033[0mframe type is null!");
			std::abort();
		}
		if (type == Frame::Type::vec || type == Frame::Type::line || type == Frame::Type::mat || type == Frame::Type::mesh)
		{
			for (uint32_t _i_cols = 0; _i_cols < n_cols; ++_i_cols)
			{
				for (uint64_t _i_rows = 0; _i_rows < n_rows; ++_i_rows)
				{
					fprintf(fp, "%10.6e\t", *(host::ptr() + _i_rows + _i_cols * n_rows));
				}
				fprintf(fp, "\n");
			}
			printf("save 2d txt file:%s\n", _file_path.c_str());
			fclose(fp);
		}
		else if (type == Frame::Type::cube || type == Frame::Type::grid)
		{
			for (uint32_t _i_slices = 0; _i_slices < n_slices; ++_i_slices)
				for (uint32_t _i_cols = 0; _i_cols < n_cols; ++_i_cols)
					for (uint64_t _i_rows = 0; _i_rows < n_rows; ++_i_rows)
						fprintf(fp, "%d,%d,%d,%10.6e\n", _i_rows, _i_cols, _i_slices, *(host::ptr() + _i_rows + _i_cols * n_rows + _i_slices * n_elem_slice));

			printf("save 3d txt file:%s\n", _file_path.c_str());
			fclose(fp);
		}
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged, MemBlock::single>::alloc(const Frame &_frame)
	{
		Frame::copy(_frame);
		host::alloc(_frame.n_elem);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged, MemBlock::single>::clear()
	{
		host::clear();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged, MemBlock::puppet>::alloc(Field<elem_type, MemType::paged, MemBlock::multiple> &_multiple_field, const Frame &_frame)
	{
		Frame::copy(_frame);
		host::size() = _frame.n_elem;
		_multiple_field.push_back(this);
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged, MemBlock::multiple>::alloc()
	{
		host::alloc();
	}
	template <typename elem_type>
	inline void Field<elem_type, MemType::paged, MemBlock::multiple>::clear()
	{
		host::clear();
	}
}
#endif