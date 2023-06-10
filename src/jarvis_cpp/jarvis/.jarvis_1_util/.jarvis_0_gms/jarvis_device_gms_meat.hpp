#ifndef JARVIS_DEVICE_GMS_MEAT
#define JARVIS_DEVICE_GMS_MEAT
#include "jarvis_device_gms_bones.hpp"
namespace jarvis
{
	inline GmsReader::GmsReader(std::string path)
	{
		gms_model_path = path;
		read_gms_data();
	}
	inline void GmsReader::read_gms_data()
	{
		FILE *fp = fopen(gms_model_path.c_str(), "rb");
		if (fp)
		{
			fread(&file_header, 240, 1, fp);
			grid.set_ndl(file_header.n_rows,
						 file_header.n_cols,
						 file_header.n_slices,
						 file_header.d_rows,
						 file_header.d_cols,
						 file_header.d_slices,
						 file_header.l_rows,
						 file_header.l_cols,
						 file_header.l_slices);
			para_num = file_header.para_num;
			data_header.alloc(para_num);
			data.alloc(para_num);
			for (int i = 0; i < para_num; i++)
			{
				data[i].alloc(grid.n_elem);
				fread(data_header.ptr(), 64, 1, fp);
				fread(data(i).ptr(), grid.n_elem * sizeof(float), 1, fp);
			}
			fclose(fp);
		}
		else
		{
			printf("read_gms_data():\033[41;37m[error]:\033[0mFile open error!");
			std::abort();
		}
	}
	inline void GmsReader::read_gms_by_order_to_ffield(Field<float> &ffield_obj)
	{
		if (para_order < para_num)
		{
			ffield_obj.alloc(grid);
			for (int k = 0; k < grid.n_slices; k++)
				for (int j = 0; j < grid.n_cols; j++)
					for (int i = 0; i < grid.n_rows; i++)
						ffield_obj(i, j, k) = data[para_order][i * grid.n_slices + j * grid.n_rows * grid.n_slices + k];
			para_order++;
		}
		else
		{
			printf("ERROR:There are no more parameters to read!");
			std::abort();
		}
	}

	inline void GmsReader::read_gms_by_order_to_fcufield(Field<float, MemType::paged_device> &fcufld_obj)
	{
		if (para_order < para_num)
		{
			fcufld_obj.alloc(grid);
			for (int j = 0; j < grid.n_cols; j++)
				for (int i = 0; i < grid.n_rows; i++)
					for (int k = 0; k < grid.n_slices; k++)
						fcufld_obj(i, j, k) = data[para_order][i * grid.n_slices + j * grid.n_rows * grid.n_slices + k];
			fcufld_obj.copy_h2d();
			cudaDeviceSynchronize();
			para_order++;
		}
		else
		{
			printf("ERROR:There are no more parameters to read!");
			std::abort();
		}
	}
	inline void GmsReader::print_info()
	{
		std::cout << "Grid Information:" << std::endl;
		std::cout << "para_order:" << para_order << ";"
				  << "para_num:" << para_num << ";" << std::endl;
		grid.print_info("GMS_MODEL:");
	}
	inline void GmsReader::clear()
	{
		for (int i = 0; i < para_num; i++)
		{
			data[i].clear();
		}
		data.clear();
	}
}
#endif
