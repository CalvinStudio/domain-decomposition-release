#ifndef JARVIS_DEVICE_GMS_BONES
#define JARVIS_DEVICE_GMS_BONES
#include ".jarvis_device_gms_header_in.h"
namespace jarvis
{
	struct GmsFileHeader
	{
		short int format;
		char name[64];
		short int year;
		char mouth;
		char day;
		char hour;
		char minutes;
		char seconds;
		char type;
		short int para_num;
		float l_rows;
		float l_cols;
		float l_slices;
		uint32_t n_rows;
		uint32_t n_cols;
		uint32_t n_slices;
		float d_rows;
		float d_cols;
		float d_slices;
		char reserved[127];
	};
	struct GmsDataHeader
	{
		char dataName[64];
		char reserved[64];
	};
	class GmsReader
	{
	public:
		Frame grid;
		GmsFileHeader file_header;
		Field<GmsDataHeader> data_header;
		Field<Field<float, MemType::paged>> data;
		GmsReader(std::string path);
		void read_gms_by_order_to_ffield(Field<float, MemType::paged> &field_obj);
		void read_gms_by_order_to_fcufield(Field<float, MemType::paged_device> &field_obj);
		void print_info();
		void clear();

	private:
		int para_order = 0;
		int para_num;
		std::string gms_model_path;
		void read_gms_data();
	};
}
#endif
