#ifndef JARVIS_HOST_FRAME_BONES
#define JARVIS_HOST_FRAME_BONES
//**********************************Developer*************************************
// 2021.04.9 BY CAIWEI CALVIN CAI
//********************************************************************************
#include ".jarvis_0_field_in.h"
namespace jarvis
{
	struct Frame
	{
		uint64_t n_elem = 0, n_rows = 0, n_elem_slice = 0;
		uint32_t n_cols = 0, n_slices = 0;
		float d_rows = 0, d_cols = 0, d_slices = 0;
		float l_rows = 0, l_cols = 0, l_slices = 0;
		float r_rows = 0, r_cols = 0, r_slices = 0;
		enum class Type
		{
			null = 0,
			vec = 1,
			mat = 2,
			cube = 3,
			line = 10,
			mesh = 20,
			grid = 30
		};
		Type type = Type::null;
		Frame() {}
		Frame(uint64_t _n_rows, uint32_t _n_cols = 1, uint32_t _n_slices = 1);
		Frame(uint64_t _n_rows, uint32_t _n_cols, uint32_t _n_slices, float _d_rows, float _d_cols, float _d_slices, float _l_rows, float _l_cols, float _l_slices);
		void set_n(uint64_t _n_rows, uint32_t _n_cols = 1, uint32_t _n_slices = 1);
		void set_nd(uint64_t _n_rows, float _d_rows);
		void set_nd(uint64_t _n_rows, uint32_t _n_cols, float _d_rows, float _d_cols);
		void set_nd(uint64_t _n_rows, uint32_t _n_cols, uint32_t _n_slices, float _d_rows, float _d_cols, float _d_slices);
		void set_ndl(uint64_t _n_rows, uint32_t _n_cols, uint32_t _n_slices, float _d_rows, float _d_cols, float _d_slices, float _l_rows, float _l_cols, float _l_slices);
		void copy(const Frame &_frame);
		void re_setndl();
		bool is_mat();
		bool is_mesh();
		bool is_cube();
		bool is_grid();
		bool operator==(const Frame &_frame);
		void print_info(string s = "");
		Frame sparse(uint32_t a);
		Type get_model_type();
		string size_info() const;
		bool is_empty() const;
		void debug_error_if_frame_is_empty(string _func_name) const;

	private:
		void set_frame_type();
	};

#define set_frame_n(frame)                        \
	uint32_t n_rows = (frame).n_rows;             \
	uint32_t n_elem_slice = (frame).n_elem_slice; \
	uint32_t n_cols = (frame).n_cols;             \
	uint32_t n_slices = (frame).n_slices

#define set_frame_nd(frame)                       \
	uint32_t n_rows = (frame).n_rows;             \
	uint32_t n_elem_slice = (frame).n_elem_slice; \
	uint32_t n_cols = (frame).n_cols;             \
	uint32_t n_slices = (frame).n_slices;         \
	float d_rows = (frame).d_rows;                \
	float d_cols = (frame).d_cols;                \
	float d_slices = (frame).d_slices

#define set_frame_ndl(frame)                      \
	uint32_t n_rows = (frame).n_rows;             \
	uint32_t n_elem_slice = (frame).n_elem_slice; \
	uint32_t n_cols = (frame).n_cols;             \
	uint32_t n_slices = (frame).n_slices;         \
	float d_rows = (frame).d_rows;                \
	float d_cols = (frame).d_cols;                \
	float d_slices = (frame).d_slices;            \
	float l_rows = (frame).l_rows;                \
	float l_cols = (frame).l_cols;                \
	float l_slices = (frame).l_slices

#define set_frame_ndlr(frame)                     \
	uint32_t n_rows = (frame).n_rows;             \
	uint32_t n_elem_slice = (frame).n_elem_slice; \
	uint32_t n_cols = (frame).n_cols;             \
	uint32_t n_slices = (frame).n_slices;         \
	float d_rows = (frame).d_rows;                \
	float d_cols = (frame).d_cols;                \
	float d_slices = (frame).d_slices;            \
	float l_rows = (frame).l_rows;                \
	float l_cols = (frame).l_cols;                \
	float l_slices = (frame).l_slices;            \
	float r_rows = (frame).r_rows;                \
	float r_cols = (frame).r_cols;                \
	float r_slices = (frame).r_slices

#define set_frame_nd_pre(pre, frame)                    \
	uint32_t pre##_n_rows = (frame).n_rows;             \
	uint32_t pre##_n_elem_slice = (frame).n_elem_slice; \
	uint32_t pre##_n_cols = (frame).n_cols;             \
	uint32_t pre##_n_slices = (frame).n_slices

#define frame_for(frame, _i_rows, _i_cols, _i_slices)                       \
	for (uint32_t _i_slices = 0; _i_slices < (frame).n_slices; ++_i_slices) \
		for (uint32_t _i_cols = 0; _i_cols < (frame).n_cols; ++_i_cols)     \
			for (uint64_t _i_rows = 0; _i_rows < (frame).n_rows; ++_i_rows)
}
#endif