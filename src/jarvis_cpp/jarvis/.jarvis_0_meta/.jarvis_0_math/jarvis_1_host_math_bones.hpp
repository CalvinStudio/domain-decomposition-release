#ifndef _JARVIS_HOST_MATH_BONES
#define _JARVIS_HOST_MATH_BONES
//**********************************Developer******************************************
// 2020.04.10 BY CAIWEI CALVIN CAI
//*************************************************************************************
#include "jarvis_0_host_math_const.hpp"
namespace jarvis
{
	template <typename numeric_type>
	numeric_type min(numeric_type _num_1, numeric_type _num_2);
	template <typename numeric_type>
	numeric_type max(numeric_type _num_1, numeric_type _num_2);
	template <typename numeric_type>
	numeric_type min(numeric_type _num_1, numeric_type _num_2, numeric_type _num_3);
	template <typename numeric_type>
	numeric_type max(numeric_type _num_1, numeric_type _num_2, numeric_type _num_3);
	template <typename numeric_type>
	numeric_type min(numeric_type *_mem_ptr, uint64_t _size);
	template <typename numeric_type>
	numeric_type max(numeric_type *_mem_ptr, uint64_t _size);
	template <typename numeric_type>
	numeric_type mean(numeric_type *_mem_ptr, uint64_t _size);
	template <typename numeric_type>
	numeric_type abs_min(numeric_type *_mem_ptr, uint64_t _size);
	template <typename numeric_type>
	numeric_type abs_max(numeric_type *_mem_ptr, uint64_t _size);
	template <typename numeric_type>
	numeric_type abs_mean(numeric_type *_mem_ptr, uint64_t _size);
	template <typename numeric_type>
	numeric_type abs_sum(numeric_type *_mem_ptr, uint64_t _size);
	template <typename numeric_type>
	void min(numeric_type *_mem_ptr, uint64_t _size, numeric_type &_min_val, uint64_t &_min_val_index);
	template <typename numeric_type>
	void max(numeric_type *_mem_ptr, uint64_t _size, numeric_type &_max_val, uint64_t &_min_val_index);
	template <typename numeric_type>
	void limit(numeric_type &_num, double _low_boundary, double _up_boundary);
	template <typename numeric_type>
	void swap(numeric_type &_num_1, numeric_type &_num_2);
	template <typename numeric_type>
	bool next_permutation(numeric_type *_mem_ptr, uint64_t _size);
	template <typename numeric_type>
	void show_permutation(numeric_type *_mem_ptr, uint64_t _size);
	//*********************************************************************************
	template <typename numeric_type>
	class Component2
	{
	public:
		numeric_type x;
		numeric_type y;
		enum class c
		{
			x = 0,
			y,
		};
	};
	//
	struct Point2D : Component2<double>
	{
		uint32_t ind;
		void print();
		double get_distance_to(Point2D a);
	};
	//
	struct uvec2 : public Component2<double>
	{
		bool is_unit()
		{
			if (x * x + y * y > 0.999999)
				return true;
			else
				return false;
		}
	};
	//
	template <typename numeric_type>
	class Component3
	{
	public:
		numeric_type x;
		numeric_type y;
		numeric_type z;
		enum class c
		{
			x = 0,
			y,
			z
		};
	};
	//
	typedef Component3<float> fcomp3;
	//
	struct Point3D : public Component3<float>
	{
		uint32_t ind;
		void print();
		void set_pos(Point3D a);
		double get_distance_to(Point3D a);
		Point3D get_pos();
	};
	//
	struct uvec3 : public Component3<float>
	{
		bool is_unit()
		{
			if (x * x + y * y + z * z > 0.999)
				return true;
			else
				return false;
		}
	};

#define jarvis_error_if_ptr_is_null(_ptr, _const_name)               \
	if (_ptr == nullptr)                                             \
	{                                                                \
		printf("%s is \033[41;37mnullptr\033[0m!!!\n", _const_name); \
		std::abort();                                                \
	}
//* time
#define tic(mark)                                 \
	clock_t clock_begin_##mark, clock_end_##mark; \
	clock_begin_##mark = clock();

#define toc(mark, string)                                                                          \
	clock_end_##mark = clock();                                                                    \
	double elapsed_time_##mark = (double)(clock_end_##mark - clock_begin_##mark) / CLOCKS_PER_SEC; \
	std::cout << string << elapsed_time_##mark << " seconds" << std::endl;
}
#endif