#ifndef _JARVIS_HOST_MATH_MEAT
#define _JARVIS_HOST_MATH_MEAT
//**********************************Developer******************************************
// 2020.04.10 BY CAIWEI CALVIN CAI
//*************************************************************************************
#include "jarvis_1_host_math_bones.hpp"
namespace jarvis
{
	template <typename numeric_type>
	numeric_type min(numeric_type _num_1, numeric_type _num_2)
	{
		if (_num_1 > _num_2)
			return _num_2;
		else
			return _num_1;
	}
	template <typename numeric_type>
	numeric_type min(numeric_type _num_1, numeric_type _num_2, numeric_type _num_3)
	{
		return min<numeric_type>(min<numeric_type>(_num_1, _num_2), _num_3);
	}
	template <typename numeric_type>
	numeric_type max(numeric_type _num_1, numeric_type _num_2)
	{
		if (_num_1 > _num_2)
			return _num_1;
		else
			return _num_2;
	};
	template <typename numeric_type>
	numeric_type max(numeric_type _num_1, numeric_type _num_2, numeric_type _num_3)
	{
		return max<numeric_type>(max<numeric_type>(_num_1, _num_2), _num_3);
	};
	template <typename numeric_type>
	numeric_type min(numeric_type *_mem_ptr, uint64_t _size)
	{
		numeric_type _min_value;
		uint64_t i;
		_min_value = _mem_ptr[0];
		for (i = 0; i < _size; ++i)
		{
			if (_min_value > _mem_ptr[i])
			{
				_min_value = _mem_ptr[i];
			}
		}
		return _min_value;
	}
	template <typename numeric_type>
	numeric_type max(numeric_type *_mem_ptr, uint64_t _size)
	{
		numeric_type _max_value;
		uint64_t i;
		_max_value = _mem_ptr[0];
		for (i = 0; i < _size; ++i)
		{
			if (_max_value < _mem_ptr[i])
			{
				_max_value = _mem_ptr[i];
			}
		}
		return _max_value;
	}
	template <typename numeric_type>
	void min(numeric_type *_mem_ptr, uint64_t _size, numeric_type &_min_val, uint64_t &_min_val_index)
	{
		uint64_t i;
		_min_val = _mem_ptr[0];
		for (i = 0; i < _size; ++i)
		{
			if (_min_val > _mem_ptr[i])
			{
				_min_val = _mem_ptr[i];
				_min_val_index = i;
			}
		}
	}
	template <typename numeric_type>
	void max(numeric_type *_mem_ptr, uint64_t _size, numeric_type &_max_val, uint64_t &_max_val_index)
	{
		_max_val = _mem_ptr[0];
		uint64_t i;
		for (i = 0; i < _size; ++i)
		{
			if (_max_val < _mem_ptr[i])
			{
				_max_val = _mem_ptr[i];
				_max_val_index = i;
			}
		}
	}
	template <typename numeric_type>
	numeric_type abs_min(numeric_type *_mem_ptr, uint64_t _size)
	{
		uint64_t i;
		numeric_type temp = std::abs(_mem_ptr[0]);
		for (i = 1; i < _size; ++i)
		{
			if (temp > _mem_ptr[i])
				temp = std::abs(_mem_ptr[i]);
		}
		return temp;
	}
	template <typename numeric_type>
	numeric_type abs_max(numeric_type *_mem_ptr, uint64_t _size)
	{
		uint64_t i;
		numeric_type temp = std::abs(_mem_ptr[0]);
		for (i = 1; i < _size; ++i)
		{
			if (temp < _mem_ptr[i])
				temp = std::abs(_mem_ptr[i]);
		}
		return temp;
	}
	template <typename numeric_type>
	numeric_type abs_mean(numeric_type *_mem_ptr, uint64_t _size)
	{
		uint64_t i;
		numeric_type _sum = 0;
		for (i = 0; i < _size; ++i)
		{
			_sum += _mem_ptr[i] * _mem_ptr[i];
		}
		return sqrt(_sum / _size);
	}
	template <typename numeric_type>
	numeric_type mean(numeric_type *_mem_ptr, uint64_t _size)
	{
		uint64_t i;
		numeric_type _sum = 0;
		for (i = 0; i < _size; ++i)
		{
			_sum += _mem_ptr[i];
		}
		return _sum / _size;
	}
	template <typename numeric_type>
	numeric_type abs_sum(numeric_type *_mem_ptr, uint64_t _size)
	{
		numeric_type abs_sum = 0;
		uint64_t i;
		for (i = 0; i < _size; ++i)
		{
			abs_sum += abs(_mem_ptr[i]);
		}
		return abs_sum;
	}
	template <typename numeric_type>
	void limit(numeric_type &_num, double _low_boundary, double _up_boundary)
	{
		if (_num < _low_boundary)
			_num = _low_boundary;
		else if (_num > _up_boundary)
			_num = _up_boundary;
	}

	template <typename numeric_type>
	void swap(numeric_type &_num_1, numeric_type &_num_2)
	{
		numeric_type temp;
		temp = _num_1;
		_num_1 = _num_2;
		_num_2 = temp;
	}
	template <typename numeric_type>
	bool next_permutation(numeric_type *_mem_ptr, uint64_t _size)
	{
		uint32_t last = _size - 1;
		uint32_t i, j, k;
		static uint32_t count = 0;
		if (count == 0)
		{
			count++;
			return true;
		}
		i = last;
		while (i > 0 && _mem_ptr[i] < _mem_ptr[i - 1])
			i--;
		if (i == 0)
			return false;
		k = i;
		for (j = last; j >= i; j--)
			if (_mem_ptr[j] > _mem_ptr[i - 1] && _mem_ptr[j] < _mem_ptr[k])
				k = j;
		swap<numeric_type>(_mem_ptr[k], _mem_ptr[i - 1]);
		for (j = last, k = i; j > k; j--, k++)
			swap<numeric_type>(_mem_ptr[j], _mem_ptr[k]);
		count++;
		return true;
	}

	template <typename numeric_type>
	void show_permutation(numeric_type *_mem_ptr, uint64_t _size)
	{
		uint32_t i;
		for (i = 0; i < _size; ++i)
			std::cout << _mem_ptr[i];
		std::cout << std::endl;
	}
	//*********************************************************************************
}
#endif