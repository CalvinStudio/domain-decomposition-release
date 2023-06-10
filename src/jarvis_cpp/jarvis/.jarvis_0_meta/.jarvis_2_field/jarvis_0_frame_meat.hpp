#ifndef JARVIS_HOST_FRAME_MEAT
#define JARVIS_HOST_FRAME_MEAT
//**********************************Developer*************************************
// 2021.04.9 BY CAIWEI CALVIN CAI
//********************************************************************************
#include "jarvis_0_frame_bones.hpp"
namespace jarvis
{
	inline Frame::Frame(uint64_t _n_rows, uint32_t _n_cols, uint32_t _n_slices)
	{
		n_elem = _n_rows * _n_cols * _n_slices;
		n_rows = _n_rows;
		n_elem_slice = _n_rows * _n_cols;
		n_cols = _n_cols;
		n_slices = _n_slices;
		//
		d_rows = 0;
		d_cols = 0;
		d_slices = 0;
		l_rows = 0;
		l_cols = 0;
		l_slices = 0;
		r_rows = 0;
		r_cols = 0;
		r_slices = 0;
		set_frame_type();
	}
	inline Frame::Frame(uint64_t _n_rows, uint32_t _n_cols, uint32_t _n_slices, float _d_rows, float _d_cols, float _d_slices, float _l_rows, float _l_cols, float _l_slices)
	{
		n_rows = _n_rows;
		n_cols = _n_cols;
		n_slices = _n_slices;
		n_elem_slice = n_rows * n_cols;
		n_elem = n_rows * n_cols * n_slices;
		//
		d_rows = _d_rows;
		d_cols = _d_cols;
		d_slices = _d_slices;
		l_rows = _l_rows;
		l_cols = _l_cols;
		l_slices = _l_slices;
		if (n_rows > 0)
		{
			r_rows = l_rows + (n_rows - 1) * d_rows;
		}
		else
		{
			l_rows = r_rows = 0;
		}
		if (n_cols > 0)
		{
			r_cols = l_cols + (n_cols - 1) * d_cols;
		}
		else
		{
			l_cols = r_cols = 0;
		}
		if (n_slices > 0)
		{
			r_slices = l_slices + (n_slices - 1) * d_slices;
		}
		else
		{
			l_slices = r_slices = 0;
		}
		set_frame_type();
	}
	inline void Frame::set_n(uint64_t _n_rows, uint32_t _n_cols, uint32_t _n_slices)
	{
		n_rows = _n_rows;
		n_cols = _n_cols;
		n_slices = _n_slices;
		n_elem_slice = _n_rows * _n_cols;
		n_elem = _n_rows * _n_cols * _n_slices;
		//
		d_rows = 0;
		d_cols = 0;
		d_slices = 0;
		l_rows = 0;
		l_cols = 0;
		l_slices = 0;
		r_rows = 0;
		r_cols = 0;
		r_slices = 0;
		set_frame_type();
	}

	inline void Frame::set_nd(uint64_t _n_rows, float _d_rows)
	{
		n_rows = _n_rows;
		n_cols = 1;
		n_slices = 1;
		n_elem_slice = n_rows * n_cols;
		n_elem = n_rows * n_cols * n_slices;
		d_rows = _d_rows;
		d_cols = 0;
		d_slices = 0;
		l_rows = 0;
		l_cols = 0;
		l_slices = 0;
		r_rows = l_rows + (n_rows - 1) * d_rows;
		r_cols = 0;
		r_slices = 0;
		set_frame_type();
	}

	inline void Frame::set_nd(uint64_t _n_rows, uint32_t _n_cols, float _d_rows, float _d_cols)
	{
		n_rows = _n_rows;
		n_cols = _n_cols;
		n_slices = 1;
		n_elem_slice = n_rows * n_cols;
		n_elem = n_rows * n_cols * n_slices;
		//
		d_rows = _d_rows;
		d_cols = _d_cols;
		d_slices = 0;
		l_rows = 0;
		l_cols = 0;
		l_slices = 0;
		r_rows = l_rows + (n_rows - 1) * d_rows;
		r_cols = l_cols + (n_cols - 1) * d_cols;
		r_slices = 0;
		set_frame_type();
	}

	inline void Frame::set_nd(uint64_t _n_rows, uint32_t _n_cols, uint32_t _n_slices, float _d_rows, float _d_cols, float _d_slices)
	{
		n_rows = _n_rows;
		n_cols = _n_cols;
		n_slices = _n_slices;
		n_elem_slice = n_rows * n_cols;
		n_elem = n_rows * n_cols * n_slices;
		//
		d_rows = _d_rows;
		d_cols = _d_cols;
		d_slices = _d_slices;
		l_rows = 0;
		l_cols = 0;
		l_slices = 0;
		r_rows = l_rows + (n_rows - 1) * d_rows;
		r_cols = l_cols + (n_cols - 1) * d_cols;
		r_slices = l_slices + (n_slices - 1) * d_slices;
		set_frame_type();
	}

	inline void Frame::set_ndl(uint64_t _n_rows, uint32_t _n_cols, uint32_t _n_slices, float _d_rows, float _d_cols, float _d_slices, float _l_rows, float _l_cols, float _l_slices)
	{
		n_rows = _n_rows;
		n_cols = _n_cols;
		n_slices = _n_slices;
		n_elem_slice = n_rows * n_cols;
		n_elem = n_rows * n_cols * n_slices;
		//
		d_rows = _d_rows;
		d_cols = _d_cols;
		d_slices = _d_slices;
		l_rows = _l_rows;
		l_cols = _l_cols;
		l_slices = _l_slices;
		if (n_rows > 0)
		{
			r_rows = l_rows + (n_rows - 1) * d_rows;
		}
		else
		{
			l_rows = r_rows = 0;
		}
		if (n_cols > 0)
		{
			r_cols = l_cols + (n_cols - 1) * d_cols;
		}
		else
		{
			l_cols = r_cols = 0;
		}
		if (n_slices > 0)
		{
			r_slices = l_slices + (n_slices - 1) * d_slices;
		}
		else
		{
			l_slices = r_slices = 0;
		}
		set_frame_type();
	}

	inline void Frame::re_setndl()
	{
		r_rows = l_rows + (n_rows - 1) * d_rows;
		r_cols = l_cols + (n_cols - 1) * d_cols;
		r_slices = l_slices + (n_slices - 1) * d_slices;
		if (n_elem_slice <= 0)
			n_elem_slice = n_rows * n_cols;
		if (n_elem <= 0)
			n_elem = n_elem_slice * n_slices;
		set_frame_type();
	}

	inline void Frame::copy(const Frame &_frame)
	{
		n_rows = _frame.n_rows;
		n_elem_slice = _frame.n_elem_slice;
		n_elem = _frame.n_elem;
		n_cols = _frame.n_cols;
		n_slices = _frame.n_slices;
		d_rows = _frame.d_rows;
		d_cols = _frame.d_cols;
		d_slices = _frame.d_slices;
		l_rows = _frame.l_rows;
		l_cols = _frame.l_cols;
		l_slices = _frame.l_slices;
		r_rows = _frame.r_rows;
		r_cols = _frame.r_cols;
		r_slices = _frame.r_slices;
		type = _frame.type;
	}

	inline void Frame::set_frame_type()
	{
		if (d_rows > 0 || d_cols > 0 || d_slices > 0)
		{
			if ((n_rows == 1 && n_cols == 1) || (n_rows == 1 && n_slices == 1) || (n_cols == 1 && n_slices == 1))
				type = Type::line;
			else if ((n_rows > 1 && n_cols > 1 && n_slices == 1) || (n_rows > 1 && n_slices > 1 && n_cols == 1) || (n_cols > 1 && n_slices > 1 && n_rows == 1))
				type = Type::mesh;
			else if (n_rows > 1 && n_cols > 1 && n_slices > 1)
				type = Type::grid;
		}
		else
		{
			if ((n_rows == 1 && n_cols == 1) || (n_rows == 1 && n_slices == 1) || (n_cols == 1 && n_slices == 1))
				type = Type::vec;
			else if ((n_rows > 1 && n_cols > 1 && n_slices == 1) || (n_rows > 1 && n_slices > 1 && n_cols == 1) || (n_cols > 1 && n_slices > 1 && n_rows == 1))
				type = Type::mat;
			else if (n_rows > 1 && n_cols > 1 && n_slices > 1)
				type = Type::cube;
		}
	}

	inline Frame::Type Frame::get_model_type()
	{
		return type;
	}

	inline bool Frame::is_mat()
	{
		if ((n_rows > 1 && n_cols > 1 && n_slices == 1) ||
			(n_rows > 1 && n_slices > 1 && n_cols == 1) ||
			(n_cols > 1 && n_slices > 1 && n_rows == 1) &&
				(d_rows == 0 && d_cols == 0 && d_slices == 0))
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	inline bool Frame::is_mesh()
	{
		if ((n_rows > 1 && d_rows > 0 && n_cols > 1 && d_cols > 0 && n_slices == 1 && d_slices == 0) ||
			(n_rows > 1 && d_rows > 0 && n_slices > 1 && d_slices > 0 && n_cols == 1 && d_cols == 0) ||
			(n_cols > 1 && d_cols > 0 && n_slices > 1 && d_slices > 0 && n_rows == 1 && d_rows == 0))
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	inline bool Frame::is_cube()
	{
		if (n_rows > 1 && n_cols > 1 && n_slices > 1 && d_rows == 0 && d_cols == 0 && d_slices == 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	inline bool Frame::is_grid()
	{
		if (n_rows > 1 && n_cols > 1 && n_slices > 1 && d_rows > 0 && d_cols > 0 && d_slices > 0)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	inline bool Frame::operator==(const Frame &_frame)
	{
		if (n_rows == _frame.n_rows && n_cols == _frame.n_cols && n_slices == _frame.n_slices &&
			d_rows == _frame.d_rows && d_cols == _frame.d_cols && d_slices == _frame.d_slices &&
			l_rows == _frame.l_rows && l_cols == _frame.l_cols && l_slices == _frame.l_slices &&
			r_rows == _frame.r_rows && r_cols == _frame.r_cols && r_slices == _frame.r_slices)
		{
			return true;
		}
		else
		{
			return false;
		}
	}

	inline void Frame::print_info(string s)
	{
		std::cout << "\033[42;37m" << s << "\033[0m" << std::endl;
		std::cout << "grid_x_num:" << setw(5) << n_rows << "; interval:" << setw(6) << d_rows << "; range:"
				  << "[" << l_rows << ", " << r_rows << "]"
				  << "\n"
				  << "grid_y_num:" << setw(5) << n_cols << "; interval:" << setw(6) << d_cols << "; range:"
				  << "[" << l_cols << ", " << r_cols << "]"
				  << "\n"
				  << "grid_z_num:" << setw(5) << n_slices << "; interval:" << setw(6) << d_slices << "; range:"
				  << "[" << l_slices << ", " << r_slices << "]"
				  << "\n";
	}

	inline Frame Frame::sparse(uint32_t a)
	{
		Frame frame_t;
		frame_t.n_rows = (n_rows - 1) / a + 1;
		frame_t.n_cols = (n_cols - 1) / a + 1;
		frame_t.n_slices = (n_slices - 1) / a + 1;
		frame_t.n_elem_slice = frame_t.n_rows * frame_t.n_cols;
		frame_t.n_elem = frame_t.n_rows * frame_t.n_cols * frame_t.n_slices;
		frame_t.d_rows = d_rows * a;
		frame_t.d_cols = d_cols * a;
		frame_t.d_slices = d_slices * a;
		frame_t.l_rows = l_rows;
		frame_t.l_cols = l_cols;
		frame_t.l_slices = l_slices;
		frame_t.r_rows = frame_t.l_rows + (frame_t.n_rows - 1) * frame_t.d_rows;
		frame_t.r_cols = frame_t.l_cols + (frame_t.n_cols - 1) * frame_t.d_cols;
		frame_t.r_slices = frame_t.l_slices + (frame_t.n_slices - 1) * frame_t.d_slices;
		return frame_t;
	}
	inline string Frame::size_info() const
	{
		return "_" + to_string(n_rows) + "_" + to_string(n_cols) + "_" + to_string(n_slices) +
			   "_" + to_string(l_rows) + "_" + to_string(l_cols) + "_" + to_string(l_slices);
	}
	inline void Frame::debug_error_if_frame_is_empty(string _func_name) const
	{
		if (this->n_elem <= 0)
		{
			std::cout << _func_name + "\033[41;37m[error]\033[0m:frame is empty!" << std::endl;
			std::abort();
		}
	}
	inline bool Frame::is_empty() const
	{
		if (n_elem == 0)
			return true;
		else
			return false;
	}
}
#endif