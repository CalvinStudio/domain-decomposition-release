#ifndef JARVIS_DEVICE_GEOMETRY_BONES
#define JARVIS_DEVICE_GEOMETRY_BONES
#include ".jarvis_device_geometry_header_in.h"
namespace jarvis
{
    struct Geometry
    {
        string shot_path;
        string recv_path;
        Field<Point3D, MemType::paged_device> shot;
        Field<Point3D, MemType::paged_device> recv;
        //
        void read_shot_line();
        void read_recv_line();
    };

    struct SubGeometry : public Geometry
    {
    public:
        bool is_have_shot = false;
        bool is_have_recv = false;
        Field<Point3D, MemType::paged_device> sub_shot;
        Field<Point3D, MemType::paged_device> sub_recv;
        //
        void initialize(string _shot_path, string _recv_path, const Frame &_sub_grid);
        //
        void get_sub_shot_and_recv(const Frame &_sub_grid);
    };
    inline void cu_cal_diff_coeff(Field<float,MemType::paged_device> &_diff_coeff, int _order_num)
    {
        _diff_coeff.alloc(Frame(_order_num));
        for (int n = 0; n < _order_num; n++)
        {
            double a1 = 1, a2 = 1, a3 = 1;
            for (int i = 0; i < _order_num; i++)
            {
                if (i != n)
                    a1 = a1 * (2 * i + 1) * (2 * i + 1);
                if (i < n)
                    a2 = a2 * ((2 * n + 1) * (2 * n + 1) - (2 * i + 1) * (2 * i + 1));
                if (i > n)
                    a3 = a3 * ((2 * i + 1) * (2 * i + 1) - (2 * n + 1) * (2 * n + 1));
            }
            a1 = a1 * std::pow(-1, n);
            _diff_coeff(n) = a1 / ((2 * n + 1) * a2 * a3);
        }
        _diff_coeff.copy_h2d();
        ;
    }
    //
    inline void cu_cal_diff_coeff_list(Field<float,MemType::paged_device> &_diff_coeff, int _order_num_start, int _order_num_end)
    {
        _diff_coeff.alloc(Frame(_order_num_end, _order_num_end - _order_num_start + 1));
        for (int i_order = _order_num_start; i_order <= _order_num_end; i_order++)
            for (int n = 0; n < i_order; n++)
            {
                double a1 = 1, a2 = 1, a3 = 1;
                for (int i = 0; i < i_order; i++)
                {
                    if (i != n)
                        a1 = a1 * (2 * i + 1) * (2 * i + 1);
                    if (i < n)
                        a2 = a2 * ((2 * n + 1) * (2 * n + 1) - (2 * i + 1) * (2 * i + 1));
                    if (i > n)
                        a3 = a3 * ((2 * i + 1) * (2 * i + 1) - (2 * n + 1) * (2 * n + 1));
                }
                a1 = a1 * std::pow(-1, n);
                _diff_coeff(n, i_order - _order_num_start) = a1 / ((2 * n + 1) * a2 * a3);
            }
        _diff_coeff.copy_h2d();
    }
}
#endif