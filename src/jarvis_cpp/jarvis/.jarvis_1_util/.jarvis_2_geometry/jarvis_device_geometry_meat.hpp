#ifndef JARVIS_DEVICE_GEOMETRY_MEAT
#define JARVIS_DEVICE_GEOMETRY_MEAT
#include "jarvis_device_geometry_bones.hpp"
namespace jarvis
{
    inline void Geometry::read_shot_line()
    {
        std::ifstream infile;
        infile.open(shot_path);
        if (infile.is_open())
        {
            infile >> shot.n_elem;
            shot.alloc(shot.n_elem);
            for (int i = 0; i < shot.n_elem; i++)
            {
                infile >> shot(i).x;
                infile >> shot(i).y;
                infile >> shot(i).z;
                shot(i).ind = i;
            }
            infile.close();
            shot.copy_h2d();
        }
        else
        {
            printf("read_shot_line(path):file open error!\n");
            std::abort();
        }
    }
    //
    inline void Geometry::read_recv_line()
    {
        std::ifstream infile;
        infile.open(recv_path);
        if (infile.is_open())
        {
            infile >> recv.n_elem;
            recv.alloc(recv.n_elem);
            for (int i = 0; i < recv.n_elem; i++)
            {
                infile >> recv(i).x;
                infile >> recv(i).y;
                infile >> recv(i).z;
                recv(i).ind = i;
            }
            infile.close();
            recv.copy_h2d();
        }
        else
        {
            printf("read_recv_line(path):file open error!\n");
            std::abort();
        }
    }

    inline void SubGeometry::initialize(string _shot_path, string _recv_path, const Frame &sub_grid)
    {
        shot_path = _shot_path;
        recv_path = _recv_path;
        read_shot_line();
        read_recv_line();
        get_sub_shot_and_recv(sub_grid);
    }
    //
    inline void SubGeometry::get_sub_shot_and_recv(const Frame &sub_grid)
    {
        int n_sub_shot = 0;
        for (int ishot = 0; ishot < shot.n_elem; ishot++)
        {
            if (shot(ishot).x >= sub_grid.l_rows && shot(ishot).x <= sub_grid.r_rows &&
                shot(ishot).y >= sub_grid.l_cols && shot(ishot).y <= sub_grid.r_cols &&
                shot(ishot).z >= sub_grid.l_slices && shot(ishot).z <= sub_grid.r_slices)
            {
                n_sub_shot++;
            }
        }
        if (n_sub_shot > 0)
        {
            is_have_shot = true;
            sub_shot.alloc(n_sub_shot);
            int sub_shot_idx = 0;
            for (int ishot = 0; ishot < shot.n_elem; ishot++)
            {
                if (shot(ishot).x >= sub_grid.l_rows && shot(ishot).x <= sub_grid.r_rows &&
                    shot(ishot).y >= sub_grid.l_cols && shot(ishot).y <= sub_grid.r_cols &&
                    shot(ishot).z >= sub_grid.l_slices && shot(ishot).z <= sub_grid.r_slices)
                {
                    sub_shot(sub_shot_idx).x = shot(ishot).x;
                    sub_shot(sub_shot_idx).y = shot(ishot).y;
                    sub_shot(sub_shot_idx).z = shot(ishot).z;
                    sub_shot(sub_shot_idx).ind = shot(ishot).ind;
                    sub_shot_idx++;
                }
            }
            sub_shot.copy_h2d();
        }

        int n_sub_recv = 0;
        for (int irecv = 0; irecv < recv.n_elem; irecv++)
        {
            if (recv(irecv).x >= sub_grid.l_rows && recv(irecv).x <= sub_grid.r_rows &&
                recv(irecv).y >= sub_grid.l_cols && recv(irecv).y <= sub_grid.r_cols &&
                recv(irecv).z >= sub_grid.l_slices && recv(irecv).z <= sub_grid.r_slices)
            {
                n_sub_recv++;
            }
        }
        if (n_sub_recv > 0)
        {
            is_have_recv = true;
            sub_recv.alloc(n_sub_recv);
            int sub_recv_idx = 0;
            for (int irecv = 0; irecv < recv.n_elem; irecv++)
            {
                if (recv(irecv).x >= sub_grid.l_rows && recv(irecv).x <= sub_grid.r_rows &&
                    recv(irecv).y >= sub_grid.l_cols && recv(irecv).y <= sub_grid.r_cols &&
                    recv(irecv).z >= sub_grid.l_slices && recv(irecv).z <= sub_grid.r_slices)
                {
                    sub_recv(sub_recv_idx).x = recv(irecv).x;
                    sub_recv(sub_recv_idx).y = recv(irecv).y;
                    sub_recv(sub_recv_idx).z = recv(irecv).z;
                    sub_recv(sub_recv_idx).ind = recv(irecv).ind;
                    sub_recv_idx++;
                }
            }
            sub_recv.copy_h2d();
        }
    }
}

#endif