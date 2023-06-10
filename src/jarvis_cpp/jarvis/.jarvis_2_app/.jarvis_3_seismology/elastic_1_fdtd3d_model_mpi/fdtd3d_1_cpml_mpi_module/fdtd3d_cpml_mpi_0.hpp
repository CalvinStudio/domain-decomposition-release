#pragma once
#ifndef _FDTD3D_CPML_MPI_BONES_HPP
#define _FDTD3D_CPML_MPI_BONES_HPP
#include ".fdtd3d_cpml_mpi_header_in.h"
namespace jarvis
{
    class mpicuCPML : public mpicuPMLBase
    {
    public:
        host_damp alpha_x, alpha_x_half;
        host_damp alpha_y, alpha_y_half;
        host_damp alpha_z, alpha_z_half;
        device_damp a_x, b_x, a_x_half, b_x_half;
        device_damp a_y, b_y, a_y_half, b_y_half;
        device_damp a_z, b_z, a_z_half, b_z_half;
        Field<float, MemType::paged, MemBlock::multiple> damp_host_memory;
        Field<float, MemType::paged_device, MemBlock::multiple> damp_device_memory;
        //
        
        //
        void cu_alloc_damp_coeff();    // 1rd
        void calculate_damp_coeff();   // 2th
        void cu_copy_damp_coeff_h2d(); // 3th

        void cuda_update_and_mpi_exchange_edge_vel();
        void cuda_update_and_mpi_exchange_edge_stress();
        //
        void clear();
    };
}
#endif