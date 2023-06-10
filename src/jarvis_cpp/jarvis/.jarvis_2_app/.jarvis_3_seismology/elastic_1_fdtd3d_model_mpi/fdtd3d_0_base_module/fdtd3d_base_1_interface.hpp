#pragma once
#ifndef _FDTD3D_BASE_INTERFACE_HPP
#define _FDTD3D_BASE_INTERFACE_HPP
#include ".fdtd3d_base_header_in.h"
#ifndef __CUDACC__
#include "../../../../../.external_header/.public_json_header.h"
#endif
namespace jarvis
{
    class GeoConst
    {
    public:
        static const uint32_t phy_fdorder_half = 4;
        static const uint32_t pml_fdorder_half = 2;
        static const uint32_t rec_width = 1;
    };
    using geo_const = GeoConst;
    //
    enum class SimulateType
    {
        pure_forward,
        rtm_forward,
        rtm_backward,
        rtm_reverse,
        fwi
    };
    struct ArgList
    {
        string proj_path;
        string model_path;
        string shot_path;
        string recv_path;
        //*
        int mpi_rank;
        int mpi_size;
        int x_divide_num;
        int y_divide_num;
        int z_divide_num;
        int slots_per_node;
        //*
        int pml_num;
        //*
        float fm;
        float T;
        float dt;
        bool is_output;
        string output_path;
        int t_snap;
        //*
        void mpi_read_arg_file(string _proj_path, int _mpi_rank, int _mpi_size)
        {
#ifndef __CUDACC__
            std::ifstream fin;
            fin.open(_proj_path + "argfile.json");
            Json::Reader reader;
            Json::Value value;
            reader.parse(fin, value, false);
            if (fin.is_open())
            {
                model_path = _proj_path + value["model_path"].asString();
                shot_path = _proj_path + value["shot_path"].asString();
                recv_path = _proj_path + value["rece_path"].asString();
                //*
                x_divide_num = value["x_divide_num"].asInt();
                y_divide_num = value["y_divide_num"].asInt();
                z_divide_num = value["z_divide_num"].asInt();
                slots_per_node = value["slots_per_node"].asInt();
                //*
                pml_num = value["pml_num"].asInt();
                //*
                fm = value["fm"].asFloat();
                T = value["T"].asFloat();
                dt = value["dt"].asFloat();
                is_output = value["is_output"].asBool();
                output_path = _proj_path + value["output_path"].asString();
                t_snap = value["t_snap"].asInt();
                //
                proj_path = _proj_path;
                mpi_rank = _mpi_rank;
                mpi_size = _mpi_size;

                if (mpi_rank == 0)
                {
                    std::cout << "-------------------------" << std::endl;
                    std::cout << "proj_path:\t" << proj_path << std::endl;
                    std::cout << "model_path:\t" << model_path << std::endl;
                    std::cout << "shot_path:\t" << shot_path << std::endl;
                    std::cout << "recv_path:\t" << recv_path << std::endl;
                    std::cout << "x_divide_num:\t" << x_divide_num << std::endl;
                    std::cout << "y_divide_num:\t" << y_divide_num << std::endl;
                    std::cout << "z_divide_num:\t" << z_divide_num << std::endl;
                    std::cout << "mpi_size:\t" << mpi_size << std::endl;
                    std::cout << "slots_per_node:\t" << slots_per_node << std::endl;
                    std::cout << "pml_num:\t" << pml_num << std::endl;
                    std::cout << "fm:     \t" << fm << std::endl;
                    std::cout << "dt:     \t" << dt << std::endl;
                    std::cout << "T:      \t" << T << std::endl;
                    std::cout << "is_output:\t" << is_output << std::endl;
                    std::cout << "output_path:\t" << output_path << std::endl;
                    std::cout << "-------------------------" << std::endl;
                }
            }
            else
            {
                std::cout << "\033[41;37marg_file open error!\033[0m" << std::endl;
                std::abort();
            }
#else
            std::cout << "mpi_read_arg_file\033[41;37m[error]:please donot compile this function with nvcc!:\033[0m" << std::endl;
            std::abort();
#endif
        }
    };

    class elastic_fdtd3d_module_Base
    {
    public:
        static jarvis_global_mpicu_Stream_t *jarvis_mpi_cuda_stream_p;
        static ArgList *arg_list_p;
    };
}
#endif