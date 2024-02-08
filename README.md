<div  align="center">    
<img src="jarvis.png" width = 50%>
</div>

## Background

3D finite-difference time-domain numerical simulation and reconstruction based on the domain decomposition technique are essential parts of high-performance computation for inverse-time migration and full-waveform inversion. However, the low GPU utilization in computing for small-sized models and the tremendous memory consumption for large-sized models may result in low computational efficiency and high memory costs. This paper proposes a contiguous memory management (CMM) method and a variable-order wave-field reconstruction (VOWR) method. The CMM allocates the memory of many small-sized arrays for MPI communications on a large-sized contiguous memory block, which aims to reduce the number of MPI communications between subdomains and improve the communication bandwidth, thus reducing the MPI time overhead and improving the GPU utilization. Meanwhile, the VOWR can flexibly set the number of layers of boundary wavefield used for source wavefield reconstruction according to the host memory capacity and accuracy requirements. It can store at least one layer of boundary wavefield, thus significantly alleviating the memory consumption of host memory. Numerical experiments show that the GPU utilization in computing for the model with a size of $120^3$ can be improved from 25\% to 90\% using the CMM method, and the VOWR method can reduce memory consumption by about 86\% while maintaining good accuracy in wave-field reconstruction.

## Introduce

This library is an easy-to-use seismic application code based on GPU and MPI programming and C++ template metaprogramming, which contains basic data types, related tools and applications.

## Software and hardware requirements

Major hardware and software list of a single cluster node in my study:

  |Hardware/Software | Model/Version                  |Number|
  | :-----            | :----                         | :---- |
  |CPU               | Intel Xeon Gold 6226R 2.9GHz   |2|
  |GPU               | NVIDIA Tesla V100S 32GB        |4|
  |Host memory       | DDR4 2933MHz 32GB             | 8|
  |Operating system  | CentOS 3.10.0-693.el7.x86\_64 | -|
  |MPI library       | NVIDIA HPC-X OpenMPI v2.8.1   | -|
  |Compilers         | GCC v7.2.0 and NVCC v10.2    | -|

Feasible configurations:

  |Hardware/Software | Model/Version                  |Number|
  | :-----           | :----                          | :----  |
  |CPU               | Intel Xeon Gold 6226R 2.9GHz   |$\geq$ 1|
  |GPU               | NVIDIA GPU (Memory $\geq$ 8GB) |$\geq$ 2|
  |Host memory       | $\geq$ 32GB                    |$\geq$ 1|
  |Operating system  | Linux                          | -|
  |MPI library       | NVIDIA HPC-X OpenMPI (Support CUDA-Aware MPI)| -|
  |Compilers         | GCC (Support C++17) and NVCC $\geq$ v7.5     | -|

## Usage

`cd domain-decomposition-release\app\.main\fdtd3d_blend_elastic_model_mpi`

`sh run_rtm_cpml_mpi.sh`

## Modify CUDA and MPI environment variables

`cd domain-decomposition-release/build/jarvis/`

`vim jarvis_0_env`

```bash
#!/bin/bash
#env
CUDA_PATH="/public/software/cuda/10.2"
NVCC_INCLUDE="-I/public/software/cuda/10.2/include"
MPI_INCLUDE="-I/public/software/mpi/hpcx/2.8.1/ompi/include"
#
NVCC_LIBRARY="-L/public/software/cuda/10.2/lib64 -lcudart"
MPI_LIBRARY="-L/public/software/mpi/hpcx/2.8.1/ompi/lib"
```

## Modify the macro parameters of the code

`cd domain-decomposition-release/src/jarvis_cpp/jarvis/.jarvis_2_app/.jarvis_3_seismology/elastic_1_fdtd3d_model_mpi/fdtd3d_0_base_module/`

`vim fdtd3d_base_1_interface.hpp`

```cpp
    class GeoConst
    {
    public:
        static const uint32_t phy_fdorder_half = 4;//Half-difference order of SPR
        static const uint32_t pml_fdorder_half = 2;//Half-difference order of SAR
        static const uint32_t rec_width = 1;//The number of layers of storage used for wavefield reconstruction
    };
```

## Modify parameters of seismic simulation 

`cd domain-decomposition-release/app/.main/fdtd3d_blend_elastic_model_mpi/fd3d_proj_small/cube100/`

`vim argfile.json`

```json
"model_path": "cube1001.gms",//Geological model contains P-wave velocity, S-wave velocity, and density.
"shot_path": "shot.txt",//Shots coordinates; The first number indicates the number of shot points.
"rece_path": "rece.txt",//Geophones coordinates
"output_path": "out/",//output path of results.
//
"is_output": false,
//
"pml_num": 10,//Absorption boundary layer number
"fm": 100,//Main frequency
"dt": 0.002,//Time sampling interval
"T": 4,//Time duration of seismic　waves
"t_snap": -1,
//*mpi parameter
//Grid Dividec
"x_divide_num": 2,//the number of divisions at the x direction
"y_divide_num": 2,//the number of divisions at the y direction
"z_divide_num": 1,//the number of divisions at the z direction
"slots_per_node": 4//Number of processes started on each compute node.
```

`vim shot.txt`

```txt
2
15 15 15
25 25 25
```

`vim rece.txt`

```txt
8
15 15 15
25 15 15
35 15 15
45 15 15
55 15 15
65 15 15
75 15 15
85 15 15
```

##  Geological model input

`cd domain-decomposition-release/src/jarvis_cpp/jarvis/.jarvis_2_app/.jarvis_3_seismology/elastic_1_fdtd3d_model_mpi/fdtd3d_0_base_module/`

`vim fdtd3d_base_3_geomodel_meat.hpp`

> The model given in the example is in gms format, which integrates the information of P-wave velocity, S-wave velocity, density, etc
```cpp
GmsReader gms(model_path);
gms.read_gms_by_order_to_ffield(phy_vp);
gms.read_gms_by_order_to_ffield(phy_vs);
gms.read_gms_by_order_to_ffield(phy_rho);
gms.clear();
```

>If you want to input a custom model, you can do so in the following way
```cpp   
phy_vp.read_raw("phy_vp.raw",Frame(676, 676, 210, 20, 20, 20, 0, 0, 0));
phy_vs.read_raw("phy_vs.raw",Frame(676, 676, 210, 20, 20, 20, 0, 0, 0));
phy_rho.read_raw("phy_rho.raw",Frame(676, 676, 210, 20, 20, 20, 0, 0, 0));
```
This code can import binary data of P-wave velocity model, S-wave velocity model and density model into variable 'phy_vp', 'phy_vs' and 'phy_rho', respectively.
Note that the index order is X, then Y, then Z.
'Frame(676, 676, 210, 20, 20, 20, 0, 0, 0)' indicates that the data body size is 676x676x210, the grid spacing is 20, and the origin coordinate is 0.

## Migration results output

`cd domain-decomposition-release/app/.main/fdtd3d_blend_elastic_model_mpi/fd3d_proj_small/cube100/out/`

`image_laplace_0_50_50_101_0.000000_0.000000_0.000000.raw`

The numbers after 'image_laplace' are the process number, Nx, Ny, Nz, origin_x, origin_y, origin_z, respectively.

These numbers make it easy for you to draw 3D images in software such as Voxler.
