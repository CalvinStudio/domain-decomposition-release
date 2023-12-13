![](jarvis.png)

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
  | | |

Feasible configurations:

  |Hardware/Software | Model/Version                  |Number|
  | :-----           | :----                          | :----  |
  |CPU               | Intel Xeon Gold 6226R 2.9GHz   |$\geq$ 1|
  |GPU               | NVIDIA GPU (Memory $\geq$ 8GB) |$\geq$ 1|
  |Host memory       | $\geq$ 32GB                    |$\geq$ 1|
  |Operating system  | Linux                          | -|
  |MPI library       | NVIDIA HPC-X OpenMPI (Support CUDA-Aware MPI)| -|
  |Compilers         | GCC (Support C++17) and NVCC $\geq$ v7.5     | -|
  | | |

## Usage

`cd domain-decomposition-release\app\.main\fdtd3d_blend_elastic_model_mpi`

`sh run_rtm_cpml_mpi.sh`
