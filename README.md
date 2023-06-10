![](jarvis.png)

## Usage
`cd domain-decomposition-release\app\.main\fdtd3d_blend_elastic_model_mpi`

`sh run_rtm_cpml_mpi.sh`

## Background
3D finite-difference time-domain numerical simulation and reconstruction based on the domain decomposition technique are essential parts of high-performance computation for inverse-time migration and full-waveform inversion. However, the low GPU utilization in computing for small-sized models and the tremendous memory consumption for large-sized models may result in low computational efficiency and high memory costs. This paper proposes a contiguous memory management (CMM) method and a variable-order wave-field reconstruction (VOWR) method. The CMM allocates the memory of many small-sized arrays for MPI communications on a large-sized contiguous memory block, which aims to reduce the number of MPI communications between subdomains and improve the communication bandwidth, thus reducing the MPI time overhead and improving the GPU utilization. Meanwhile, the VOWR can flexibly set the number of layers of boundary wavefield used for source wavefield reconstruction according to the host memory capacity and accuracy requirements. It can store at least one layer of boundary wavefield, thus significantly alleviating the memory consumption of host memory. Numerical experiments show that the GPU utilization in computing for the model with a size of $120^3$ can be improved from 25\% to 90\% using the CMM method, and the VOWR method can reduce memory consumption by about 86\% while maintaining good accuracy in wave-field reconstruction.

## Introduce
This library is an easy-to-use seismic application code based on GPU and MPI programming and C++ template metaprogramming, which contains basic data types, related tools and applications.