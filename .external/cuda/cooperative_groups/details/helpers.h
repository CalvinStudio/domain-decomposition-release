 /* Copyright 1993-2016 NVIDIA Corporation.  All rights reserved.
  *
  * NOTICE TO LICENSEE:
  *
  * The source code and/or documentation ("Licensed Deliverables") are
  * subject to NVIDIA intellectual property rights under U.S. and
  * international Copyright laws.
  *
  * The Licensed Deliverables contained herein are PROPRIETARY and
  * CONFIDENTIAL to NVIDIA and are being provided under the terms and
  * conditions of a form of NVIDIA software license agreement by and
  * between NVIDIA and Licensee ("License Agreement") or electronically
  * accepted by Licensee.  Notwithstanding any terms or conditions to
  * the contrary in the License Agreement, reproduction or disclosure
  * of the Licensed Deliverables to any third party without the express
  * written consent of NVIDIA is prohibited.
  *
  * NOTWITHSTANDING ANY TERMS OR CONDITIONS TO THE CONTRARY IN THE
  * LICENSE AGREEMENT, NVIDIA MAKES NO REPRESENTATION ABOUT THE
  * SUITABILITY OF THESE LICENSED DELIVERABLES FOR ANY PURPOSE.  THEY ARE
  * PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF ANY KIND.
  * NVIDIA DISCLAIMS ALL WARRANTIES WITH REGARD TO THESE LICENSED
  * DELIVERABLES, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY,
  * NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
  * NOTWITHSTANDING ANY TERMS OR CONDITIONS TO THE CONTRARY IN THE
  * LICENSE AGREEMENT, IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY
  * SPECIAL, INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OR ANY
  * DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
  * WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
  * ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
  * OF THESE LICENSED DELIVERABLES.
  *
  * U.S. Government End Users.  These Licensed Deliverables are a
  * "commercial item" as that term is defined at 48 C.F.R. 2.101 (OCT
  * 1995), consisting of "commercial computer software" and "commercial
  * computer software documentation" as such terms are used in 48
  * C.F.R. 12.212 (SEPT 1995) and are provided to the U.S. Government
  * only as a commercial end item.  Consistent with 48 C.F.R.12.212 and
  * 48 C.F.R. 227.7202-1 through 227.7202-4 (JUNE 1995), all
  * U.S. Government End Users acquire the Licensed Deliverables with
  * only those rights set forth herein.
  *
  * Any use of the Licensed Deliverables in individual and commercial
  * software must include, in the user documentation and internal
  * comments to the code, the above Disclaimer and U.S. Government End
  * Users Notice.
  */

#ifndef _COOPERATIVE_GROUPS_HELPERS_H_
# define _COOPERATIVE_GROUPS_HELPERS_H_

#include "info.h"
#include "sync.h"

_CG_BEGIN_NAMESPACE

namespace details {
#ifdef _CG_CPP11_FEATURES
    template <typename Ty> struct _is_float_or_half          : public _CG_STL_NAMESPACE::is_floating_point<Ty> {};
# ifdef _CG_HAS_FP16_COLLECTIVE
    template <>            struct _is_float_or_half<__half>  : public _CG_STL_NAMESPACE::true_type {};
    template <>            struct _is_float_or_half<__half2> : public _CG_STL_NAMESPACE::true_type {};
# endif
    template <typename Ty>
    using  is_float_or_half = _is_float_or_half<typename _CG_STL_NAMESPACE::remove_cv<Ty>::type>;
#endif

    _CG_STATIC_QUALIFIER unsigned long long vec3_to_linear(dim3 index, dim3 nIndex) {
        return ((unsigned long long)index.z * nIndex.y * nIndex.x) +
               ((unsigned long long)index.y * nIndex.x) +
                (unsigned long long)index.x;
    }

    namespace tile {
        template <unsigned int TileCount, unsigned int TileMask, unsigned int LaneMask, unsigned int ShiftCount>
        struct __tile_helpers{
            _CG_STATIC_CONST_DECL unsigned int tileCount = TileCount;
            _CG_STATIC_CONST_DECL unsigned int tileMask = TileMask;
            _CG_STATIC_CONST_DECL unsigned int laneMask = LaneMask;
            _CG_STATIC_CONST_DECL unsigned int shiftCount = ShiftCount;
        };

        template <unsigned int> struct tile_helpers;
        template <> struct tile_helpers<32> : public __tile_helpers<1,  0xFFFFFFFF, 0x1F, 5> {};
        template <> struct tile_helpers<16> : public __tile_helpers<2,  0x0000FFFF, 0x0F, 4> {};
        template <> struct tile_helpers<8>  : public __tile_helpers<4,  0x000000FF, 0x07, 3> {};
        template <> struct tile_helpers<4>  : public __tile_helpers<8,  0x0000000F, 0x03, 2> {};
        template <> struct tile_helpers<2>  : public __tile_helpers<16, 0x00000003, 0x01, 1> {};
        template <> struct tile_helpers<1>  : public __tile_helpers<32, 0x00000001, 0x00, 0> {};
    };

    namespace multi_grid {
        struct multi_grid_functions;
    };

    namespace grid {
        _CG_STATIC_QUALIFIER void sync(unsigned int *bar) {
            unsigned int expected = gridDim.x * gridDim.y * gridDim.z;

            details::sync_grids(expected, bar);
        }

        _CG_STATIC_QUALIFIER unsigned long long size()
        {
            return ((unsigned long long)blockDim.z * gridDim.z) *
                   ((unsigned long long)blockDim.y * gridDim.y) *
                   ((unsigned long long)blockDim.x * gridDim.x);
        }

        _CG_STATIC_QUALIFIER unsigned long long thread_rank()
        {
            unsigned long long blkIdx = vec3_to_linear(blockIdx, gridDim);
            return (blkIdx * (blockDim.x * blockDim.y * blockDim.z)) + vec3_to_linear(threadIdx, blockDim);
        }

        _CG_STATIC_QUALIFIER dim3 grid_dim()
        {
            return (dim3(gridDim.x, gridDim.y, gridDim.z));
        }
    };


#if defined(_CG_HAS_MULTI_GRID_GROUP)

    namespace multi_grid {
        _CG_STATIC_QUALIFIER unsigned long long get_intrinsic_handle()
        {
            return (cudaCGGetIntrinsicHandle(cudaCGScopeMultiGrid));
        }

        _CG_STATIC_QUALIFIER void sync(const unsigned long long handle)
        {
            cudaError_t err = cudaCGSynchronize(handle, 0);
        }

        _CG_STATIC_QUALIFIER unsigned int size(const unsigned long long handle)
        {
            unsigned int numThreads = 0;
            cudaCGGetSize(&numThreads, NULL, handle);
            return numThreads;
        }

        _CG_STATIC_QUALIFIER unsigned int thread_rank(const unsigned long long handle)
        {
            unsigned int threadRank = 0;
            cudaCGGetRank(&threadRank, NULL, handle);
            return threadRank;
        }

        _CG_STATIC_QUALIFIER unsigned int grid_rank(const unsigned long long handle)
        {
            unsigned int gridRank = 0;
            cudaCGGetRank(NULL, &gridRank, handle);
            return gridRank;
        }

        _CG_STATIC_QUALIFIER unsigned int num_grids(const unsigned long long handle)
        {
            unsigned int numGrids = 0;
            cudaCGGetSize(NULL, &numGrids, handle);
            return numGrids;
        }

# ifdef _CG_CPP11_FEATURES
        struct multi_grid_functions {
            decltype(multi_grid::get_intrinsic_handle) *get_intrinsic_handle;
            decltype(multi_grid::sync) *sync;
            decltype(multi_grid::size) *size;
            decltype(multi_grid::thread_rank) *thread_rank;
            decltype(multi_grid::grid_rank) *grid_rank;
            decltype(multi_grid::num_grids) *num_grids;
        };

        template <typename = void>
        _CG_STATIC_QUALIFIER const multi_grid_functions* load_grid_intrinsics() {
            __constant__ static const multi_grid_functions mgf {
                &multi_grid::get_intrinsic_handle,
                &multi_grid::sync,
                &multi_grid::size,
                &multi_grid::thread_rank,
                &multi_grid::grid_rank,
                &multi_grid::num_grids
            };

            return &mgf;
        }
# endif
    };
#endif

    namespace cta {

        _CG_STATIC_QUALIFIER void sync()
        {
            __barrier_sync(0);
        }

        _CG_STATIC_QUALIFIER unsigned long long size()
        {
            return ((unsigned long long)blockDim.x * blockDim.y * blockDim.z);
        }

        _CG_STATIC_QUALIFIER unsigned long long thread_rank()
        {
            return vec3_to_linear(threadIdx, blockDim);
        }

        _CG_STATIC_QUALIFIER dim3 group_index()
        {
            return (dim3(blockIdx.x, blockIdx.y, blockIdx.z));
        }

        _CG_STATIC_QUALIFIER dim3 thread_index()
        {
            return (dim3(threadIdx.x, threadIdx.y, threadIdx.z));
        }

        _CG_STATIC_QUALIFIER dim3 block_dim()
        {
            return (dim3(blockDim.x, blockDim.y, blockDim.z));
        }

    };

    _CG_STATIC_QUALIFIER unsigned int laneid()
    {
        unsigned int laneid;
        asm volatile("mov.u32 %0, %%laneid;" : "=r"(laneid));
        return laneid;
    }

    _CG_STATIC_QUALIFIER unsigned int warpsz()
    {
        unsigned int warpSize;
        asm volatile("mov.u32 %0, WARP_SZ;" : "=r"(warpSize));
        return warpSize;
    }

    _CG_STATIC_QUALIFIER unsigned int lanemask32_eq()
    {
        unsigned int lanemask32_eq;
        asm volatile("mov.u32 %0, %%lanemask_eq;" : "=r"(lanemask32_eq));
        return (lanemask32_eq);
    }

    _CG_STATIC_QUALIFIER unsigned int lanemask32_lt()
    {
        unsigned int lanemask32_lt;
        asm volatile("mov.u32 %0, %%lanemask_lt;" : "=r"(lanemask32_lt));
        return (lanemask32_lt);
    }

    _CG_STATIC_QUALIFIER void abort()
    {
        _CG_ABORT();
    }

    template <typename Ty>
    _CG_QUALIFIER void assert_if_not_arithmetic() {
#ifdef _CG_CPP11_FEATURES
        static_assert(
            _CG_STL_NAMESPACE::is_integral<Ty>::value ||
            details::is_float_or_half<Ty>::value,
            "Error: Ty is neither integer or float"
        );
#endif
    }

}; // !Namespace internal

_CG_END_NAMESPACE

#endif /* !_COOPERATIVE_GROUPS_HELPERS_H_ */
