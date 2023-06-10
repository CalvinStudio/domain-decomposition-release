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

#ifndef _CG_REDUCE_H_
#define _CG_REDUCE_H_

#include "info.h"
#include "coalesced_reduce.h"
#include "functional.h"

_CG_BEGIN_NAMESPACE

namespace details {

    template <class Ty>
    using _redux_is_add_supported = _CG_STL_NAMESPACE::integral_constant<
            bool,
            _CG_STL_NAMESPACE::is_integral<Ty>::value && (sizeof(Ty) <= 4)>;

    template <class Ty>
    using redux_is_add_supported = _redux_is_add_supported<Ty>;

    // A specialization for 64 bit logical operations is possible
    // but for now only accelerate 32 bit bitwise ops
    template <class Ty>
    using redux_is_logical_supported = redux_is_add_supported<Ty>;

    // Base operator support case
    template <class TyOp, class Ty> struct _redux_op_supported                 : public _CG_STL_NAMESPACE::false_type {};
#ifdef _CG_HAS_OP_REDUX
    template <class Ty> struct _redux_op_supported<cooperative_groups::plus<Ty>, Ty>          : public redux_is_add_supported<Ty> {};
    template <class Ty> struct _redux_op_supported<cooperative_groups::less<Ty>, Ty>          : public redux_is_add_supported<Ty> {};
    template <class Ty> struct _redux_op_supported<cooperative_groups::greater<Ty>, Ty>       : public redux_is_add_supported<Ty> {};
    template <class Ty> struct _redux_op_supported<cooperative_groups::bit_and<Ty>, Ty>       : public redux_is_logical_supported<Ty> {};
    template <class Ty> struct _redux_op_supported<cooperative_groups::bit_or<Ty>, Ty>        : public redux_is_logical_supported<Ty> {};
    template <class Ty> struct _redux_op_supported<cooperative_groups::bit_xor<Ty>, Ty>       : public redux_is_logical_supported<Ty> {};
    template <class Ty> struct _redux_op_supported<const cooperative_groups::plus<Ty>, Ty>    : public redux_is_add_supported<Ty> {};
    template <class Ty> struct _redux_op_supported<const cooperative_groups::less<Ty>, Ty>    : public redux_is_add_supported<Ty> {};
    template <class Ty> struct _redux_op_supported<const cooperative_groups::greater<Ty>, Ty> : public redux_is_add_supported<Ty> {};
    template <class Ty> struct _redux_op_supported<const cooperative_groups::bit_and<Ty>, Ty> : public redux_is_logical_supported<Ty> {};
    template <class Ty> struct _redux_op_supported<const cooperative_groups::bit_or<Ty>, Ty>  : public redux_is_logical_supported<Ty> {};
    template <class Ty> struct _redux_op_supported<const cooperative_groups::bit_xor<Ty>, Ty> : public redux_is_logical_supported<Ty> {};
#endif

    template <template <class> class TyOp, class Ty>
    using redux_op_supported = _redux_op_supported<
            typename _CG_STL_NAMESPACE::remove_cv<TyOp<Ty>>::type,
            Ty>;

    // Group support for all reduce operations
    template <class TyGroup> struct _reduce_group_supported : public _CG_STL_NAMESPACE::false_type {};

    template <unsigned int Sz, typename TyPar>
    struct _reduce_group_supported<cooperative_groups::thread_block_tile<Sz, TyPar>> : public _CG_STL_NAMESPACE::true_type {};
    template <> struct _reduce_group_supported<cooperative_groups::coalesced_group>  : public _CG_STL_NAMESPACE::true_type {};

    template <typename TyGroup>
    using reduce_group_supported = _reduce_group_supported<_CG_STL_NAMESPACE::remove_cv_t<TyGroup>>;

    template <template <class> class TyOp>
    _CG_STATIC_QUALIFIER int pick_redux(int mask, int val);
    template <template <class> class TyOp>
    _CG_STATIC_QUALIFIER unsigned int pick_redux(int mask, unsigned int val);

#ifdef _CG_HAS_OP_REDUX
    template <> _CG_QUALIFIER int pick_redux<cooperative_groups::plus>(int mask, int val) {
        return __reduce_add_sync(mask, val);
    }
    template <> _CG_QUALIFIER int pick_redux<cooperative_groups::less>(int mask, int val) {
        return __reduce_min_sync(mask, val);
    }
    template <> _CG_QUALIFIER int pick_redux<cooperative_groups::greater>(int mask, int val) {
        return __reduce_max_sync(mask, val);
    }
    template <> _CG_QUALIFIER int pick_redux<cooperative_groups::bit_and>(int mask, int val) {
        return __reduce_and_sync(mask, val);
    }
    template <> _CG_QUALIFIER int pick_redux<cooperative_groups::bit_xor>(int mask, int val) {
        return __reduce_xor_sync(mask, val);
    }
    template <> _CG_QUALIFIER int pick_redux<cooperative_groups::bit_or>(int mask, int val) {
        return __reduce_or_sync(mask, val);
    }

    template <> _CG_QUALIFIER unsigned int pick_redux<cooperative_groups::plus>(int mask, unsigned int val) {
        return __reduce_add_sync(mask, val);
    }
    template <> _CG_QUALIFIER unsigned int pick_redux<cooperative_groups::less>(int mask, unsigned int val) {
        return __reduce_min_sync(mask, val);
    }
    template <> _CG_QUALIFIER unsigned int pick_redux<cooperative_groups::greater>(int mask, unsigned int val) {
        return __reduce_max_sync(mask, val);
    }
    template <> _CG_QUALIFIER unsigned int pick_redux<cooperative_groups::bit_and>(int mask, unsigned int val) {
        return __reduce_and_sync(mask, val);
    }
    template <> _CG_QUALIFIER unsigned int pick_redux<cooperative_groups::bit_xor>(int mask, unsigned int val) {
        return __reduce_xor_sync(mask, val);
    }
    template <> _CG_QUALIFIER unsigned int pick_redux<cooperative_groups::bit_or>(int mask, unsigned int val) {
        return __reduce_or_sync(mask, val);
    }
#endif


    template <typename TyVal, bool = _CG_STL_NAMESPACE::is_unsigned<TyVal>::value>
    struct _accelerated_op;

    // Signed type redux intrinsic dispatch
    template <typename TyVal>
    struct _accelerated_op<TyVal, false> {
        template <template <class> class TyOp>
        _CG_STATIC_QUALIFIER TyVal redux(int mask, TyVal val) {
            return static_cast<TyVal>(pick_redux<TyOp>(mask, static_cast<int>(val)));
        }
    };

    // Unsigned type redux intrinsic dispatch
    template <typename TyVal>
    struct _accelerated_op<TyVal, true> {
        template <template <class> class TyOp>
        _CG_STATIC_QUALIFIER TyVal redux(int mask, TyVal val) {
            return static_cast<TyVal>(pick_redux<TyOp>(mask, static_cast<unsigned int>(val)));
        }
    };

    template <typename TyVal>
    using accelerated_op = _accelerated_op<TyVal>;


    class _reduce {
        template <class Ty, template <class> class TyOp>
        using redux_is_usable = typename _CG_STL_NAMESPACE::enable_if<redux_op_supported<TyOp, Ty>::value, void>::type*;

        template <class Ty, template <class> class TyOp>
        using redux_is_not_usable = typename _CG_STL_NAMESPACE::enable_if<!redux_op_supported<TyOp, Ty>::value, void>::type*;

    public:
        // Dispatch to redux if the combination of op and args are supported
        template <
            class Ty,
            typename TyGroup,
            template <class> class TyOp,
            redux_is_usable<Ty, TyOp> = nullptr>
        _CG_STATIC_QUALIFIER Ty run(const TyGroup& group, Ty&& val, TyOp<Ty>&&) {
            // Retrieve the mask for the group and dispatch to redux
            return accelerated_op<Ty>::template redux<TyOp>(group.get_mask(), val);
        }

        // Fallback shuffle sync reduction
        template <
            class Ty,
            typename TyGroup,
            template <class> class TyOp,
            redux_is_not_usable<Ty, TyOp> = nullptr>
        _CG_STATIC_QUALIFIER Ty run(const TyGroup& group, Ty&& val, TyOp<Ty>&& op) {
            //Dispatch to fallback shuffle sync accelerated reduction
            return coalesced_reduce(group, _CG_STL_NAMESPACE::forward<Ty>(val), _CG_STL_NAMESPACE::forward<TyOp<Ty>>(op));
        }
    };

    template <typename TyLhs, typename TyRhs>
    using is_op_type_same = _CG_STL_NAMESPACE::is_same<
        typename _CG_STL_NAMESPACE::remove_reference<
            typename _CG_STL_NAMESPACE::remove_cv<TyLhs>::type
        >::type, 
        typename _CG_STL_NAMESPACE::remove_reference<
            typename _CG_STL_NAMESPACE::remove_cv<TyRhs>::type
        >::type
    >;
} // details

template <typename TyVal, typename TyArg, template <class> class TyOp, typename TyGroup>
_CG_QUALIFIER TyVal reduce(const TyGroup& group, TyArg&& val, TyOp<TyVal>&& op) {
    static_assert(details::reduce_group_supported<TyGroup>::value, "This group does not exclusively represent a tile");
    static_assert(details::is_op_type_same<TyVal, TyArg>::value, "Operator and argument types differ");

    return details::_reduce::run(group, _CG_STL_NAMESPACE::forward<TyVal>(val), _CG_STL_NAMESPACE::forward<TyOp<TyVal>>(op));
}

template <typename TyVal, typename TyFn, typename TyGroup, typename TyRet = _CG_STL_NAMESPACE::remove_reference_t<_CG_STL_NAMESPACE::remove_cv_t<TyVal>>>
_CG_QUALIFIER TyRet reduce(const TyGroup& group, TyVal&& val, TyFn&& op) {
    static_assert(details::reduce_group_supported<TyGroup>::value, "This group does not exclusively represent a tile");

    return details::coalesced_reduce(group, _CG_STL_NAMESPACE::forward<TyVal>(val), _CG_STL_NAMESPACE::forward<TyFn>(op));
}

_CG_END_NAMESPACE

#endif // _CG_REDUCE_H_
