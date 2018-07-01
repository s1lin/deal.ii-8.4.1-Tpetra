/*=============================================================================
    Copyright (c) 2001-2011 Joel de Guzman

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

    This is an auto-generated file. Do not edit!
==============================================================================*/
namespace boost {
    namespace fusion {
        struct void_;
        namespace result_of {
            template<
                    typename T0 = void_, typename T1 = void_, typename T2 = void_, typename T3 = void_, typename T4 = void_, typename T5 = void_, typename T6 = void_, typename T7 = void_, typename T8 = void_, typename T9 = void_, typename T10 = void_, typename T11 = void_, typename T12 = void_, typename T13 = void_, typename T14 = void_, typename T15 = void_, typename T16 = void_, typename T17 = void_, typename T18 = void_, typename T19 = void_, typename T20 = void_, typename T21 = void_, typename T22 = void_, typename T23 = void_, typename T24 = void_, typename T25 = void_, typename T26 = void_, typename T27 = void_, typename T28 = void_, typename T29 = void_, typename T30 = void_, typename T31 = void_, typename T32 = void_, typename T33 = void_, typename T34 = void_, typename T35 = void_, typename T36 = void_, typename T37 = void_, typename T38 = void_, typename T39 = void_, typename Extra = void_
            >
            struct make_deque;
            template<>
            struct make_deque<> {
                typedef deque<> type;
            };
        }
        BOOST_FUSION_GPU_ENABLED inline deque<>

        make_deque() {
            return deque<>();
        }

        namespace result_of {
            template<typename T0>
            struct make_deque<T0, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
                typedef deque<typename detail::as_fusion_element<T0>::type> type;
            };
        }
        template<typename T0>
        BOOST_FUSION_GPU_ENABLED
        inline deque<typename detail::as_fusion_element<T0>::type>
        make_deque(T0
        const& _0) {
        return
        deque<typename detail::as_fusion_element<T0>::type>(
                _0);
    }
    namespace result_of {
        template<typename T0, typename T1>
        struct make_deque<T0, T1, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
            typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type> type;
        };
    }
    template<typename T0, typename T1>
    BOOST_FUSION_GPU_ENABLED
    inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type>
    make_deque(T0
    const& _0 ,
    T1 const &_1
    ) {
    return
    deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type>(
            _0, _1
    );
}
namespace result_of {
    template<typename T0, typename T1, typename T2>
    struct make_deque<T0, T1, T2, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type> type;
    };
}
template<typename T0, typename T1, typename T2>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type>(
        _0, _1, _2
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3>
    struct make_deque<T0, T1, T2, T3, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type>(
        _0, _1, _2, _3
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4>
    struct make_deque<T0, T1, T2, T3, T4, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type>(
        _0, _1, _2, _3, _4
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
    struct make_deque<T0, T1, T2, T3, T4, T5, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type>(
        _0, _1, _2, _3, _4, _5
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type>(
        _0, _1, _2, _3, _4, _5, _6
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28 ,
T29 const &_29
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28 ,
T29 const &_29, T30
const& _30)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28 ,
T29 const &_29, T30
const& _30 ,
T31 const &_31
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28 ,
T29 const &_29, T30
const& _30 ,
T31 const &_31, T32
const& _32)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, T33, void_, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28 ,
T29 const &_29, T30
const& _30 ,
T31 const &_31, T32
const& _32 ,
T33 const &_33
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, T33, T34, void_, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28 ,
T29 const &_29, T30
const& _30 ,
T31 const &_31, T32
const& _32 ,
T33 const &_33, T34
const& _34)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34, typename T35>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, T33, T34, T35, void_, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34, typename T35>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28 ,
T29 const &_29, T30
const& _30 ,
T31 const &_31, T32
const& _32 ,
T33 const &_33, T34
const& _34 ,
T35 const &_35
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34, _35
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34, typename T35, typename T36>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, T33, T34, T35, T36, void_, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34, typename T35, typename T36>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28 ,
T29 const &_29, T30
const& _30 ,
T31 const &_31, T32
const& _32 ,
T33 const &_33, T34
const& _34 ,
T35 const &_35, T36
const& _36)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34, _35, _36
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, T33, T34, T35, T36, T37, void_, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type, typename detail::as_fusion_element<T37>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type, typename detail::as_fusion_element<T37>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28 ,
T29 const &_29, T30
const& _30 ,
T31 const &_31, T32
const& _32 ,
T33 const &_33, T34
const& _34 ,
T35 const &_35, T36
const& _36 ,
T37 const &_37
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type, typename detail::as_fusion_element<T37>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34, _35, _36, _37
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, T33, T34, T35, T36, T37, T38, void_, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type, typename detail::as_fusion_element<T37>::type, typename detail::as_fusion_element<T38>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type, typename detail::as_fusion_element<T37>::type, typename detail::as_fusion_element<T38>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28 ,
T29 const &_29, T30
const& _30 ,
T31 const &_31, T32
const& _32 ,
T33 const &_33, T34
const& _34 ,
T35 const &_35, T36
const& _36 ,
T37 const &_37, T38
const& _38)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type, typename detail::as_fusion_element<T37>::type, typename detail::as_fusion_element<T38>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34, _35, _36, _37, _38
);
}
namespace result_of {
    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39>
    struct make_deque<T0, T1, T2, T3, T4, T5, T6, T7, T8, T9, T10, T11, T12, T13, T14, T15, T16, T17, T18, T19, T20, T21, T22, T23, T24, T25, T26, T27, T28, T29, T30, T31, T32, T33, T34, T35, T36, T37, T38, T39, void_> {
        typedef deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type, typename detail::as_fusion_element<T37>::type, typename detail::as_fusion_element<T38>::type, typename detail::as_fusion_element<T39>::type> type;
    };
}
template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5, typename T6, typename T7, typename T8, typename T9, typename T10, typename T11, typename T12, typename T13, typename T14, typename T15, typename T16, typename T17, typename T18, typename T19, typename T20, typename T21, typename T22, typename T23, typename T24, typename T25, typename T26, typename T27, typename T28, typename T29, typename T30, typename T31, typename T32, typename T33, typename T34, typename T35, typename T36, typename T37, typename T38, typename T39>
BOOST_FUSION_GPU_ENABLED
inline deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type, typename detail::as_fusion_element<T37>::type, typename detail::as_fusion_element<T38>::type, typename detail::as_fusion_element<T39>::type>
make_deque(T0
const& _0 ,
T1 const &_1, T2
const& _2 ,
T3 const &_3, T4
const& _4 ,
T5 const &_5, T6
const& _6 ,
T7 const &_7, T8
const& _8 ,
T9 const &_9, T10
const& _10 ,
T11 const &_11, T12
const& _12 ,
T13 const &_13, T14
const& _14 ,
T15 const &_15, T16
const& _16 ,
T17 const &_17, T18
const& _18 ,
T19 const &_19, T20
const& _20 ,
T21 const &_21, T22
const& _22 ,
T23 const &_23, T24
const& _24 ,
T25 const &_25, T26
const& _26 ,
T27 const &_27, T28
const& _28 ,
T29 const &_29, T30
const& _30 ,
T31 const &_31, T32
const& _32 ,
T33 const &_33, T34
const& _34 ,
T35 const &_35, T36
const& _36 ,
T37 const &_37, T38
const& _38 ,
T39 const &_39
)
{
return
deque<typename detail::as_fusion_element<T0>::type, typename detail::as_fusion_element<T1>::type, typename detail::as_fusion_element<T2>::type, typename detail::as_fusion_element<T3>::type, typename detail::as_fusion_element<T4>::type, typename detail::as_fusion_element<T5>::type, typename detail::as_fusion_element<T6>::type, typename detail::as_fusion_element<T7>::type, typename detail::as_fusion_element<T8>::type, typename detail::as_fusion_element<T9>::type, typename detail::as_fusion_element<T10>::type, typename detail::as_fusion_element<T11>::type, typename detail::as_fusion_element<T12>::type, typename detail::as_fusion_element<T13>::type, typename detail::as_fusion_element<T14>::type, typename detail::as_fusion_element<T15>::type, typename detail::as_fusion_element<T16>::type, typename detail::as_fusion_element<T17>::type, typename detail::as_fusion_element<T18>::type, typename detail::as_fusion_element<T19>::type, typename detail::as_fusion_element<T20>::type, typename detail::as_fusion_element<T21>::type, typename detail::as_fusion_element<T22>::type, typename detail::as_fusion_element<T23>::type, typename detail::as_fusion_element<T24>::type, typename detail::as_fusion_element<T25>::type, typename detail::as_fusion_element<T26>::type, typename detail::as_fusion_element<T27>::type, typename detail::as_fusion_element<T28>::type, typename detail::as_fusion_element<T29>::type, typename detail::as_fusion_element<T30>::type, typename detail::as_fusion_element<T31>::type, typename detail::as_fusion_element<T32>::type, typename detail::as_fusion_element<T33>::type, typename detail::as_fusion_element<T34>::type, typename detail::as_fusion_element<T35>::type, typename detail::as_fusion_element<T36>::type, typename detail::as_fusion_element<T37>::type, typename detail::as_fusion_element<T38>::type, typename detail::as_fusion_element<T39>::type>(
        _0, _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12, _13, _14, _15, _16, _17, _18, _19, _20, _21, _22, _23, _24, _25, _26, _27, _28, _29, _30, _31, _32, _33, _34, _35, _36, _37, _38, _39
);
}
}}
