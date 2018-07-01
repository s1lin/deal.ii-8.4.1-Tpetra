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
                    typename K0 = void_, typename K1 = void_, typename K2 = void_, typename K3 = void_, typename K4 = void_, typename K5 = void_, typename K6 = void_, typename K7 = void_, typename K8 = void_, typename K9 = void_, typename K10 = void_, typename K11 = void_, typename K12 = void_, typename K13 = void_, typename K14 = void_, typename K15 = void_, typename K16 = void_, typename K17 = void_, typename K18 = void_, typename K19 = void_, typename K20 = void_, typename K21 = void_, typename K22 = void_, typename K23 = void_, typename K24 = void_, typename K25 = void_, typename K26 = void_, typename K27 = void_, typename K28 = void_, typename K29 = void_, typename K30 = void_, typename K31 = void_, typename K32 = void_, typename K33 = void_, typename K34 = void_, typename K35 = void_, typename K36 = void_, typename K37 = void_, typename K38 = void_, typename K39 = void_, typename K40 = void_, typename K41 = void_, typename K42 = void_, typename K43 = void_, typename K44 = void_, typename K45 = void_, typename K46 = void_, typename K47 = void_, typename K48 = void_, typename K49 = void_, typename D0 = void_, typename D1 = void_, typename D2 = void_, typename D3 = void_, typename D4 = void_, typename D5 = void_, typename D6 = void_, typename D7 = void_, typename D8 = void_, typename D9 = void_, typename D10 = void_, typename D11 = void_, typename D12 = void_, typename D13 = void_, typename D14 = void_, typename D15 = void_, typename D16 = void_, typename D17 = void_, typename D18 = void_, typename D19 = void_, typename D20 = void_, typename D21 = void_, typename D22 = void_, typename D23 = void_, typename D24 = void_, typename D25 = void_, typename D26 = void_, typename D27 = void_, typename D28 = void_, typename D29 = void_, typename D30 = void_, typename D31 = void_, typename D32 = void_, typename D33 = void_, typename D34 = void_, typename D35 = void_, typename D36 = void_, typename D37 = void_, typename D38 = void_, typename D39 = void_, typename D40 = void_, typename D41 = void_, typename D42 = void_, typename D43 = void_, typename D44 = void_, typename D45 = void_, typename D46 = void_, typename D47 = void_, typename D48 = void_, typename D49 = void_, typename Extra = void_
            >
            struct map_tie;
            template<>
            struct map_tie<> {
                typedef map<> type;
            };
        }
        BOOST_FUSION_GPU_ENABLED inline map<>

        map_tie() {
            return map<>();
        }

        namespace result_of {
            template<
                    typename K0, typename D0
            >
            struct map_tie<K0, D0, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
                typedef map <fusion::pair<K0, typename add_reference<D0>::type>> type;
            };
        }
        template<
                typename K0, typename D0
        >
        BOOST_FUSION_GPU_ENABLED
        inline map<fusion::pair<K0, typename add_reference<D0>::type> >
        map_tie(D0
        & _0) {
        return

        map<fusion::pair<K0, typename add_reference<D0>::type> > (
        fusion::pair_tie<K0>(_0));
    }
    namespace result_of {
        template<
                typename K0, typename K1, typename D0, typename D1
        >
        struct map_tie<K0, K1, D0, D1, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
            typedef map <fusion::pair<K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>> type;
        };
    }
    template<
            typename K0, typename K1, typename D0, typename D1
    >
    BOOST_FUSION_GPU_ENABLED
    inline map<fusion::pair<K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type> >
    map_tie(D0
    & _0 ,
    D1 &_1
    ) {
    return

    map<fusion::pair<K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type> >(
            fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1)

    );
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename D0, typename D1, typename D2
    >
    struct map_tie<K0, K1, K2, D0, D1, D2, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map <fusion::pair<K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>> type;
    };
}
template<
        typename K0, typename K1, typename K2, typename D0, typename D1, typename D2
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename D0, typename D1, typename D2, typename D3
    >
    struct map_tie<K0, K1, K2, K3, D0, D1, D2, D3, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename D0, typename D1, typename D2, typename D3
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(_3)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename D0, typename D1, typename D2, typename D3, typename D4
    >
    struct map_tie<K0, K1, K2, K3, K4, D0, D1, D2, D3, D4, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename D0, typename D1, typename D2, typename D3, typename D4
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, D0, D1, D2, D3, D4, D5, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, D0, D1, D2, D3, D4, D5, D6, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, D0, D1, D2, D3, D4, D5, D6, D7, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(_7)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, D0, D1, D2, D3, D4, D5, D6, D7, D8, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(_11)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>, fusion::pair<K39, typename add_reference<D39>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38 ,
D39 &_39
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38), fusion::pair_tie<K39>(
        _39)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, K40, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, D40, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>, fusion::pair<K39, typename add_reference<D39>::type>, fusion::pair<K40, typename add_reference<D40>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38 ,
D39 &_39, D40
& _40)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38), fusion::pair_tie<K39>(
        _39), fusion::pair_tie<K40>(_40)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, K40, K41, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, D40, D41, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>, fusion::pair<K39, typename add_reference<D39>::type>, fusion::pair<K40, typename add_reference<D40>::type>, fusion::pair<K41, typename add_reference<D41>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38 ,
D39 &_39, D40
& _40 ,
D41 &_41
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38), fusion::pair_tie<K39>(
        _39), fusion::pair_tie<K40>(_40), fusion::pair_tie<K41>(_41)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, K40, K41, K42, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, D40, D41, D42, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>, fusion::pair<K39, typename add_reference<D39>::type>, fusion::pair<K40, typename add_reference<D40>::type>, fusion::pair<K41, typename add_reference<D41>::type>, fusion::pair<K42, typename add_reference<D42>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38 ,
D39 &_39, D40
& _40 ,
D41 &_41, D42
& _42)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38), fusion::pair_tie<K39>(
        _39), fusion::pair_tie<K40>(_40), fusion::pair_tie<K41>(_41), fusion::pair_tie<K42>(_42)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, K40, K41, K42, K43, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, D40, D41, D42, D43, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>, fusion::pair<K39, typename add_reference<D39>::type>, fusion::pair<K40, typename add_reference<D40>::type>, fusion::pair<K41, typename add_reference<D41>::type>, fusion::pair<K42, typename add_reference<D42>::type>, fusion::pair<K43, typename add_reference<D43>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38 ,
D39 &_39, D40
& _40 ,
D41 &_41, D42
& _42 ,
D43 &_43
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38), fusion::pair_tie<K39>(
        _39), fusion::pair_tie<K40>(_40), fusion::pair_tie<K41>(_41), fusion::pair_tie<K42>(_42), fusion::pair_tie<K43>(
        _43)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, K40, K41, K42, K43, K44, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, D40, D41, D42, D43, D44, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>, fusion::pair<K39, typename add_reference<D39>::type>, fusion::pair<K40, typename add_reference<D40>::type>, fusion::pair<K41, typename add_reference<D41>::type>, fusion::pair<K42, typename add_reference<D42>::type>, fusion::pair<K43, typename add_reference<D43>::type>, fusion::pair<K44, typename add_reference<D44>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38 ,
D39 &_39, D40
& _40 ,
D41 &_41, D42
& _42 ,
D43 &_43, D44
& _44)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38), fusion::pair_tie<K39>(
        _39), fusion::pair_tie<K40>(_40), fusion::pair_tie<K41>(_41), fusion::pair_tie<K42>(_42), fusion::pair_tie<K43>(
        _43), fusion::pair_tie<K44>(_44)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename K45, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44, typename D45
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, K40, K41, K42, K43, K44, K45, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, D40, D41, D42, D43, D44, D45, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> , fusion::pair<K45, typename add_reference<D45>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename K45, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44, typename D45
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>, fusion::pair<K39, typename add_reference<D39>::type>, fusion::pair<K40, typename add_reference<D40>::type>, fusion::pair<K41, typename add_reference<D41>::type>, fusion::pair<K42, typename add_reference<D42>::type>, fusion::pair<K43, typename add_reference<D43>::type>, fusion::pair<K44, typename add_reference<D44>::type>, fusion::pair<K45, typename add_reference<D45>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38 ,
D39 &_39, D40
& _40 ,
D41 &_41, D42
& _42 ,
D43 &_43, D44
& _44 ,
D45 &_45
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> , fusion::pair<K45, typename add_reference<D45>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38), fusion::pair_tie<K39>(
        _39), fusion::pair_tie<K40>(_40), fusion::pair_tie<K41>(_41), fusion::pair_tie<K42>(_42), fusion::pair_tie<K43>(
        _43), fusion::pair_tie<K44>(_44), fusion::pair_tie<K45>(_45)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename K45, typename K46, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44, typename D45, typename D46
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, K40, K41, K42, K43, K44, K45, K46, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, D40, D41, D42, D43, D44, D45, D46, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> , fusion::pair<K45, typename add_reference<D45>::type> , fusion::pair<K46, typename add_reference<D46>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename K45, typename K46, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44, typename D45, typename D46
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>, fusion::pair<K39, typename add_reference<D39>::type>, fusion::pair<K40, typename add_reference<D40>::type>, fusion::pair<K41, typename add_reference<D41>::type>, fusion::pair<K42, typename add_reference<D42>::type>, fusion::pair<K43, typename add_reference<D43>::type>, fusion::pair<K44, typename add_reference<D44>::type>, fusion::pair<K45, typename add_reference<D45>::type>, fusion::pair<K46, typename add_reference<D46>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38 ,
D39 &_39, D40
& _40 ,
D41 &_41, D42
& _42 ,
D43 &_43, D44
& _44 ,
D45 &_45, D46
& _46)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> , fusion::pair<K45, typename add_reference<D45>::type> , fusion::pair<K46, typename add_reference<D46>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38), fusion::pair_tie<K39>(
        _39), fusion::pair_tie<K40>(_40), fusion::pair_tie<K41>(_41), fusion::pair_tie<K42>(_42), fusion::pair_tie<K43>(
        _43), fusion::pair_tie<K44>(_44), fusion::pair_tie<K45>(_45), fusion::pair_tie<K46>(_46)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename K45, typename K46, typename K47, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44, typename D45, typename D46, typename D47
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, K40, K41, K42, K43, K44, K45, K46, K47, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, D40, D41, D42, D43, D44, D45, D46, D47, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> , fusion::pair<K45, typename add_reference<D45>::type> , fusion::pair<K46, typename add_reference<D46>::type> , fusion::pair<K47, typename add_reference<D47>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename K45, typename K46, typename K47, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44, typename D45, typename D46, typename D47
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>, fusion::pair<K39, typename add_reference<D39>::type>, fusion::pair<K40, typename add_reference<D40>::type>, fusion::pair<K41, typename add_reference<D41>::type>, fusion::pair<K42, typename add_reference<D42>::type>, fusion::pair<K43, typename add_reference<D43>::type>, fusion::pair<K44, typename add_reference<D44>::type>, fusion::pair<K45, typename add_reference<D45>::type>, fusion::pair<K46, typename add_reference<D46>::type>, fusion::pair<K47, typename add_reference<D47>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38 ,
D39 &_39, D40
& _40 ,
D41 &_41, D42
& _42 ,
D43 &_43, D44
& _44 ,
D45 &_45, D46
& _46 ,
D47 &_47
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> , fusion::pair<K45, typename add_reference<D45>::type> , fusion::pair<K46, typename add_reference<D46>::type> , fusion::pair<K47, typename add_reference<D47>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38), fusion::pair_tie<K39>(
        _39), fusion::pair_tie<K40>(_40), fusion::pair_tie<K41>(_41), fusion::pair_tie<K42>(_42), fusion::pair_tie<K43>(
        _43), fusion::pair_tie<K44>(_44), fusion::pair_tie<K45>(_45), fusion::pair_tie<K46>(_46), fusion::pair_tie<K47>(
        _47)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename K45, typename K46, typename K47, typename K48, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44, typename D45, typename D46, typename D47, typename D48
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, K40, K41, K42, K43, K44, K45, K46, K47, K48, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, D40, D41, D42, D43, D44, D45, D46, D47, D48, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> , fusion::pair<K45, typename add_reference<D45>::type> , fusion::pair<K46, typename add_reference<D46>::type> , fusion::pair<K47, typename add_reference<D47>::type> , fusion::pair<K48, typename add_reference<D48>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename K45, typename K46, typename K47, typename K48, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44, typename D45, typename D46, typename D47, typename D48
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>, fusion::pair<K39, typename add_reference<D39>::type>, fusion::pair<K40, typename add_reference<D40>::type>, fusion::pair<K41, typename add_reference<D41>::type>, fusion::pair<K42, typename add_reference<D42>::type>, fusion::pair<K43, typename add_reference<D43>::type>, fusion::pair<K44, typename add_reference<D44>::type>, fusion::pair<K45, typename add_reference<D45>::type>, fusion::pair<K46, typename add_reference<D46>::type>, fusion::pair<K47, typename add_reference<D47>::type>, fusion::pair<K48, typename add_reference<D48>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38 ,
D39 &_39, D40
& _40 ,
D41 &_41, D42
& _42 ,
D43 &_43, D44
& _44 ,
D45 &_45, D46
& _46 ,
D47 &_47, D48
& _48)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> , fusion::pair<K45, typename add_reference<D45>::type> , fusion::pair<K46, typename add_reference<D46>::type> , fusion::pair<K47, typename add_reference<D47>::type> , fusion::pair<K48, typename add_reference<D48>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38), fusion::pair_tie<K39>(
        _39), fusion::pair_tie<K40>(_40), fusion::pair_tie<K41>(_41), fusion::pair_tie<K42>(_42), fusion::pair_tie<K43>(
        _43), fusion::pair_tie<K44>(_44), fusion::pair_tie<K45>(_45), fusion::pair_tie<K46>(_46), fusion::pair_tie<K47>(
        _47), fusion::pair_tie<K48>(_48)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename K45, typename K46, typename K47, typename K48, typename K49, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44, typename D45, typename D46, typename D47, typename D48, typename D49
    >
    struct map_tie<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, K40, K41, K42, K43, K44, K45, K46, K47, K48, K49, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, D40, D41, D42, D43, D44, D45, D46, D47, D48, D49, void_> {
        typedef map<fusion::pair <
                    K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> , fusion::pair<K45, typename add_reference<D45>::type> , fusion::pair<K46, typename add_reference<D46>::type> , fusion::pair<K47, typename add_reference<D47>::type> , fusion::pair<K48, typename add_reference<D48>::type> , fusion::pair<K49, typename add_reference<D49>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename K40, typename K41, typename K42, typename K43, typename K44, typename K45, typename K46, typename K47, typename K48, typename K49, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39, typename D40, typename D41, typename D42, typename D43, typename D44, typename D45, typename D46, typename D47, typename D48, typename D49
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename add_reference<D0>::type>, fusion::pair<K1, typename add_reference<D1>::type>, fusion::pair<K2, typename add_reference<D2>::type>, fusion::pair<K3, typename add_reference<D3>::type>, fusion::pair<K4, typename add_reference<D4>::type>, fusion::pair<K5, typename add_reference<D5>::type>, fusion::pair<K6, typename add_reference<D6>::type>, fusion::pair<K7, typename add_reference<D7>::type>, fusion::pair<K8, typename add_reference<D8>::type>, fusion::pair<K9, typename add_reference<D9>::type>, fusion::pair<K10, typename add_reference<D10>::type>, fusion::pair<K11, typename add_reference<D11>::type>, fusion::pair<K12, typename add_reference<D12>::type>, fusion::pair<K13, typename add_reference<D13>::type>, fusion::pair<K14, typename add_reference<D14>::type>, fusion::pair<K15, typename add_reference<D15>::type>, fusion::pair<K16, typename add_reference<D16>::type>, fusion::pair<K17, typename add_reference<D17>::type>, fusion::pair<K18, typename add_reference<D18>::type>, fusion::pair<K19, typename add_reference<D19>::type>, fusion::pair<K20, typename add_reference<D20>::type>, fusion::pair<K21, typename add_reference<D21>::type>, fusion::pair<K22, typename add_reference<D22>::type>, fusion::pair<K23, typename add_reference<D23>::type>, fusion::pair<K24, typename add_reference<D24>::type>, fusion::pair<K25, typename add_reference<D25>::type>, fusion::pair<K26, typename add_reference<D26>::type>, fusion::pair<K27, typename add_reference<D27>::type>, fusion::pair<K28, typename add_reference<D28>::type>, fusion::pair<K29, typename add_reference<D29>::type>, fusion::pair<K30, typename add_reference<D30>::type>, fusion::pair<K31, typename add_reference<D31>::type>, fusion::pair<K32, typename add_reference<D32>::type>, fusion::pair<K33, typename add_reference<D33>::type>, fusion::pair<K34, typename add_reference<D34>::type>, fusion::pair<K35, typename add_reference<D35>::type>, fusion::pair<K36, typename add_reference<D36>::type>, fusion::pair<K37, typename add_reference<D37>::type>, fusion::pair<K38, typename add_reference<D38>::type>, fusion::pair<K39, typename add_reference<D39>::type>, fusion::pair<K40, typename add_reference<D40>::type>, fusion::pair<K41, typename add_reference<D41>::type>, fusion::pair<K42, typename add_reference<D42>::type>, fusion::pair<K43, typename add_reference<D43>::type>, fusion::pair<K44, typename add_reference<D44>::type>, fusion::pair<K45, typename add_reference<D45>::type>, fusion::pair<K46, typename add_reference<D46>::type>, fusion::pair<K47, typename add_reference<D47>::type>, fusion::pair<K48, typename add_reference<D48>::type>, fusion::pair<K49, typename add_reference<D49>::type>
>
map_tie(D0
& _0 ,
D1 &_1, D2
& _2 ,
D3 &_3, D4
& _4 ,
D5 &_5, D6
& _6 ,
D7 &_7, D8
& _8 ,
D9 &_9, D10
& _10 ,
D11 &_11, D12
& _12 ,
D13 &_13, D14
& _14 ,
D15 &_15, D16
& _16 ,
D17 &_17, D18
& _18 ,
D19 &_19, D20
& _20 ,
D21 &_21, D22
& _22 ,
D23 &_23, D24
& _24 ,
D25 &_25, D26
& _26 ,
D27 &_27, D28
& _28 ,
D29 &_29, D30
& _30 ,
D31 &_31, D32
& _32 ,
D33 &_33, D34
& _34 ,
D35 &_35, D36
& _36 ,
D37 &_37, D38
& _38 ,
D39 &_39, D40
& _40 ,
D41 &_41, D42
& _42 ,
D43 &_43, D44
& _44 ,
D45 &_45, D46
& _46 ,
D47 &_47, D48
& _48 ,
D49 &_49
)
{
return map<fusion::pair <
           K0, typename add_reference<D0>::type> , fusion::pair<K1, typename add_reference<D1>::type> , fusion::pair<K2, typename add_reference<D2>::type> , fusion::pair<K3, typename add_reference<D3>::type> , fusion::pair<K4, typename add_reference<D4>::type> , fusion::pair<K5, typename add_reference<D5>::type> , fusion::pair<K6, typename add_reference<D6>::type> , fusion::pair<K7, typename add_reference<D7>::type> , fusion::pair<K8, typename add_reference<D8>::type> , fusion::pair<K9, typename add_reference<D9>::type> , fusion::pair<K10, typename add_reference<D10>::type> , fusion::pair<K11, typename add_reference<D11>::type> , fusion::pair<K12, typename add_reference<D12>::type> , fusion::pair<K13, typename add_reference<D13>::type> , fusion::pair<K14, typename add_reference<D14>::type> , fusion::pair<K15, typename add_reference<D15>::type> , fusion::pair<K16, typename add_reference<D16>::type> , fusion::pair<K17, typename add_reference<D17>::type> , fusion::pair<K18, typename add_reference<D18>::type> , fusion::pair<K19, typename add_reference<D19>::type> , fusion::pair<K20, typename add_reference<D20>::type> , fusion::pair<K21, typename add_reference<D21>::type> , fusion::pair<K22, typename add_reference<D22>::type> , fusion::pair<K23, typename add_reference<D23>::type> , fusion::pair<K24, typename add_reference<D24>::type> , fusion::pair<K25, typename add_reference<D25>::type> , fusion::pair<K26, typename add_reference<D26>::type> , fusion::pair<K27, typename add_reference<D27>::type> , fusion::pair<K28, typename add_reference<D28>::type> , fusion::pair<K29, typename add_reference<D29>::type> , fusion::pair<K30, typename add_reference<D30>::type> , fusion::pair<K31, typename add_reference<D31>::type> , fusion::pair<K32, typename add_reference<D32>::type> , fusion::pair<K33, typename add_reference<D33>::type> , fusion::pair<K34, typename add_reference<D34>::type> , fusion::pair<K35, typename add_reference<D35>::type> , fusion::pair<K36, typename add_reference<D36>::type> , fusion::pair<K37, typename add_reference<D37>::type> , fusion::pair<K38, typename add_reference<D38>::type> , fusion::pair<K39, typename add_reference<D39>::type> , fusion::pair<K40, typename add_reference<D40>::type> , fusion::pair<K41, typename add_reference<D41>::type> , fusion::pair<K42, typename add_reference<D42>::type> , fusion::pair<K43, typename add_reference<D43>::type> , fusion::pair<K44, typename add_reference<D44>::type> , fusion::pair<K45, typename add_reference<D45>::type> , fusion::pair<K46, typename add_reference<D46>::type> , fusion::pair<K47, typename add_reference<D47>::type> , fusion::pair<K48, typename add_reference<D48>::type> , fusion::pair<K49, typename add_reference<D49>::type> >(
fusion::pair_tie<K0>(_0), fusion::pair_tie<K1>(_1), fusion::pair_tie<K2>(_2), fusion::pair_tie<K3>(
        _3), fusion::pair_tie<K4>(_4), fusion::pair_tie<K5>(_5), fusion::pair_tie<K6>(_6), fusion::pair_tie<K7>(
        _7), fusion::pair_tie<K8>(_8), fusion::pair_tie<K9>(_9), fusion::pair_tie<K10>(_10), fusion::pair_tie<K11>(
        _11), fusion::pair_tie<K12>(_12), fusion::pair_tie<K13>(_13), fusion::pair_tie<K14>(_14), fusion::pair_tie<K15>(
        _15), fusion::pair_tie<K16>(_16), fusion::pair_tie<K17>(_17), fusion::pair_tie<K18>(_18), fusion::pair_tie<K19>(
        _19), fusion::pair_tie<K20>(_20), fusion::pair_tie<K21>(_21), fusion::pair_tie<K22>(_22), fusion::pair_tie<K23>(
        _23), fusion::pair_tie<K24>(_24), fusion::pair_tie<K25>(_25), fusion::pair_tie<K26>(_26), fusion::pair_tie<K27>(
        _27), fusion::pair_tie<K28>(_28), fusion::pair_tie<K29>(_29), fusion::pair_tie<K30>(_30), fusion::pair_tie<K31>(
        _31), fusion::pair_tie<K32>(_32), fusion::pair_tie<K33>(_33), fusion::pair_tie<K34>(_34), fusion::pair_tie<K35>(
        _35), fusion::pair_tie<K36>(_36), fusion::pair_tie<K37>(_37), fusion::pair_tie<K38>(_38), fusion::pair_tie<K39>(
        _39), fusion::pair_tie<K40>(_40), fusion::pair_tie<K41>(_41), fusion::pair_tie<K42>(_42), fusion::pair_tie<K43>(
        _43), fusion::pair_tie<K44>(_44), fusion::pair_tie<K45>(_45), fusion::pair_tie<K46>(_46), fusion::pair_tie<K47>(
        _47), fusion::pair_tie<K48>(_48), fusion::pair_tie<K49>(_49)
);
}
}}
