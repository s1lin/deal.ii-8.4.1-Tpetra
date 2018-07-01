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
                    typename K0 = void_, typename K1 = void_, typename K2 = void_, typename K3 = void_, typename K4 = void_, typename K5 = void_, typename K6 = void_, typename K7 = void_, typename K8 = void_, typename K9 = void_, typename K10 = void_, typename K11 = void_, typename K12 = void_, typename K13 = void_, typename K14 = void_, typename K15 = void_, typename K16 = void_, typename K17 = void_, typename K18 = void_, typename K19 = void_, typename K20 = void_, typename K21 = void_, typename K22 = void_, typename K23 = void_, typename K24 = void_, typename K25 = void_, typename K26 = void_, typename K27 = void_, typename K28 = void_, typename K29 = void_, typename K30 = void_, typename K31 = void_, typename K32 = void_, typename K33 = void_, typename K34 = void_, typename K35 = void_, typename K36 = void_, typename K37 = void_, typename K38 = void_, typename K39 = void_, typename D0 = void_, typename D1 = void_, typename D2 = void_, typename D3 = void_, typename D4 = void_, typename D5 = void_, typename D6 = void_, typename D7 = void_, typename D8 = void_, typename D9 = void_, typename D10 = void_, typename D11 = void_, typename D12 = void_, typename D13 = void_, typename D14 = void_, typename D15 = void_, typename D16 = void_, typename D17 = void_, typename D18 = void_, typename D19 = void_, typename D20 = void_, typename D21 = void_, typename D22 = void_, typename D23 = void_, typename D24 = void_, typename D25 = void_, typename D26 = void_, typename D27 = void_, typename D28 = void_, typename D29 = void_, typename D30 = void_, typename D31 = void_, typename D32 = void_, typename D33 = void_, typename D34 = void_, typename D35 = void_, typename D36 = void_, typename D37 = void_, typename D38 = void_, typename D39 = void_, typename Extra = void_
            >
            struct make_map;
            template<>
            struct make_map<> {
                typedef map<> type;
            };
        }
        BOOST_FUSION_GPU_ENABLED inline map<>

        make_map() {
            return map<>();
        }

        namespace result_of {
            template<
                    typename K0, typename D0
            >
            struct make_map<K0, D0, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
                typedef map <fusion::pair<K0, typename detail::as_fusion_element<D0>::type>> type;
            };
        }
        template<
                typename K0, typename D0
        >
        BOOST_FUSION_GPU_ENABLED
        inline map<fusion::pair<K0, typename detail::as_fusion_element<D0>::type> >
        make_map(D0
        const& _0) {
        return

        map<fusion::pair<K0, typename detail::as_fusion_element<D0>::type> > (
        fusion::make_pair<K0>(_0));
    }
    namespace result_of {
        template<
                typename K0, typename K1, typename D0, typename D1
        >
        struct make_map<K0, K1, D0, D1, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
            typedef map <fusion::pair<K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>> type;
        };
    }
    template<
            typename K0, typename K1, typename D0, typename D1
    >
    BOOST_FUSION_GPU_ENABLED
    inline map<fusion::pair<K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type> >
    make_map(D0
    const& _0 ,
    D1 const &_1
    ) {
    return

    map<fusion::pair<K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type> >(
            fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1)

    );
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename D0, typename D1, typename D2
    >
    struct make_map<K0, K1, K2, D0, D1, D2, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map <fusion::pair<K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>> type;
    };
}
template<
        typename K0, typename K1, typename K2, typename D0, typename D1, typename D2
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename D0, typename D1, typename D2, typename D3
    >
    struct make_map<K0, K1, K2, K3, D0, D1, D2, D3, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename D0, typename D1, typename D2, typename D3
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(_3)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename D0, typename D1, typename D2, typename D3, typename D4
    >
    struct make_map<K0, K1, K2, K3, K4, D0, D1, D2, D3, D4, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename D0, typename D1, typename D2, typename D3, typename D4
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5
    >
    struct make_map<K0, K1, K2, K3, K4, K5, D0, D1, D2, D3, D4, D5, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, D0, D1, D2, D3, D4, D5, D6, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, D0, D1, D2, D3, D4, D5, D6, D7, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(_7)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, D0, D1, D2, D3, D4, D5, D6, D7, D8, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(_14)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(_17)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(_20)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(_23)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(_26)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>, fusion::pair<K29, typename detail::as_fusion_element<D29>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28 ,
D29 const &_29
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28), fusion::make_pair<K29>(_29)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>, fusion::pair<K29, typename detail::as_fusion_element<D29>::type>, fusion::pair<K30, typename detail::as_fusion_element<D30>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28 ,
D29 const &_29, D30
const& _30)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28), fusion::make_pair<K29>(
        _29), fusion::make_pair<K30>(_30)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>, fusion::pair<K29, typename detail::as_fusion_element<D29>::type>, fusion::pair<K30, typename detail::as_fusion_element<D30>::type>, fusion::pair<K31, typename detail::as_fusion_element<D31>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28 ,
D29 const &_29, D30
const& _30 ,
D31 const &_31
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28), fusion::make_pair<K29>(
        _29), fusion::make_pair<K30>(_30), fusion::make_pair<K31>(_31)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>, fusion::pair<K29, typename detail::as_fusion_element<D29>::type>, fusion::pair<K30, typename detail::as_fusion_element<D30>::type>, fusion::pair<K31, typename detail::as_fusion_element<D31>::type>, fusion::pair<K32, typename detail::as_fusion_element<D32>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28 ,
D29 const &_29, D30
const& _30 ,
D31 const &_31, D32
const& _32)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28), fusion::make_pair<K29>(
        _29), fusion::make_pair<K30>(_30), fusion::make_pair<K31>(_31), fusion::make_pair<K32>(_32)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>, fusion::pair<K29, typename detail::as_fusion_element<D29>::type>, fusion::pair<K30, typename detail::as_fusion_element<D30>::type>, fusion::pair<K31, typename detail::as_fusion_element<D31>::type>, fusion::pair<K32, typename detail::as_fusion_element<D32>::type>, fusion::pair<K33, typename detail::as_fusion_element<D33>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28 ,
D29 const &_29, D30
const& _30 ,
D31 const &_31, D32
const& _32 ,
D33 const &_33
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28), fusion::make_pair<K29>(
        _29), fusion::make_pair<K30>(_30), fusion::make_pair<K31>(_31), fusion::make_pair<K32>(
        _32), fusion::make_pair<K33>(_33)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>, fusion::pair<K29, typename detail::as_fusion_element<D29>::type>, fusion::pair<K30, typename detail::as_fusion_element<D30>::type>, fusion::pair<K31, typename detail::as_fusion_element<D31>::type>, fusion::pair<K32, typename detail::as_fusion_element<D32>::type>, fusion::pair<K33, typename detail::as_fusion_element<D33>::type>, fusion::pair<K34, typename detail::as_fusion_element<D34>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28 ,
D29 const &_29, D30
const& _30 ,
D31 const &_31, D32
const& _32 ,
D33 const &_33, D34
const& _34)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28), fusion::make_pair<K29>(
        _29), fusion::make_pair<K30>(_30), fusion::make_pair<K31>(_31), fusion::make_pair<K32>(
        _32), fusion::make_pair<K33>(_33), fusion::make_pair<K34>(_34)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, void_, void_, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> , fusion::pair<K35, typename detail::as_fusion_element<D35>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>, fusion::pair<K29, typename detail::as_fusion_element<D29>::type>, fusion::pair<K30, typename detail::as_fusion_element<D30>::type>, fusion::pair<K31, typename detail::as_fusion_element<D31>::type>, fusion::pair<K32, typename detail::as_fusion_element<D32>::type>, fusion::pair<K33, typename detail::as_fusion_element<D33>::type>, fusion::pair<K34, typename detail::as_fusion_element<D34>::type>, fusion::pair<K35, typename detail::as_fusion_element<D35>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28 ,
D29 const &_29, D30
const& _30 ,
D31 const &_31, D32
const& _32 ,
D33 const &_33, D34
const& _34 ,
D35 const &_35
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> , fusion::pair<K35, typename detail::as_fusion_element<D35>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28), fusion::make_pair<K29>(
        _29), fusion::make_pair<K30>(_30), fusion::make_pair<K31>(_31), fusion::make_pair<K32>(
        _32), fusion::make_pair<K33>(_33), fusion::make_pair<K34>(_34), fusion::make_pair<K35>(_35)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, void_, void_, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> , fusion::pair<K35, typename detail::as_fusion_element<D35>::type> , fusion::pair<K36, typename detail::as_fusion_element<D36>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>, fusion::pair<K29, typename detail::as_fusion_element<D29>::type>, fusion::pair<K30, typename detail::as_fusion_element<D30>::type>, fusion::pair<K31, typename detail::as_fusion_element<D31>::type>, fusion::pair<K32, typename detail::as_fusion_element<D32>::type>, fusion::pair<K33, typename detail::as_fusion_element<D33>::type>, fusion::pair<K34, typename detail::as_fusion_element<D34>::type>, fusion::pair<K35, typename detail::as_fusion_element<D35>::type>, fusion::pair<K36, typename detail::as_fusion_element<D36>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28 ,
D29 const &_29, D30
const& _30 ,
D31 const &_31, D32
const& _32 ,
D33 const &_33, D34
const& _34 ,
D35 const &_35, D36
const& _36)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> , fusion::pair<K35, typename detail::as_fusion_element<D35>::type> , fusion::pair<K36, typename detail::as_fusion_element<D36>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28), fusion::make_pair<K29>(
        _29), fusion::make_pair<K30>(_30), fusion::make_pair<K31>(_31), fusion::make_pair<K32>(
        _32), fusion::make_pair<K33>(_33), fusion::make_pair<K34>(_34), fusion::make_pair<K35>(
        _35), fusion::make_pair<K36>(_36)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, void_, void_, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> , fusion::pair<K35, typename detail::as_fusion_element<D35>::type> , fusion::pair<K36, typename detail::as_fusion_element<D36>::type> , fusion::pair<K37, typename detail::as_fusion_element<D37>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>, fusion::pair<K29, typename detail::as_fusion_element<D29>::type>, fusion::pair<K30, typename detail::as_fusion_element<D30>::type>, fusion::pair<K31, typename detail::as_fusion_element<D31>::type>, fusion::pair<K32, typename detail::as_fusion_element<D32>::type>, fusion::pair<K33, typename detail::as_fusion_element<D33>::type>, fusion::pair<K34, typename detail::as_fusion_element<D34>::type>, fusion::pair<K35, typename detail::as_fusion_element<D35>::type>, fusion::pair<K36, typename detail::as_fusion_element<D36>::type>, fusion::pair<K37, typename detail::as_fusion_element<D37>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28 ,
D29 const &_29, D30
const& _30 ,
D31 const &_31, D32
const& _32 ,
D33 const &_33, D34
const& _34 ,
D35 const &_35, D36
const& _36 ,
D37 const &_37
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> , fusion::pair<K35, typename detail::as_fusion_element<D35>::type> , fusion::pair<K36, typename detail::as_fusion_element<D36>::type> , fusion::pair<K37, typename detail::as_fusion_element<D37>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28), fusion::make_pair<K29>(
        _29), fusion::make_pair<K30>(_30), fusion::make_pair<K31>(_31), fusion::make_pair<K32>(
        _32), fusion::make_pair<K33>(_33), fusion::make_pair<K34>(_34), fusion::make_pair<K35>(
        _35), fusion::make_pair<K36>(_36), fusion::make_pair<K37>(_37)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, void_, void_, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> , fusion::pair<K35, typename detail::as_fusion_element<D35>::type> , fusion::pair<K36, typename detail::as_fusion_element<D36>::type> , fusion::pair<K37, typename detail::as_fusion_element<D37>::type> , fusion::pair<K38, typename detail::as_fusion_element<D38>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>, fusion::pair<K29, typename detail::as_fusion_element<D29>::type>, fusion::pair<K30, typename detail::as_fusion_element<D30>::type>, fusion::pair<K31, typename detail::as_fusion_element<D31>::type>, fusion::pair<K32, typename detail::as_fusion_element<D32>::type>, fusion::pair<K33, typename detail::as_fusion_element<D33>::type>, fusion::pair<K34, typename detail::as_fusion_element<D34>::type>, fusion::pair<K35, typename detail::as_fusion_element<D35>::type>, fusion::pair<K36, typename detail::as_fusion_element<D36>::type>, fusion::pair<K37, typename detail::as_fusion_element<D37>::type>, fusion::pair<K38, typename detail::as_fusion_element<D38>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28 ,
D29 const &_29, D30
const& _30 ,
D31 const &_31, D32
const& _32 ,
D33 const &_33, D34
const& _34 ,
D35 const &_35, D36
const& _36 ,
D37 const &_37, D38
const& _38)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> , fusion::pair<K35, typename detail::as_fusion_element<D35>::type> , fusion::pair<K36, typename detail::as_fusion_element<D36>::type> , fusion::pair<K37, typename detail::as_fusion_element<D37>::type> , fusion::pair<K38, typename detail::as_fusion_element<D38>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28), fusion::make_pair<K29>(
        _29), fusion::make_pair<K30>(_30), fusion::make_pair<K31>(_31), fusion::make_pair<K32>(
        _32), fusion::make_pair<K33>(_33), fusion::make_pair<K34>(_34), fusion::make_pair<K35>(
        _35), fusion::make_pair<K36>(_36), fusion::make_pair<K37>(_37), fusion::make_pair<K38>(_38)
);
}
namespace result_of {
    template<
            typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39
    >
    struct make_map<K0, K1, K2, K3, K4, K5, K6, K7, K8, K9, K10, K11, K12, K13, K14, K15, K16, K17, K18, K19, K20, K21, K22, K23, K24, K25, K26, K27, K28, K29, K30, K31, K32, K33, K34, K35, K36, K37, K38, K39, D0, D1, D2, D3, D4, D5, D6, D7, D8, D9, D10, D11, D12, D13, D14, D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25, D26, D27, D28, D29, D30, D31, D32, D33, D34, D35, D36, D37, D38, D39, void_> {
        typedef map<fusion::pair <
                    K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> , fusion::pair<K35, typename detail::as_fusion_element<D35>::type> , fusion::pair<K36, typename detail::as_fusion_element<D36>::type> , fusion::pair<K37, typename detail::as_fusion_element<D37>::type> , fusion::pair<K38, typename detail::as_fusion_element<D38>::type> , fusion::pair<K39, typename detail::as_fusion_element<D39>::type> > type;
    };
}
template<
        typename K0, typename K1, typename K2, typename K3, typename K4, typename K5, typename K6, typename K7, typename K8, typename K9, typename K10, typename K11, typename K12, typename K13, typename K14, typename K15, typename K16, typename K17, typename K18, typename K19, typename K20, typename K21, typename K22, typename K23, typename K24, typename K25, typename K26, typename K27, typename K28, typename K29, typename K30, typename K31, typename K32, typename K33, typename K34, typename K35, typename K36, typename K37, typename K38, typename K39, typename D0, typename D1, typename D2, typename D3, typename D4, typename D5, typename D6, typename D7, typename D8, typename D9, typename D10, typename D11, typename D12, typename D13, typename D14, typename D15, typename D16, typename D17, typename D18, typename D19, typename D20, typename D21, typename D22, typename D23, typename D24, typename D25, typename D26, typename D27, typename D28, typename D29, typename D30, typename D31, typename D32, typename D33, typename D34, typename D35, typename D36, typename D37, typename D38, typename D39
>
BOOST_FUSION_GPU_ENABLED
inline map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type>, fusion::pair<K1, typename detail::as_fusion_element<D1>::type>, fusion::pair<K2, typename detail::as_fusion_element<D2>::type>, fusion::pair<K3, typename detail::as_fusion_element<D3>::type>, fusion::pair<K4, typename detail::as_fusion_element<D4>::type>, fusion::pair<K5, typename detail::as_fusion_element<D5>::type>, fusion::pair<K6, typename detail::as_fusion_element<D6>::type>, fusion::pair<K7, typename detail::as_fusion_element<D7>::type>, fusion::pair<K8, typename detail::as_fusion_element<D8>::type>, fusion::pair<K9, typename detail::as_fusion_element<D9>::type>, fusion::pair<K10, typename detail::as_fusion_element<D10>::type>, fusion::pair<K11, typename detail::as_fusion_element<D11>::type>, fusion::pair<K12, typename detail::as_fusion_element<D12>::type>, fusion::pair<K13, typename detail::as_fusion_element<D13>::type>, fusion::pair<K14, typename detail::as_fusion_element<D14>::type>, fusion::pair<K15, typename detail::as_fusion_element<D15>::type>, fusion::pair<K16, typename detail::as_fusion_element<D16>::type>, fusion::pair<K17, typename detail::as_fusion_element<D17>::type>, fusion::pair<K18, typename detail::as_fusion_element<D18>::type>, fusion::pair<K19, typename detail::as_fusion_element<D19>::type>, fusion::pair<K20, typename detail::as_fusion_element<D20>::type>, fusion::pair<K21, typename detail::as_fusion_element<D21>::type>, fusion::pair<K22, typename detail::as_fusion_element<D22>::type>, fusion::pair<K23, typename detail::as_fusion_element<D23>::type>, fusion::pair<K24, typename detail::as_fusion_element<D24>::type>, fusion::pair<K25, typename detail::as_fusion_element<D25>::type>, fusion::pair<K26, typename detail::as_fusion_element<D26>::type>, fusion::pair<K27, typename detail::as_fusion_element<D27>::type>, fusion::pair<K28, typename detail::as_fusion_element<D28>::type>, fusion::pair<K29, typename detail::as_fusion_element<D29>::type>, fusion::pair<K30, typename detail::as_fusion_element<D30>::type>, fusion::pair<K31, typename detail::as_fusion_element<D31>::type>, fusion::pair<K32, typename detail::as_fusion_element<D32>::type>, fusion::pair<K33, typename detail::as_fusion_element<D33>::type>, fusion::pair<K34, typename detail::as_fusion_element<D34>::type>, fusion::pair<K35, typename detail::as_fusion_element<D35>::type>, fusion::pair<K36, typename detail::as_fusion_element<D36>::type>, fusion::pair<K37, typename detail::as_fusion_element<D37>::type>, fusion::pair<K38, typename detail::as_fusion_element<D38>::type>, fusion::pair<K39, typename detail::as_fusion_element<D39>::type>
>
make_map(D0
const& _0 ,
D1 const &_1, D2
const& _2 ,
D3 const &_3, D4
const& _4 ,
D5 const &_5, D6
const& _6 ,
D7 const &_7, D8
const& _8 ,
D9 const &_9, D10
const& _10 ,
D11 const &_11, D12
const& _12 ,
D13 const &_13, D14
const& _14 ,
D15 const &_15, D16
const& _16 ,
D17 const &_17, D18
const& _18 ,
D19 const &_19, D20
const& _20 ,
D21 const &_21, D22
const& _22 ,
D23 const &_23, D24
const& _24 ,
D25 const &_25, D26
const& _26 ,
D27 const &_27, D28
const& _28 ,
D29 const &_29, D30
const& _30 ,
D31 const &_31, D32
const& _32 ,
D33 const &_33, D34
const& _34 ,
D35 const &_35, D36
const& _36 ,
D37 const &_37, D38
const& _38 ,
D39 const &_39
)
{
return map<fusion::pair <
           K0, typename detail::as_fusion_element<D0>::type> , fusion::pair<K1, typename detail::as_fusion_element<D1>::type> , fusion::pair<K2, typename detail::as_fusion_element<D2>::type> , fusion::pair<K3, typename detail::as_fusion_element<D3>::type> , fusion::pair<K4, typename detail::as_fusion_element<D4>::type> , fusion::pair<K5, typename detail::as_fusion_element<D5>::type> , fusion::pair<K6, typename detail::as_fusion_element<D6>::type> , fusion::pair<K7, typename detail::as_fusion_element<D7>::type> , fusion::pair<K8, typename detail::as_fusion_element<D8>::type> , fusion::pair<K9, typename detail::as_fusion_element<D9>::type> , fusion::pair<K10, typename detail::as_fusion_element<D10>::type> , fusion::pair<K11, typename detail::as_fusion_element<D11>::type> , fusion::pair<K12, typename detail::as_fusion_element<D12>::type> , fusion::pair<K13, typename detail::as_fusion_element<D13>::type> , fusion::pair<K14, typename detail::as_fusion_element<D14>::type> , fusion::pair<K15, typename detail::as_fusion_element<D15>::type> , fusion::pair<K16, typename detail::as_fusion_element<D16>::type> , fusion::pair<K17, typename detail::as_fusion_element<D17>::type> , fusion::pair<K18, typename detail::as_fusion_element<D18>::type> , fusion::pair<K19, typename detail::as_fusion_element<D19>::type> , fusion::pair<K20, typename detail::as_fusion_element<D20>::type> , fusion::pair<K21, typename detail::as_fusion_element<D21>::type> , fusion::pair<K22, typename detail::as_fusion_element<D22>::type> , fusion::pair<K23, typename detail::as_fusion_element<D23>::type> , fusion::pair<K24, typename detail::as_fusion_element<D24>::type> , fusion::pair<K25, typename detail::as_fusion_element<D25>::type> , fusion::pair<K26, typename detail::as_fusion_element<D26>::type> , fusion::pair<K27, typename detail::as_fusion_element<D27>::type> , fusion::pair<K28, typename detail::as_fusion_element<D28>::type> , fusion::pair<K29, typename detail::as_fusion_element<D29>::type> , fusion::pair<K30, typename detail::as_fusion_element<D30>::type> , fusion::pair<K31, typename detail::as_fusion_element<D31>::type> , fusion::pair<K32, typename detail::as_fusion_element<D32>::type> , fusion::pair<K33, typename detail::as_fusion_element<D33>::type> , fusion::pair<K34, typename detail::as_fusion_element<D34>::type> , fusion::pair<K35, typename detail::as_fusion_element<D35>::type> , fusion::pair<K36, typename detail::as_fusion_element<D36>::type> , fusion::pair<K37, typename detail::as_fusion_element<D37>::type> , fusion::pair<K38, typename detail::as_fusion_element<D38>::type> , fusion::pair<K39, typename detail::as_fusion_element<D39>::type> >(
fusion::make_pair<K0>(_0), fusion::make_pair<K1>(_1), fusion::make_pair<K2>(_2), fusion::make_pair<K3>(
        _3), fusion::make_pair<K4>(_4), fusion::make_pair<K5>(_5), fusion::make_pair<K6>(_6), fusion::make_pair<K7>(
        _7), fusion::make_pair<K8>(_8), fusion::make_pair<K9>(_9), fusion::make_pair<K10>(_10), fusion::make_pair<K11>(
        _11), fusion::make_pair<K12>(_12), fusion::make_pair<K13>(_13), fusion::make_pair<K14>(
        _14), fusion::make_pair<K15>(_15), fusion::make_pair<K16>(_16), fusion::make_pair<K17>(
        _17), fusion::make_pair<K18>(_18), fusion::make_pair<K19>(_19), fusion::make_pair<K20>(
        _20), fusion::make_pair<K21>(_21), fusion::make_pair<K22>(_22), fusion::make_pair<K23>(
        _23), fusion::make_pair<K24>(_24), fusion::make_pair<K25>(_25), fusion::make_pair<K26>(
        _26), fusion::make_pair<K27>(_27), fusion::make_pair<K28>(_28), fusion::make_pair<K29>(
        _29), fusion::make_pair<K30>(_30), fusion::make_pair<K31>(_31), fusion::make_pair<K32>(
        _32), fusion::make_pair<K33>(_33), fusion::make_pair<K34>(_34), fusion::make_pair<K35>(
        _35), fusion::make_pair<K36>(_36), fusion::make_pair<K37>(_37), fusion::make_pair<K38>(
        _38), fusion::make_pair<K39>(_39)
);
}
}}
