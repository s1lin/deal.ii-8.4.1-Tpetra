
// Copyright (C) 2008-2011 Daniel James.
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_UNORDERED_FWD_HPP_INCLUDED
#define BOOST_UNORDERED_FWD_HPP_INCLUDED

#include <boost/config.hpp>

#if defined(BOOST_HAS_PRAGMA_ONCE)
#pragma once
#endif

namespace boost {
    namespace unordered {
        struct piecewise_construct_t {
        };
        const piecewise_construct_t piecewise_construct = piecewise_construct_t();
    }
}

#endif
