/*=============================================================================
    Copyright (c) 2001-2014 Joel de Guzman
    Copyright (c) 2001-2011 Hartmut Kaiser

    Distributed under the Boost Software License, Version 1.0. (See accompanying
    file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
=============================================================================*/
#if !defined(SPIRIT_LIST_MARCH_24_2007_1031AM)
#define SPIRIT_LIST_MARCH_24_2007_1031AM

#if defined(_MSC_VER)
#pragma once
#endif

#include <boost/spirit/home/x3/core/parser.hpp>
#include <boost/spirit/home/x3/support/traits/container_traits.hpp>
#include <boost/spirit/home/x3/support/traits/attribute_of.hpp>
#include <boost/spirit/home/x3/core/detail/parse_into_container.hpp>

namespace boost {
    namespace spirit {
        namespace x3 {
            template<typename Left, typename Right>
            struct list : binary_parser<Left, Right, list<Left, Right>> {
                typedef binary_parser <Left, Right, list<Left, Right>> base_type;
                static bool const handles_container = true;
                static bool const has_attribute = true;

                list(Left const &left, Right const &right)
                        : base_type(left, right) {}

                template<typename Iterator, typename Context, typename RContext, typename Attribute>
                bool parse(Iterator &first, Iterator const &last, Context const &context, RContext &rcontext,
                           Attribute &attr) const {
                    // in order to succeed we need to match at least one element
                    if (!detail::parse_into_container(
                            this->left, first, last, context, rcontext, attr))
                        return false;

                    Iterator save = first;
                    while (this->right.parse(first, last, context, rcontext, unused)
                           && detail::parse_into_container(
                            this->left, first, last, context, rcontext, attr)) {
                        save = first;
                    }

                    first = save;
                    return true;
                }
            };

            template<typename Left, typename Right>
            inline list<
                    typename extension::as_parser<Left>::value_type, typename extension::as_parser<Right>::value_type>
            operator%(Left const &left, Right const &right) {
                return {as_parser(left), as_parser(right)};
            }
        }
    }
}

namespace boost {
    namespace spirit {
        namespace x3 {
            namespace traits {
                template<typename Left, typename Right, typename Context>
                struct attribute_of<x3::list<Left, Right>, Context>
                        : traits::build_container<
                                typename attribute_of<Left, Context>::type> {
                };
            }
        }
    }
}

#endif
