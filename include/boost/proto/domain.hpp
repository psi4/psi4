#ifndef BOOST_PP_IS_ITERATING
    ///////////////////////////////////////////////////////////////////////////////
    /// \file domain.hpp
    /// Contains definition of domain\<\> class template and helpers for
    /// defining domains with a generator and a grammar for controlling
    /// operator overloading.
    //
    //  Copyright 2008 Eric Niebler. Distributed under the Boost
    //  Software License, Version 1.0. (See accompanying file
    //  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

    #ifndef BOOST_PROTO_DOMAIN_HPP_EAN_02_13_2007
    #define BOOST_PROTO_DOMAIN_HPP_EAN_02_13_2007

    #include <boost/ref.hpp>
    #include <boost/mpl/bool.hpp>
    #include <boost/mpl/assert.hpp>
    #include <boost/preprocessor/cat.hpp>
    #include <boost/preprocessor/facilities/intercept.hpp>
    #include <boost/preprocessor/iteration/iterate.hpp>
    #include <boost/preprocessor/repetition/repeat.hpp>
    #include <boost/preprocessor/repetition/enum_trailing_params.hpp>
    #include <boost/preprocessor/repetition/enum_trailing_binary_params.hpp>
    #include <boost/preprocessor/repetition/enum_params.hpp>
    #include <boost/preprocessor/arithmetic/inc.hpp>
    #include <boost/preprocessor/arithmetic/dec.hpp>
    #include <boost/preprocessor/arithmetic/add.hpp>
    #include <boost/preprocessor/control/expr_if.hpp>
    #include <boost/proto/proto_fwd.hpp>
    #include <boost/proto/generate.hpp>

    #ifdef _MSC_VER
    #define BOOST_PROTO_DISABLE_MSVC_C4584 __pragma(warning(disable: 4584))
    #else
    #define BOOST_PROTO_DISABLE_MSVC_C4584 
    #endif

    namespace boost { namespace proto
    {

        namespace detail
        {
            struct not_a_generator
            {};

            struct not_a_grammar
            {};

            struct not_a_domain
            {};

            template<typename Super>
            struct super_domain : super_domain<typename Super::super>
            {
                typedef Super super;
                using super_domain<typename Super::super>::test;
                super_domain test(super_domain);
                super_domain test(super_domain<default_domain> const &);
                super *& get_super();
            };

            template<>
            struct super_domain<not_a_domain>
            {
                typedef not_a_domain super;
                super_domain test(...);
                super *& get_super();
            };

            template<>
            struct super_domain<default_domain>
            {
                typedef default_domain super;
                template<typename T> T test(T);
                super *& get_super();
            };

            template<typename T, int N> char (&select_domain(T*&))[N];
            template<typename T, int N> char (&select_domain(...))[1];

            template<
                int Index
                BOOST_PP_ENUM_TRAILING_BINARY_PARAMS(
                    BOOST_PROTO_MAX_ARITY
                  , typename D
                  , = void BOOST_PP_INTERCEPT
                )
            >
            struct select_nth
            {
                BOOST_MPL_ASSERT_MSG((false), PROTO_DOMAIN_MISMATCH, (select_nth));
                typedef not_a_domain type;
            };

            template<typename Void = void>
            struct deduce_domain0
            {
                typedef default_domain type;
            };

            #define BOOST_PP_ITERATION_PARAMS_1 (3, (1, BOOST_PROTO_MAX_ARITY, <boost/proto/domain.hpp>))
            #include BOOST_PP_ITERATE()
        }

        namespace domainns_
        {
            /// \brief For use in defining domain tags to be used
            /// with \c proto::extends\<\>. A \e Domain associates
            /// an expression type with a \e Generator, and optionally
            /// a \e Grammar.
            ///
            /// The Generator determines how new expressions in the
            /// domain are constructed. Typically, a generator wraps
            /// all new expressions in a wrapper that imparts
            /// domain-specific behaviors to expressions within its
            /// domain. (See \c proto::extends\<\>.)
            ///
            /// The Grammar determines whether a given expression is
            /// valid within the domain, and automatically disables
            /// any operator overloads which would cause an invalid
            /// expression to be created. By default, the Grammar
            /// parameter defaults to the wildcard, \c proto::_, which
            /// makes all expressions valid within the domain.
            ///
            /// Example:
            /// \code
            /// template<typename Expr>
            /// struct MyExpr;
            ///
            /// struct MyGrammar
            ///   : or_< terminal<_>, plus<MyGrammar, MyGrammar> >
            /// {};
            ///
            /// // Define MyDomain, in which all expressions are
            /// // wrapped in MyExpr<> and only expressions that
            /// // conform to MyGrammar are allowed.
            /// struct MyDomain
            ///   : domain<generator<MyExpr>, MyGrammar>
            /// {};
            ///
            /// // Use MyDomain to define MyExpr
            /// template<typename Expr>
            /// struct MyExpr
            ///   : extends<Expr, MyExpr<Expr>, MyDomain>
            /// {
            ///     // ...
            /// };
            /// \endcode
            ///
            BOOST_PROTO_DISABLE_MSVC_C4584
            template<
                typename Generator // = default_generator
              , typename Grammar   // = proto::_
              , typename Super     // = detail::not_a_domain
            >
            struct domain
              : detail::super_domain<Super>
              , Generator
            {
                typedef Generator proto_generator;
                typedef Grammar proto_grammar;

                /// INTERNAL ONLY
                typedef void proto_is_domain_;
            };

            /// \brief The domain expressions have by default, if
            /// \c proto::extends\<\> has not been used to associate
            /// a domain with an expression.
            ///
            struct default_domain
              : domain<>
            {};

            /// \brief A pseudo-domain for use in functions and
            /// metafunctions that require a domain parameter. It
            /// indicates that the domain of the parent node should
            /// be inferred from the domains of the child nodes.
            ///
            /// \attention \c deduce_domain is not itself a valid domain.
            ///
            struct deduce_domain
              : domain<detail::not_a_generator, detail::not_a_grammar, detail::not_a_domain>
            {};
        }

        namespace result_of
        {
            /// A metafunction that returns \c mpl::true_
            /// if the type \c T is the type of a Proto domain;
            /// \c mpl::false_ otherwise. If \c T inherits from
            /// \c proto::domain\<\>, \c is_domain\<T\> is
            /// \c mpl::true_.
            template<typename T, typename Void  /* = void*/>
            struct is_domain
              : mpl::false_
            {};

            /// INTERNAL ONLY
            ///
            template<typename T>
            struct is_domain<T, typename T::proto_is_domain_>
              : mpl::true_
            {};

            /// A metafunction that returns the domain of
            /// a given type. If \c T is a Proto expression
            /// type, it returns that expression's associated
            /// domain. If not, it returns
            /// \c proto::default_domain.
            template<typename T, typename Void /* = void*/>
            struct domain_of
            {
                typedef default_domain type;
            };

            /// INTERNAL ONLY
            ///
            template<typename T>
            struct domain_of<T, typename T::proto_is_expr_>
            {
                typedef typename T::proto_domain type;
            };

            /// INTERNAL ONLY
            ///
            template<typename T>
            struct domain_of<T &, void>
            {
                typedef typename domain_of<T>::type type;
            };

            /// INTERNAL ONLY
            ///
            template<typename T>
            struct domain_of<boost::reference_wrapper<T>, void>
            {
                typedef typename domain_of<T>::type type;
            };

            /// INTERNAL ONLY
            ///
            template<typename T>
            struct domain_of<boost::reference_wrapper<T> const, void>
            {
                typedef typename domain_of<T>::type type;
            };
        }
    }}

    #endif

#else

    #define N BOOST_PP_ITERATION()

            template<BOOST_PP_ENUM_PARAMS(BOOST_PROTO_MAX_ARITY, typename T)>
            struct select_nth<BOOST_PP_DEC(N), BOOST_PP_ENUM_PARAMS(BOOST_PROTO_MAX_ARITY, T)>
            {
                typedef BOOST_PP_CAT(T, BOOST_PP_DEC(N)) type;
            };

            template<BOOST_PP_ENUM_PARAMS(N, typename D)>
            struct BOOST_PP_CAT(common_domain, N)
            {
                #define M0(Z, M, DATA)                                                              \
                static detail::super_domain<BOOST_PP_CAT(D, M)> & BOOST_PP_CAT(d, M);               \
                /**/
                BOOST_PP_REPEAT(N, M0, ~)
                #undef M0

                enum e
                {
                    value0 = 0
                    #define M0(Z, M, DATA)                                                          \
                        BOOST_PP_EXPR_IF(M, .test)(BOOST_PP_CAT(d, M))                              \
                    /**/
                    #define M1(Z, M, DATA)                                                          \
                      , BOOST_PP_CAT(value, BOOST_PP_INC(M)) =                                      \
                        BOOST_PP_CAT(value, M) ?                                                    \
                        BOOST_PP_CAT(value, M) :                                                    \
                        sizeof(detail::select_domain<BOOST_PP_CAT(D, M), BOOST_PP_ADD(M, 2)>(       \
                            BOOST_PP_REPEAT_ ## Z(N, M0, ~).get_super()                             \
                        )) - 1                                                                      \
                    /**/
                    BOOST_PP_REPEAT(N, M1, ~)
                    #undef M1
                    #undef M0
                  , value = BOOST_PP_CAT(value, N) - 1
                };

                typedef typename select_nth<value BOOST_PP_ENUM_TRAILING_PARAMS(N, D)>::type type;
                //typedef BOOST_TYPEOF_TPL(d0.test(d1).test(d2).base()) type;
            };

            template<BOOST_PP_ENUM_PARAMS(N, typename E)>
            struct BOOST_PP_CAT(deduce_domain, N)
              : BOOST_PP_CAT(common_domain, N)<
                    BOOST_PP_ENUM_BINARY_PARAMS(N, typename domain_of<E, >::type BOOST_PP_INTERCEPT)
                >
            {};

    #undef N

#endif
