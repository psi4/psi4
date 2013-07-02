/*
 This file is part of MADNESS.

 Copyright (C) 2007,2010 Oak Ridge National Laboratory

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

 For more information please contact:

 Robert J. Harrison
 Oak Ridge National Laboratory
 One Bethel Valley Road
 P.O. Box 2008, MS-6367

 email: harrisonrj@ornl.gov
 tel:   865-241-3937
 fax:   865-572-0680


 $Id $
 */

#ifndef MADNESS_WORLD_TASKFN_H__INCLUDED
#define MADNESS_WORLD_TASKFN_H__INCLUDED

//namespace madness {

    namespace {

        template <typename T>
        struct ArgCountHelper {
            static const unsigned int value = 1;
        };

        template <>
        struct ArgCountHelper<void> {
            static const unsigned int value = 0;
        };

        // Counts the number of arguments that will be given to a task function
        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T, typename a7T, typename a8T, typename a9T>
        struct ArgCount : public std::integral_constant<unsigned int, ArgCountHelper<a1T>::value
                + ArgCountHelper<a2T>::value + ArgCountHelper<a3T>::value
                + ArgCountHelper<a4T>::value + ArgCountHelper<a5T>::value
                + ArgCountHelper<a6T>::value + ArgCountHelper<a7T>::value
                + ArgCountHelper<a8T>::value + ArgCountHelper<a9T>::value> {
        };

        template <typename T>
        struct value_type {
            typedef typename std::remove_cv<T>::type type;
        };

        template <typename T>
        struct value_type<Future<T> > {
            typedef T type;
        };

        template <typename T>
        struct value_type<T&> {
            typedef typename value_type<T>::type type;
        };

        template <typename T>
        struct value_type<T*> {
            typedef T* type;
        };

        template <typename T, std::size_t N>
        struct value_type<T[N]> {
            typedef T* type;
        };

        template <typename T, std::size_t N>
        struct value_type<const T[N]> {
            typedef const T* type;
        };

        template <std::size_t N>
        struct value_type<char[N]> {
            typedef std::string type;
        };

        template <std::size_t N>
        struct value_type<const char[N]> {
            typedef std::string type;
        };

        template <>
        struct value_type<void> {
            typedef void type;
        };

        template <typename fnT>
        struct task_result_type {
            typedef typename remove_fcvr<typename result_of<fnT>::type>::type resultT;
            typedef Future<resultT> futureT;
        };

    } // namespace

    /// Wrap a callable object and its arguments into a task function

    /// The callable object may have up to 10 arguments
    template <typename fnT, typename arg1_type = void, typename arg2_type = void,
            typename arg3_type = void, typename arg4_type = void, typename arg5_type = void,
            typename arg6_type = void, typename arg7_type = void, typename arg8_type = void,
            typename arg9_type = void>
    struct TaskFn : public TaskInterface {
    private:
        /// This class type
        typedef TaskFn<fnT, arg1_type, arg2_type, arg3_type, arg4_type, arg5_type,
                arg6_type, arg7_type, arg8_type, arg9_type> TaskFn_;

    public:

        typedef fnT functionT; ///< The task function type
        typedef typename task_result_type<fnT>::resultT resultT;
        ///< The result type of the function
        typedef typename task_result_type<fnT>::futureT futureT;

        // argument value typedefs
        typedef typename value_type<arg1_type>::type arg1T; ///< Argument 1 type
        typedef typename value_type<arg2_type>::type arg2T; ///< Argument 2 type
        typedef typename value_type<arg3_type>::type arg3T; ///< Argument 3 type
        typedef typename value_type<arg4_type>::type arg4T; ///< Argument 4 type
        typedef typename value_type<arg5_type>::type arg5T; ///< Argument 5 type
        typedef typename value_type<arg6_type>::type arg6T; ///< Argument 6 type
        typedef typename value_type<arg7_type>::type arg7T; ///< Argument 7 type
        typedef typename value_type<arg8_type>::type arg8T; ///< Argument 8 type
        typedef typename value_type<arg9_type>::type arg9T; ///< Argument 9 type

        static const unsigned int arity = ArgCount<arg1_type, arg2_type, arg3_type, arg4_type,
                arg5_type, arg6_type, arg7_type, arg8_type, arg9_type>::value;
        ///< The number of arguments given for the function
        ///< \note This may not match the arity of the function
        ///< if it has default parameter values

    private:
        futureT result_; ///< The task Future result
        const functionT func_; ///< The task function

        // If the value of the argument is known at the time the
        // Note: The type argNT for argN, where N  is > arity should be void

        Future<arg1T> arg1_;///< Argument 1 that will be given to the function
        Future<arg2T> arg2_;///< Argument 2 that will be given to the function
        Future<arg3T> arg3_;///< Argument 3 that will be given to the function
        Future<arg4T> arg4_;///< Argument 4 that will be given to the function
        Future<arg5T> arg5_;///< Argument 5 that will be given to the function
        Future<arg6T> arg6_;///< Argument 6 that will be given to the function
        Future<arg7T> arg7_;///< Argument 7 that will be given to the function
        Future<arg8T> arg8_;///< Argument 8 that will be given to the function
        Future<arg9T> arg9_;///< Argument 9 that will be given to the function

        // These functions are here because we have to differentiate the call
        // based on the number of arguments passed to the function and the
        // return type.

        template <unsigned int N, typename R>
        typename enable_if_c<(!std::is_void<R>::value) && (N == 0u)>::type runner() {
            result_.set(func_());
        }

        template <unsigned int N, typename R>
        typename enable_if_c<(!std::is_void<R>::value) && (N == 1u)>::type runner() {
            result_.set(func_(arg1_));
        }

        template <unsigned int N, typename R>
        typename enable_if_c<(!std::is_void<R>::value) && (N == 2u)>::type runner() {
            result_.set(func_(arg1_, arg2_));
        }

        template <unsigned int N, typename R>
        typename enable_if_c<(!std::is_void<R>::value) && (N == 3u)>::type runner() {
            result_.set(func_(arg1_, arg2_, arg3_));
        }

        template <unsigned int N, typename R>
        typename enable_if_c<(!std::is_void<R>::value) && (N == 4u)>::type runner() {
            result_.set(func_(arg1_, arg2_, arg3_, arg4_));
        }

        template <unsigned int N, typename R>
        typename enable_if_c<(!std::is_void<R>::value) && (N == 5u)>::type runner() {
            result_.set(func_(arg1_, arg2_, arg3_, arg4_, arg5_));
        }

        template <unsigned int N, typename R>
        typename enable_if_c<(!std::is_void<R>::value) && (N == 6u)>::type runner() {
            result_.set(func_(arg1_, arg2_, arg3_, arg4_, arg5_, arg6_));
        }

        template <unsigned int N, typename R>
        typename enable_if_c<(!std::is_void<R>::value) && (N == 7u)>::type runner() {
            result_.set(func_(arg1_, arg2_, arg3_, arg4_, arg5_, arg6_, arg7_));
        }

        template <unsigned int N, typename R>
        typename enable_if_c<(!std::is_void<R>::value) && (N == 8u)>::type runner() {
            result_.set(func_(arg1_, arg2_, arg3_, arg4_, arg5_, arg6_, arg7_, arg8_));
        }

        template <unsigned int N, typename R>
        typename enable_if_c<(!std::is_void<R>::value) && (N == 9u)>::type runner() {
            result_.set(func_(arg1_, arg2_, arg3_, arg4_, arg5_, arg6_, arg7_, arg8_, arg9_));
        }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 0u)>::type runner() {
            result_.set(func_());
        }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 1u)>::type runner() {
            func_(arg1_);
            result_.set();
        }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 2u)>::type runner() {
            func_(arg1_, arg2_);
            result_.set();
        }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 3u)>::type runner() {
            func_(arg1_, arg2_, arg3_);
            result_.set();
        }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 4u)>::type runner() {
            func_(arg1_, arg2_, arg3_, arg4_);
            result_.set();
        }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 5u)>::type runner() {
            func_(arg1_, arg2_, arg3_, arg4_, arg5_);
            result_.set();
        }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 6u)>::type runner() {
            func_(arg1_, arg2_, arg3_, arg4_, arg5_, arg6_);
            result_.set();
        }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 7u)>::type runner() {
            func_(arg1_, arg2_, arg3_, arg4_, arg5_, arg6_, arg7_);
            result_.set();
        }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 8u)>::type runner() {
            func_(arg1_, arg2_, arg3_, arg4_, arg5_, arg6_, arg7_, arg8_);
            result_.set();
        }

        template <unsigned int N, typename R>
        typename enable_if_c<std::is_void<R>::value && (N == 9u)>::type runner() {
            func_(arg1_, arg2_, arg3_, arg4_, arg5_, arg6_, arg7_, arg8_, arg9_);
            result_.set();
        }

        /// Register non-ready future as a dependency

        /// \tparam T The type of the future to check
        /// \param fut The future to check
        template <typename T>
        inline void check_dependency(Future<T>& fut) {
            if(!fut.probe()) {
                DependencyInterface::inc();
                fut.register_callback(this);
            }
        }

        /// Future<void> is always ready => no op
        inline void check_dependency(Future<void>&) { }

        /// Future<void> is always ready => no op
        inline void check_dependency(Future<Void>&) { }

        /// Check dependencies and register callbacks where necessary
        void check_dependencies() {
            // Check only what we need
            switch(arity) {
                case 9u:
                    check_dependency(arg9_);
                case 8u:
                    check_dependency(arg8_);
                case 7u:
                    check_dependency(arg7_);
                case 6u:
                    check_dependency(arg6_);
                case 5u:
                    check_dependency(arg5_);
                case 4u:
                    check_dependency(arg4_);
                case 3u:
                    check_dependency(arg3_);
                case 2u:
                    check_dependency(arg2_);
                case 1u:
                    check_dependency(arg1_);
                default:
                    ;
            }
        }

        // Copies are not allowed.
        TaskFn(const TaskFn_&);
        TaskFn_ operator=(TaskFn_&);

    public:

        TaskFn(const futureT& result, functionT func, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(), arg2_(),
            arg3_(), arg4_(), arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 0u);
            check_dependencies();
        }

        template <typename a1T>
        TaskFn(const futureT& result, functionT func, const a1T& a1,
                const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(),
            arg3_(), arg4_(), arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 1u);
            check_dependencies();
        }

        template <typename a1T, typename a2T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const TaskAttributes& attr = TaskAttributes()) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(), arg4_(), arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 2u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(), arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 3u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 4u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(a5), arg6_(), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 5u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(a5), arg6_(a6), arg7_(), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 6u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T, typename a7T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const a7T& a7, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(a5), arg6_(a6), arg7_(a7), arg8_(), arg9_()
        {
            MADNESS_ASSERT(arity == 7u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T, typename a7T, typename a8T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6, const a7T& a7,
                const a8T& a8, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(a5), arg6_(a6), arg7_(a7), arg8_(a8), arg9_()
        {
            MADNESS_ASSERT(arity == 8u);
            check_dependencies();
        }

        template <typename a1T, typename a2T, typename a3T, typename a4T, typename a5T,
                typename a6T, typename a7T, typename a8T, typename a9T>
        TaskFn(const futureT& result, functionT func, const a1T& a1, const a2T& a2,
                const a3T& a3, const a4T& a4, const a5T& a5, const a6T& a6,
                const a7T& a7, const a8T& a8, const a9T& a9, const TaskAttributes& attr) :
            TaskInterface(attr), result_(result), func_(func), arg1_(a1), arg2_(a2),
            arg3_(a3), arg4_(a4), arg5_(a5), arg6_(a6), arg7_(a7), arg8_(a8), arg9_(a9)
        {
            MADNESS_ASSERT(arity == 9u);
            check_dependencies();
        }

        virtual ~TaskFn() { }

        void run(World&) {
            runner<arity, resultT> ();
        }

        const futureT& result() const { return result_; }

    protected:
        // Intel does not like using so explicitly instantiate it.
        void run(const TaskThreadEnv& env) {
            TaskInterface::run(env);
        }
    }; // class TaskFn

//} // namespace madness


#endif // MADNESS_WORLD_TASKFN_H__INCLUDED
