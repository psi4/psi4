namespace psi{ namespace libmatrix {

    namespace detail {
        template <typename matrix_type>
        struct create_matrix {        
            static matrix_type create(const std::string& name, const Dimension& m, const Dimension& n) {
                return matrix_type (name, m, n);
            }
        };
        
        // In the case that the generic template isn't sufficient to create the matrix
        // write a specialization of it, e.g:
        //
        // template <>
        // struct create_matrix<libmints_matrix_wrapper> {
        //     typedef libmints_matrix_wrapper matrix_type;
        //     matrix_type create(const std::string& name, const Dimension& m, const Dimension& n) {
        //         return matrix_type(name, m, n);
        //     }
        // }

        // If an implementation is not available due to configuration, derive from this class
        // so that compiler errors will be generated if the developer attempts to use them.
        class is_not_available
        {
        protected:
            is_not_available() {}
            ~is_not_available() {}
        private:  // emphasize the following members are private
            is_not_available( const is_not_available& );
            const is_not_available& operator=( const is_not_available& );
        };
    }

    // Macro for defining unavailable matrix wrappers
    #define UNAVAILABLE_MATRIX(type) struct type : public detail::is_not_available { \
        void print() const; \
        void fill(double); \
        void add(const type&); \
    private: \
        type(); \
        type(const std::string&, const Dimension&, const Dimension&); \
    };

}} // end namespaces

