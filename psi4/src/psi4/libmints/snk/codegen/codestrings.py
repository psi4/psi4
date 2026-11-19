#
#   For strings used in printing code
#

def get_codestrings(hrr, need_coeff_ratios, same_center):
    args_str_T = (
            "        T* integrals,              //!< overwritten, cartesian integrals for this grid batch\n" +
            ("        T* base_integrals,         //!< overwritten, base [LA+LB|0] integrals for HRR\n" if hrr else "") +
            "        const T* v_scaled_boys,    //!< Obara-Saika S-type seed integrals (primitive index slow)\n" +
            "        const T* v_p,              //!< (a + b) where a/b is the exponent of multishell A/B respectively\n" +
            "        const T* v_PGx,            //!< (Px - Gx); P is the GPT combined center; G is a gridpoint (primitive index slow)\n" +
            "        const T* v_PGy,            //!< (Py - Gy); P is the GPT combined center; G is a gridpoint (primitive index slow)\n" +
            "        const T* v_PGz,            //!< (Pz - Gz); P is the GPT combined center; G is a gridpoint (primitive index slow)"
        )
    args_str_f = (
            "        float* integrals, \n" +
            ("        float* base_integrals, \n" if hrr else "") +
            "        const float* v_scaled_boys, \n" +
            "        const float* v_p, \n" +
            "        const float* v_PGx, \n" +
            "        const float* v_PGy, \n" +
            "        const float* v_PGz, "
        )
    args_str_d = (
            "        double* integrals, \n" +
            ("        double* base_integrals, \n" if hrr else "") +
            "        const double* v_scaled_boys, \n" +
            "        const double* v_p, \n" +
            "        const double* v_PGx, \n" +
            "        const double* v_PGy, \n" +
            "        const double* v_PGz, "
        )
    if not same_center:
        args_str_T += (
                "\n        const T* v_PAx,            //!< (Px - Ax); P is the GPT combined center; A is the center of multishell A"
                "\n        const T* v_PAy,            //!< (Px - Ax); P is the GPT combined center; A is the center of multishell A"
                "\n        const T* v_PAz,            //!< (Px - Ax); P is the GPT combined center; A is the center of multishell A"
                "\n        T ABx,                     //!< (Ax - Bx); A/B is the center of multishell A/B, respectively"
                "\n        T ABy,                     //!< (Ax - Bx); A/B is the center of multishell A/B, respectively"
                "\n        T ABz,                     //!< (Ax - Bx); A/B is the center of multishell A/B, respectively"
            )
        args_str_f += (
                "\n        const float* v_PAx, "
                "\n        const float* v_PAy, "
                "\n        const float* v_PAz, "
                "\n        float ABx, "
                "\n        float ABy, "
                "\n        float ABz, "
            )
        args_str_d += (
                "\n        const double* v_PAx, "
                "\n        const double* v_PAy, "
                "\n        const double* v_PAz, "
                "\n        double ABx, "
                "\n        double ABy, "
                "\n        double ABz, "
            )
    if need_coeff_ratios:
        args_str_T += (
                "\n        const T* v_ratios,         //!< ratios of contraction coefficients to max L coefficient (primitive index slow)"
            )
        args_str_f += (
                "\n        const float* v_ratios, "
            )
        args_str_d += (
                "\n        const double* v_ratios, "
            )

    header_str = args_str_T
    header_str += (
            "\n        int nprim,                 //!< number of primitive shellpairs in this contracted shellpair"
            "\n        int gb_size);              //!< number of gridpoints in this grid batch"
        )
    args_str_T += (
            "\n        int nprim,                 //!< number of primitive shellpairs in this contracted shellpair"
            "\n        int gb_size)               //!< number of gridpoints in this grid batch"
        )
    args_str_f += (
            "\n        int nprim, "
            "\n        int gb_size)"
        )
    args_str_d += (
            "\n        int nprim, "
            "\n        int gb_size)"
        )
    return args_str_T, args_str_f, args_str_d, header_str
