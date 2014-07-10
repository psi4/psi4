#ifndef _ABSTRACT_OPS_H_
#define _ABSTRACT_OPS_H_

/* abstract operations, 'regular' (reg) and 'complex' (cpl) */

#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

static double __elem_op_var;
static double __elem_op_var2;

/* assignment e.g. a = b */
#define assign_reg(a,b) (a) = (b)
#define assign_cpl(a,b) (a).real = (b).real; \
                        (a).imag = (b).imag

/* assignment of zero e.g. a = 0 */
#define assign_zero_reg(a) (a) = 0
#define assign_zero_cpl(a) (a).real = 0; \
                           (a).imag = 0

/* assignment of a sum of two values e.g. a = b + c */
#define assign_add_reg(a,b,c) (a) = ((b) + (c))
#define assign_add_cpl(a,b,c) (a).real = ((b).real + (c).real); \
                              (a).imag = ((b).imag + (c).imag)

/* assignment of a product of two values e.g. a = b * c */
#define assign_mul_reg(a,b,c) (a) = ((b) * (c))
#define assign_mul_cpl(a,b,c) (a).real = ((b).real*(c).real-(b).imag*(c).imag);\
                              (a).imag = ((b).real*(c).imag+(b).imag*(c).real)

/* assignment of a product of two valus e.g. a = b * c */
#define assign_mul_constant_reg(a,b,c) (a) = ((b) * (c))
#define assign_mul_constant_cpl(a,b,c) (a).real = ((b) * (c).real); \
                                       (a).imag = ((b) * (c).imag)

/* assignment of a product of two valus e.g. a = b * c */
#define assign_mul_reg(a,b,c) (a) = ((b) * (c))
#define assign_mul_cpl(a,b,c) (a).real = ((b).real*(c).real-(b).imag*(c).imag); \
                              (a).imag = ((b).real*(c).imag+(b).imag*(c).real)

/* assignment of a quotient of two valus e.g. a = b / c */
#if 0
#define assign_div_reg(a,b,c) (a) = ((b) / (c))
#define assign_div_cpl(a,b,c) (a).real = (((b).real*(c).real+(b).imag*(c).imag) \
                                         /((c).real*(c).real+(c).imag*(c).imag)); \
                              (a).imag = (((b).imag*(c).real-(b).real*(c).imag) \
                                         /((c).real*(c).real+(c).imag*(c).imag))
#else
#define assign_div_reg(a,b,c) (a) = ((b) / (c))
#define assign_div_cpl(a,b,c) __elem_op_var = ((c).real*(c).real+(c).imag*(c).imag); \
                              (a).real = (((b).real*(c).real+(b).imag*(c).imag) \
                                         /__elem_op_var; \
                              (a).imag = (((b).imag*(c).real-(b).real*(c).imag) \
                                         /__elem_op_var
#endif

/* in-place assignment of a sum e.g. a = a + b written a += b */
#define add_assign_reg(a,b) (a) += (b)
#define add_assign_cpl(a,b) (a).real += (b).real; \
                            (a).imag += (b).imag

/* not equal to zero e.g. a != 0 */
#define neq_zero_reg(a) (0 != (a))
#define neq_zero_cpl(a) (0 != (a).real || 0 != (a).imag)

/* equal to zero e.g. a == 0 */
#define eq_zero_reg(a) (0 == (a))
#define eq_zero_cpl(a) (0 == (a).real && 0 == (a).imag)

/* equality e.g. a == b */
#define eq_reg(a,b) ((a) == (b))
#define eq_cpl(a,b) ((a).real == (b).real && (a).imag == (b).imag)

/* absolute value */
#define abs_reg(a) ((a) < 0 ? -(a) : (a))
#define abs_cpl(a) (a) = sqrt((a).real*(a).real+(a).imag*(a).imag)

/* assignment of a maximum of two values e.g. if(b > c) a = b else a = c */
#define assign_max_reg(a,b,c) (a) = (b) > (c) ? (b) : (c)
#define assign_max_cpl(a,b,c) __elem_op_var = ((b).real*(b).real+(b).imag*(b).imag); \
                              __elem_op_var2 = ((c).real*(c).real+(c).imag*(c).imag); \
                              (a).real = __elem_op_var > __elem_op_var2 \
                                       ? (b).real : (c).real; \
                              (a).imag = __elem_op_var > __elem_op_var2 \
                                       ? (b).imag : (c).imag

/* assignment of a miniimum of two values e.g. if(b > c) a = b else a = c */
#define assign_min_reg(a,b,c) (a) = (b) < (c) ? (b) : (c)
#define assign_min_cpl(a,b,c)  __elem_op_var = ((b).real*(b).real+(b).imag*(b).imag); \
                               __elem_op_var2 = ((c).real*(c).real+(c).imag*(c).imag); \
                              (a).real = __elem_op_var < __elem_op_var2 \
                                       ? (b).real : (c).real; \
                              (a).imag = __elem_op_var < __elem_op_var2 \
                                       ? (b).imag : (c).imag

/* assignment of an absolute value e.g. a = |b| */
#define assign_abs_reg(a,b) (a) = abs_reg(b)
/* Note: absolute value of a complex number is usually sqrt(x*x + y*y) but this
 * can lead to overflows and/or underflows. Instead, we use the well-known
 * hypot solution:
 *       double hypot(double x,double y)
 *       {
 *           double t;
 *           x = abs(x);
 *           y = abs(y);
 *           t = min(x,y);
 *           x = max(x,y);
 *           y = t;
 *           return x*sqrt(1+(y/x)*(y/x));
 *       }
 */
#if HAVE_HYPOT
#define assign_abs_cpl(a,b) (a).real = hypot((b).real, (b).imag); \
                            (a).imag = 0.0
#else
#define assign_abs_cpl(a,b)                                         \
if (abs_reg((b).real) >= abs_reg((b).imag)) {                       \
    (a).real = abs_reg((b).real) *                                  \
               sqrt(1.0 + ((b).imag/(b).real)*((b).imag/(b).real)); \
} else {                                                            \
    (a).real = abs_reg((b).imag) *                                  \
               sqrt(1.0 + ((b).real/(b).imag)*((b).real/(b).imag)); \
}                                                                   \
(a).imag = 0.0
#endif


/* assignment of a random value */
#define sign (1.0 * rand() / RAND_MAX > 1.0 ? 1.0 : -1.0)
#define assign_rand_reg(a,val) (a) = 1.0 * (val) * rand() / RAND_MAX * sign
#define assign_rand_cpl(a,val) assign_rand_reg((a).real, (val).real); \
                               assign_rand_reg((a).imag, (val).imag)

#endif /* _ABSTRACT_OPS_H_ */
