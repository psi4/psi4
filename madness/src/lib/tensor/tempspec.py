
typelist = ["int","long","double","float","double_complex","float_complex"]

complex_typelist = ["double_complex","float_complex"]

gemm_typelist = ["double","float","double_complex","float_complex"]

transform_typelist = [["double", "double"],
                      ["float", "float"],
                      ["double_complex", "double_complex"],
                      ["float_complex", "float_complex"],
                      ["double_complex", "double"],
                      ["double"        ,"double_complex"],
                      ["float_complex", "float"],
                      ["float"         ,"float_complex"]]

f = open("tensor_spec.h","w")
for t in typelist:
    f.write("\n// Instantiations for %s\n" % t)
    f.write("template class Tensor<%s>;\n" % t)
    f.write("template class SliceTensor<%s>;\n" % t)
    f.write("template std::ostream& operator << (std::ostream& s, const Tensor<%s>& t);\n" % t)
    f.write("template Tensor<%s> copy(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor<%s> outer(const Tensor<%s>& left, const Tensor<%s>& right);\n" % (t,t,t))
    f.write("template Tensor< Tensor<%s>::scalar_type > abs(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor<%s> transpose(const Tensor<%s>& t);\n" % (t,t))

f.write("\n// Instantiations for inner, inner_result and transform \n")
for t,q in transform_typelist:
    f.write("template void inner_result(const Tensor<%s>& left, const Tensor<%s>& right,\n" % (t,q))
    f.write("                           long k0, long k1, Tensor< TensorResultType<%s,%s>::type >& result);\n" % (t,q))
    f.write("template Tensor<TensorResultType<%s,%s>::type> inner(const Tensor<%s>& left, const Tensor<%s>& right,\n" % (t,q,t,q))
    f.write("                         long k0, long k1);\n")
    f.write("template Tensor<TensorResultType<%s,%s>::type> transform(const Tensor<%s>& t, const Tensor<%s>& c);\n" % (t,q,t,q))
    f.write("template Tensor<TensorResultType<%s,%s>::type> general_transform(const Tensor<%s>& t, const Tensor<%s> c[]);\n" % (t,q,t,q))
    f.write("template Tensor<TensorResultType<%s,%s>::type>& fast_transform(const Tensor<%s>& t, const Tensor<%s>& c, Tensor< TensorResultType<%s,%s>::type >& result, Tensor< TensorResultType<%s,%s>::type >& work);\n" % (t,q,t,q,t,q,t,q))


f.write("\n// Instantiations only for complex types\n")
for t in complex_typelist:
    f.write("\n// Instantiations for %s" % t)
    f.write("template Tensor< Tensor<%s>::scalar_type > arg(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor< Tensor<%s>::scalar_type > real(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor< Tensor<%s>::scalar_type > imag(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor<%s> conj(const Tensor<%s>& t);\n" % (t,t))
    f.write("template Tensor<%s> conj_transpose(const Tensor<%s>& t);\n" % (t,t))

#for t in typelist:
#    for q in typelist:
#        f.write("template TENSOR_RESULT_TYPE(%s,%s) Tensor<%s>::trace_conj(const Tensor<%s>& q);\n" % (t,q,t,q))

f.close()
    
## f.write("\n\nPUT THESE AT THE BOTTOM OF mxm.cc\n"
## for t in typelist:
##     f.write("template void mTxm<%s>(long dimi, long dimj, long dimk, %s * c, const %s *a, const %s *b);" % (t,t,t,t)
##     f.write("template void mxmT<%s>(long dimi, long dimj, long dimk, %s * c, const %s *a, const %s *b);" % (t,t,t,t)
##     f.write("template void mTxmT<%s>(long dimi, long dimj, long dimk, %s * c, const %s *a, const %s *b);" % (t,t,t,t)
##     f.write("template void mxm<%s>(long dimi, long dimj, long dimk, %s * c, const %s *a, const %s *b);" % (t,t,t,t)
    

##f.write("\n\nPUT THESE AT THE BOTTOM OF tensoriter.cc\n"
f = open("tensoriter_spec.h","w")

a = []
for t in typelist:
    for q in typelist:
        for r in typelist:
            a.append("template class TensorIterator<%s,%s,%s>;\n" % (t,q,r))
a.sort()

prev = a[0]
f.write(prev)
for cur in a[1:]:
    if cur != prev:
        f.write(cur)
        prev = cur
        
f.close()
