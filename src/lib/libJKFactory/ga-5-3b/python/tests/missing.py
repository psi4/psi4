"""reports how much of numpy as been overridden or imported by gain"""
import numpy as np
from ga4py import gain
import inspect

np_function_count = 0
ov_function_count = 0
ov_function_miss = []
np_class_count = 0
ov_class_count = 0
ov_class_miss = []
np_ufunc_count = 0
ov_ufunc_count = 0
ov_ufunc_miss = []
for name in dir(np):
    np_obj = getattr(np, name)
    override = False
    if hasattr(gain, name):
        #if not me: print "gain override exists for: %s" % name
        override = True
    #else:
    #    setattr(gain, name, getattr(np, name))
    if type(np_obj) is type(np.add):
        np_ufunc_count += 1
        if override:
            ov_ufunc_count += 1
        else:
            ov_ufunc_miss.append(name)
    elif inspect.isfunction(np_obj):
        np_function_count += 1
        if override:
            ov_function_count += 1
        else:
            ov_function_miss.append(name)
    elif inspect.isclass(np_obj):
        np_class_count += 1
        if override:
            ov_class_count += 1
        else:
            ov_class_miss.append(name)
print "%d/%d numpy functions overridden by gain" % (
        ov_function_count,np_function_count)
if ov_function_count != np_function_count:
    print "missing functions"
    print ov_function_miss
print "%d/%d numpy classes overridden by gain" % (
        ov_class_count,np_class_count)
if ov_class_count != np_class_count:
    print "missing classes"
    print ov_class_miss
print "%d/%d numpy ufuncs overridden by gain" % (
        ov_ufunc_count,np_ufunc_count)
if ov_ufunc_count != np_ufunc_count:
    print "missing ufuncs"
    print ov_ufunc_miss
