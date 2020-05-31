# Restricted HF reference, singlet states
add_subdirectory(h2o-adc1)
add_subdirectory(h2o-adc2)
add_subdirectory(h2o-adc2x)
add_subdirectory(h2o-adc3)
add_subdirectory(h2o-cvs-adc1)
add_subdirectory(h2o-cvs-adc2)
add_subdirectory(h2o-cvs-adc2x)
add_subdirectory(h2o-cvs-adc3)
add_subdirectory(h2o-fc-adc2)
add_subdirectory(h2o-fv-adc2x)

# Restricted HF reference, other spin
add_subdirectory(h2o-cvs-adc2x-triplets)
add_subdirectory(h2o-adc2-any)

# Unrestricted HF reference
add_subdirectory(cn-adc1)
add_subdirectory(cn-adc2)
add_subdirectory(hf-adc2-spin-flip)

# Specific properties
add_subdirectory(methyloxirane-rotation)

if(ENABLE_cppe)
    add_subdirectory(formaldehyde-pe-adc2)
endif()
