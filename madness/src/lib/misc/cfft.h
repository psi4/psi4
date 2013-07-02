//   fft.h - declaration of class
//   of fast Fourier transform - FFT
//
//   The code is property of LIBROW
//   You can use it on your own
//   When utilizing credit LIBROW site

#ifndef MADNESS_MISC_CFFT_H__INCLUDED
#define MADNESS_MISC_CFFT_H__INCLUDED

//   Include complex numbers header
#include <complex>

typedef std::complex<double> double_complex;

class CFFT
{
public:
    //   FORWARD FOURIER TRANSFORM
    //     Input  - input data
    //     Output - transform result
    //     N      - length of both input data and result
    static bool Forward(const double_complex *const Input, double_complex *const Output, const unsigned int N);
    
    //   FORWARD FOURIER TRANSFORM, INPLACE VERSION
    //     Data - both input data and output
    //     N    - length of input data
    static bool Forward(double_complex *const Data, const unsigned int N);
    
    //   INVERSE FOURIER TRANSFORM
    //     Input  - input data
    //     Output - transform result
    //     N      - length of both input data and result
    //     Scale  - if to scale result
    static bool Inverse(const double_complex *const Input, double_complex *const Output, const unsigned int N, const bool Scale = true);
    
    //   INVERSE FOURIER TRANSFORM, INPLACE VERSION
    //     Data  - both input data and output
    //     N     - length of both input data and result
    //     Scale - if to scale result
    static bool Inverse(double_complex *const Data, const unsigned int N, const bool Scale = true);
    
protected:
    //   Rearrange function and its inplace version
    static void Rearrange(const double_complex *const Input, double_complex *const Output, const unsigned int N);
    static void Rearrange(double_complex *const Data, const unsigned int N);
    
    //   FFT implementation
    static void Perform(double_complex *const Data, const unsigned int N, const bool Inverse = false);
    
    //   Scaling of inverse FFT result
    static void Scale(double_complex *const Data, const unsigned int N);
};

#endif // MADNESS_MISC_CFFT_H__INCLUDED
