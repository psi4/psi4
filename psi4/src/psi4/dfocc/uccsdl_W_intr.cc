/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "dfocc.h"

namespace psi{
namespace dfoccwave {

// ccsdl_Z_ints
void DFOCC::uccsdl_ZMBEJ_AAAA()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // AAAA Block
    // Z_MBEJ =  <MB||EJ> - \sum_{NF} t_NJ^BF <MN||EF> + \sum_{nf} t_Jn^Bf <Mn|Ef>     (91) 
    // Z_MBEJ =  <MB||EJ>
    // Z(ME,JB) = (ME|JB) - <ME|JB>
    Z = std::make_shared<Tensor2d>("Z (ME|JB)", naoccA, navirA, naoccA, navirA);
    Z->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (IJ|AB)", naoccA, naoccA, navirA, navirA);
    L->gemm(true, false, bQijA, bQabA, 1.0, 0.0);
    Z->sort(1324, L, -1.0, 1.0);
    L.reset();
    // Z_MBEJ -= \sum_{NF} t_NJ^BF <MN||EF>
    // <MN||EF> = <MN|EF> - <MN|FE> = (ME|FN) - <ME|FN>
    // t <JN|BF> = t (NF|JB) (sort: 2413)
    J = std::make_shared<Tensor2d>("J (ME|NF)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (ME|NF)", naoccA, navirA, naoccA, navirA);
    K->sort(1432, J, -1.0, 0.0);
    K->axpy(J,1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T_ <NF|JB>", naoccA, navirA, naoccA, navirA);
    T->sort(1423, T2, 1.0, 0.0);
    T2.reset();
    Z->gemm(false, false, K, T, -1.0, 1.0);
    K.reset();
    T.reset();
    // Z_MBEJ += \sum_{nf} t_Jn^Bf <Mn|Ef> 
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <nf|JB>", naoccB, navirB, naoccA, navirA);
    T->sort(2413, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (ME|nf)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    Z->gemm(false, false, J, T, 1.0, 1.0);
    Z->write(psio_, PSIF_DFOCC_AMPS);
//Z->print();
    Z.reset();
    J.reset();
    T.reset();
} //ccsdl_ZMBEJ_AAAA()

void DFOCC::uccsdl_Zmbej_BBBB()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // BBBB Block
    // Z_mbej =  <mb||ej> - \sum_{nf} t_nj^bf <mn||ef> + \sum_{NF} t_Nj^Fb <Nm|Fe>     (92)
    // Z_mbej =  <mb||ej>
    // Z(me,jb) = (me|jb) - <me|jb>
    Z = std::make_shared<Tensor2d>("Z (me|jb)", naoccB, navirB, naoccB, navirB);
    Z->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    L = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ij|ab)", naoccB, naoccB, navirB, navirB);
    L->gemm(true, false, bQijB, bQabB, 1.0, 0.0);
    Z->sort(1324, L, -1.0, 1.0);
    L.reset();
    // Z_mbej -=  \sum_{nf} t_nj^bf <mn||ef> 
    J = std::make_shared<Tensor2d>("J (me|nf)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (me|nf)", naoccB, navirB, naoccB, navirB);
    K->sort(1432, J, -1.0, 0.0);
    K->axpy(J,1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T_ <nf|jb>", naoccB, navirB, naoccB, navirB);
    T->sort(1423, T2, 1.0, 0.0);
    T2.reset();
    Z->gemm(false, false, K, T, -1.0, 1.0);
    K.reset();
    T.reset();
    // Z_mbej +=  \sum_{NF} t_Nj^Fb <Nm|Fe>
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <NF|jb>", naoccA, navirA, naoccB, navirB);
    T->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (me|NF)", naoccB, navirB, naoccA, navirA);
    J->gemm(true, false, bQiaB, bQiaA, 1.0, 0.0);
    Z->gemm(false, false, J, T, 1.0, 1.0);
    Z->write(psio_, PSIF_DFOCC_AMPS);
//Z->print();
    Z.reset();
    J.reset();
    T.reset();
} // ccsdl_Zmbej_BBBB()

void DFOCC::uccsdl_ZMbEj_ABAB()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // ABAB Block
    // Z_MbEj =  <Mb|Ej> - \sum_{nf} t_nj^bf <Mn|Ef> + \sum_{NF} t_Nj^Fb <MN||EF>     (93)
    // Z_MbEj =  <Mb|Ej> 
    // Z(ME,jb) = (ME|jb)
    Z = std::make_shared<Tensor2d>("Z (ME|jb)", naoccA, navirA, naoccB, navirB);
    Z->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    // Z_MbEj -= \sum_{nf} t_nj^bf <Mn|Ef>
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T_ <nf|jb>", naoccB, navirB, naoccB, navirB);
    T->sort(1423, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (ME|nf)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    Z->gemm(false, false, J, T, -1.0, 1.0);
    J.reset();
    T.reset();
    // Z_MbEj += \sum_{NF} t_Nj^Fb <MN||EF>
    J = std::make_shared<Tensor2d>("J (ME|NF)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (ME|NF)", naoccA, navirA, naoccA, navirA);
    K->sort(1432, J, -1.0, 0.0);
    K->axpy(J,1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <NF|jb>", naoccA, navirA, naoccB, navirB);
    T->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    Z->gemm(false, false, K, T, 1.0, 1.0);
//Z->print();
    Z->write(psio_, PSIF_DFOCC_AMPS);
    Z.reset();
    K.reset();
    T.reset();
} // ccsdl_ZMbEj_ABAB()

void DFOCC::uccsdl_ZmBeJ_BABA()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // BABA Block
    // Z_mBeJ =  <Bm|Je> - \sum_{NF} t_NJ^BF <Nm|Fe> + \sum_{nf} t_Jn^Bf <mn||ef>      (94)  
    // Z_mBeJ =  <Bm|Je>
    // <Bm|Je> = <mB|eJ> = (me|BJ) = (me|JB) 
    Z = std::make_shared<Tensor2d>("Z (me|JB)", naoccB, navirB, naoccA, navirA);
    Z->gemm(true, false, bQiaB, bQiaA, 1.0, 0.0);
    // Z_mBeJ -= \sum_{NF} t_NJ^BF <Nm|Fe>  
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T_ <NF|JB>", naoccA, navirA, naoccA, navirA);
    T->sort(1423, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (me|NF)", naoccB, navirB, naoccA, navirA);
    J->gemm(true, false, bQiaB, bQiaA, 1.0, 0.0);//UB
    Z->gemm(false, false, J, T, -1.0, 1.0);
    J.reset();
    T.reset();
    // Z_mBeJ += \sum_{nf} t_Jn^Bf <mn||ef>
    J = std::make_shared<Tensor2d>("J (me|nf)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (me|nf)", naoccB, navirB, naoccB, navirB);
    K->sort(1432, J, -1.0, 0.0);
    K->axpy(J,1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <nf|JB>", naoccB, navirB, naoccA, navirA);
    T->sort(2413, T2, 1.0, 0.0);
    T2.reset();
    Z->gemm(false, false, K, T, 1.0, 1.0);
    Z->write(psio_, PSIF_DFOCC_AMPS);
//Z->print();
    Z.reset();
    K.reset();
    T.reset();
} // ccsdl_ZmBeJ_BABA()

void DFOCC::uccsdl_ZMbeJ_ABBA()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // ABBA Block
    // Z_MbeJ =  - <Mb|Je> + \sum_{nF} t_Jn^Fb <Mn|Fe>     (95)
    // Z_MbeJ =  - <Mb|Je>
    // W(Me,Jb) = - <Me|Jb> = -(MJ|be) = Y (Me,Jb) (sort: 1423)
    J = std::make_shared<Tensor2d>("Int (MJ|be)", naoccA, naoccA, navirB, navirB);
    J->gemm(true, false, bQijA, bQabB, 1.0, 0.0);
    Z = std::make_shared<Tensor2d>("Z (Me|Jb)", naoccA, navirB, naoccA, navirB);
    Z->sort(1423, J, -1.0, 0.0);
    J.reset();
    // Z_MbeJ +=  \sum_{nF} t_Jn^Fb <Mn|Fe>   
    J = std::make_shared<Tensor2d>("J (MF|ne)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (Me|nF)", naoccA, navirB, naoccB, navirA);
    X->sort(1432, J, 1.0, 0.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <nF|Jb>", naoccB, navirA, naoccA, navirB);
    T->sort(2314, T2, 1.0, 0.0);
    T2.reset();
    Z->gemm(false, false, X, T, 1.0, 1.0);
    Z->write(psio_, PSIF_DFOCC_AMPS);
//Z->print();
    Z.reset();
    X.reset();
    T.reset(); 
}  // ccsdl_ZMbeJ_ABBA() 

void DFOCC::uccsdl_ZmBEj_BAAB()
{
    SharedTensor2d J, W, I, K, X, Y, T, Z, L;
    SharedTensor2d T2, Tau, T2new;
    // BAAB Block
    // Z_mBEj =  - <Bm|Ej> - \sum_{Nf} t_Nj^Bf <Nm|Ef>     (96) 
    // Z_mBEj =  - <Bm|Ej>
    // W(mE,jB) = - <mE|jB> = -(EB|mj)
    J = std::make_shared<Tensor2d>("Int (EB|mj)", navirA, navirA, naoccB, naoccB);
    J->gemm(true, false, bQabA, bQijB, 1.0, 0.0);
    Z = std::make_shared<Tensor2d>("Z (mE|jB)", naoccB, navirA, naoccB, navirA);
    Z->sort(3142, J, -1.0, 0.0);
    J.reset();
    // Z_mBEj -= \sum_{Nf} t_Nj^Bf <Nm|Ef> 
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <Nf|jB>", naoccA, navirB, naoccB, navirA);
    T->sort(1423, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (NE|mf)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (mE|Nf)", naoccB, navirA, naoccA, navirB);
    X->sort(3214, J, -1.0, 0.0);
    J.reset();
    Z->gemm(false, false, X, T, -1.0, 1.0);
    Z->write(psio_, PSIF_DFOCC_AMPS);
//Z->print();
    Z.reset();
    X.reset();
    T.reset();
}// ccsdl_ZmBEj_BAAB()

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DFOCC::uccsdl_WMBEJ_AAAA()
{
    SharedTensor2d W, X, K, T, Z, T1;

    // W_MBEJ = Z_MBEJ + \sum_{Q} t_JB^Qp b_ME^Q + \sum_{Q} t_BE^Q (t_JM^Q + b_JM^Q) - \sum_{Q} t_JM^Q b_BE^Q   (97)
    // W_MBEJ = Z_MBEJ  
    uccsdl_ZMBEJ_AAAA();
    Z = std::make_shared<Tensor2d>("Z (ME|JB)", naoccA, navirA, naoccA, navirA);
    Z->read(psio_, PSIF_DFOCC_AMPS);
    W = std::make_shared<Tensor2d>("WL (ME|JB)", naoccA, navirA, naoccA, navirA);
    W->copy(Z); 
    Z.reset();
    // W_MBEJ += \sum_{Q} t_JB^Qp b_ME^Q 
    T1 = std::make_shared<Tensor2d>("T1p (Q|IA)", nQ, naoccA, navirA);
    T1->read(psio_, PSIF_DFOCC_AMPS);
    W->gemm(true, false, bQiaA, T1, 1.0, 1.0);
    T1.reset();
    // W_MBEJ += \sum_{Q} t_BE^Q (t_JM^Q + b_JM^Q)
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->copy(bQijA);
    K->add(T);
    T.reset();
    T = std::make_shared<Tensor2d>("T1 (Q|AB)", nQ, navirA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (JM|BE)", naoccA, naoccA, navirA, navirA);
    X->gemm(true, false, K, T, 1.0, 0.0);
    T.reset();
    K.reset();
    // W_MBEJ -= \sum_{Q} t_JM^Q b_BE^Q 
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQabA, -1.0, 1.0);   
    T.reset(); 
    W->sort(2413, X, 1.0, 1.0);
    W->write(psio_, PSIF_DFOCC_AMPS);
//W->print();
    W.reset();
    X.reset();
} // end ccsdl_WMBEJ_AAAA 


void DFOCC::uccsdl_Wmbej_BBBB()
{
    SharedTensor2d W, X, K, T, Z, T1;

    // W_mbej = Z_mbej + \sum_{Q} t_jb^Qp b_me^Q + \sum_{Q} t_be^Q (t_jm^Q + b_jm^Q) - \sum_{Q} t_jm^Q b_be^Q   (98)
    // W_mbej = Z_mbej
    uccsdl_Zmbej_BBBB();
    Z = std::make_shared<Tensor2d>("Z (me|jb)", naoccB, navirB, naoccB, navirB);
    Z->read(psio_, PSIF_DFOCC_AMPS);
    W = std::make_shared<Tensor2d>("WL (me|jb)", naoccB, navirB, naoccB, navirB);
    W->copy(Z);
    Z.reset();    
    // W_mbej += \sum_{Q} t_jb^Qp b_me^Q 
    T1 = std::make_shared<Tensor2d>("T1p (Q|ia)", nQ, naoccB, navirB);
    T1->read(psio_, PSIF_DFOCC_AMPS);
    W->gemm(true, false, bQiaB, T1, 1.0, 1.0);
    T1.reset();
    // W_mbej += \sum_{Q} t_be^Q (t_jm^Q + b_jm^Q)
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB);
    K->copy(bQijB);
    K->add(T);
    T.reset();
    T = std::make_shared<Tensor2d>("T1 (Q|ab)", nQ, navirB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (jm|be)", naoccB, naoccB, navirB, navirB);
    X->gemm(true, false, K, T, 1.0, 0.0);
    T.reset();
    K.reset();
    // W_mbej -= \sum_{Q} t_jm^Q b_be^Q 
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQabB, -1.0, 1.0);
    T.reset();
    W->sort(2413, X, 1.0, 1.0);
    W->write(psio_, PSIF_DFOCC_AMPS);
//W->print();
    W.reset();
    X.reset();
} // end ccsdl_Wmbej_BBBB


void DFOCC::uccsdl_WMbEj_ABAB()
{
    SharedTensor2d  W, T1, Z;

    // W_MbEj = Z_MbEj + \sum_{Q} t_jb^Qp b_ME^Q   (99) 
    // W_MbEj = Z_MbEj
    uccsdl_ZMbEj_ABAB();
    Z = std::make_shared<Tensor2d>("Z (ME|jb)", naoccA, navirA, naoccB, navirB);
    Z->read(psio_, PSIF_DFOCC_AMPS);
    W = std::make_shared<Tensor2d>("WL (ME|jb)", naoccA, navirA, naoccB, navirB);
    W->copy(Z);   
    Z.reset();   
    // W_MbEj += \sum_{Q} t_jb^Qp b_ME^Q   
    T1 = std::make_shared<Tensor2d>("T1p (Q|ia)", nQ, naoccB, navirB);
    T1->read(psio_, PSIF_DFOCC_AMPS);
    W->gemm(true, false, bQiaA, T1, 1.0, 1.0);
    T1.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
//W->print();
    W.reset();
} // end ccsdl_WMbEj_ABAB

void DFOCC::uccsdl_WmBeJ_BABA()
{
    SharedTensor2d  W, T1, Z;
    // W_mBeJ = Z_mBeJ + \sum_{Q} t_JB^Qp b_me^Q  (100)
    // W_mBeJ = Z_mBeJ 
    uccsdl_ZmBeJ_BABA();
    Z = std::make_shared<Tensor2d>("Z (me|JB)", naoccB, navirB, naoccA, navirA);
    Z->read(psio_, PSIF_DFOCC_AMPS);
    W = std::make_shared<Tensor2d>("WL (me|JB)", naoccB, navirB, naoccA, navirA);
    W->copy(Z);
    Z.reset();
    // W_mBeJ += \sum_{Q} t_JB^Qp b_me^Q 
    T1 = std::make_shared<Tensor2d>("T1p (Q|IA)", nQ, naoccA, navirA);
    T1->read(psio_, PSIF_DFOCC_AMPS);
    W->gemm(true, false, bQiaB, T1, 1.0, 1.0);
    T1.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
//W->print();
    W.reset();
} // end ccsdl_WmBeJ_BABA

void DFOCC::uccsdl_WMbeJ_ABBA()
{
    SharedTensor2d W, Z, T, K, X;
    // W_MbeJ = Z_MbeJ + \sum_{Q} t_be^Q (t_JM^Q + b_JM^Q) - \sum_{Q} t_JM^Q b_be^Q   (101) 
    // W_MbeJ = Z_MbeJ 
    uccsdl_ZMbeJ_ABBA();
    Z = std::make_shared<Tensor2d>("Z (Me|Jb)", naoccA, navirB, naoccA, navirB);
    Z->read(psio_, PSIF_DFOCC_AMPS);
    W = std::make_shared<Tensor2d>("WL (Me|Jb)", naoccA, navirB, naoccA, navirB);
    W->copy(Z);
    Z.reset();
    // W_MbeJ += \sum_{Q} t_be^Q (t_JM^Q + b_JM^Q)
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->copy(bQijA);
    K->axpy(T, 1.0);
    T.reset();
    T = std::make_shared<Tensor2d>("T1 (Q|ab)", nQ, navirB, navirB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (JM|be)", naoccA, naoccA, navirB, navirB);
    X->gemm(true, false, K, T, 1.0, 0.0);
    T.reset();
    K.reset();
    // W_MbeJ -= \sum_{Q} t_JM^Q b_be^Q
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQabB, -1.0, 1.0);
    T.reset();
    W->sort(2413, X, 1.0, 1.0);
    X.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
//W->print();
    W.reset();

//// 
/*
    X = std::make_shared<Tensor2d>("WL (Me|Jb)", naoccA, navirB, naoccA, navirB);
    X->read(psio_, PSIF_DFOCC_AMPS);
    SharedTensor2d Y = std::make_shared<Tensor2d>("WL (ME|jb)", naoccA, navirA, naoccB, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    X->add(Y);
    Y.reset();
    X->print();
    X.reset();
*/

} // end ccsdl_WMbeJ_ABBA


void DFOCC::uccsdl_WmBEj_BAAB()
{
    SharedTensor2d W, Z, T, K, X;

    // W_mBEj = Z_mBEj + \sum_{Q} t_BE^Q (t_jm^Q + b_jm^Q) - \sum_{Q} t_jm^Q b_BE^Q   (102) 
    // W_mBEj = Z_mBEj 
    uccsdl_ZmBEj_BAAB();
    Z = std::make_shared<Tensor2d>("Z (mE|jB)", naoccB, navirA, naoccB, navirA);
    Z->read(psio_, PSIF_DFOCC_AMPS);
    W = std::make_shared<Tensor2d>("WL (mE|jB)", naoccB, navirA, naoccB, navirA);
    W->copy(Z);
    Z.reset();
    // W_mBEj += \sum_{Q} t_BE^Q (t_jm^Q + b_jm^Q) 
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB);
    K->copy(bQijB);
    K->axpy(T, 1.0);
    T.reset();
    T = std::make_shared<Tensor2d>("T1 (Q|AB)", nQ, navirA, navirA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (jm|BE)", naoccB, naoccB, navirA, navirA);
    X->gemm(true, false, K, T, 1.0, 0.0);
    T.reset();
    K.reset();
    // W_mBEj -= \sum_{Q} t_jm^Q b_BE^Q  
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQabA, -1.0, 1.0);
    T.reset();
    W->sort(2413, X, 1.0, 1.0);
    W->write(psio_, PSIF_DFOCC_AMPS);
//W->print();
    W.reset();
    X.reset();
} // end ccsdl_WmBEj_BAAB

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void DFOCC::uccsdl_WMNIE_AAAA()
{
    SharedTensor2d T, K, X, W;
    
    // W_MNIE = P_(MN) \sum_{Q} (t_IM^Q + b_IM^Q) * b_NE^Q   (103) 
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->copy(bQijA);
    K->add(T);
    T.reset();
    X = std::make_shared<Tensor2d>("X (IM|NE)", naoccA, naoccA, naoccA, navirA);
    X->gemm(true, false, K, bQiaA, 1.0, 0.0);
    K.reset();
    W = std::make_shared<Tensor2d>("WL (MN|IE)", naoccA, naoccA, naoccA, navirA);
    W->sort(2314, X, 1.0, 0.0);
    W->sort(3214, X, -1.0, 1.0);
    W->write(psio_, PSIF_DFOCC_AMPS);
//W->print();
    W.reset();
    X.reset();
} // end ccsdl_WMNIE_AAAA


void DFOCC::uccsdl_Wmnie_BBBB()
{
    SharedTensor2d T, K, X, W;

    // W_MNIE = P_(mn) \sum_{Q} (t_im^Q + b_im^Q) * b_ne^Q   (104) 
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB);
    K->copy(bQijB);
    K->add(T);
    T.reset();
    X = std::make_shared<Tensor2d>("X (im|ne)", naoccB, naoccB, naoccB, navirB);
    X->gemm(true, false, K, bQiaB, 1.0, 0.0);
    K.reset();
    W = std::make_shared<Tensor2d>("WL (mn|ie)", naoccB, naoccB, naoccB, navirB);
    W->sort(2314, X, 1.0, 0.0);
    W->sort(3214, X, -1.0, 1.0);    
    W->write(psio_, PSIF_DFOCC_AMPS);
//W->print();
    W.reset();
    X.reset();
} // end ccsdl_Wmnie_BBBB

void DFOCC::uccsdl_WMnIe_ABAB()
{
    SharedTensor2d T, K, X, W;

    // W_MnIe = \sum_{Q} (t_IM^Q + b_IM^Q) * b_ne^Q   (105) 
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|IJ)", nQ, naoccA, naoccA);
    K->copy(bQijA);
    K->add(T);
    T.reset();
    X = std::make_shared<Tensor2d>("X (IM|ne)", naoccA, naoccA, naoccB, navirB);
    X->gemm(true, false, K, bQiaB, 1.0, 0.0);
    K.reset();
    W = std::make_shared<Tensor2d>("WL (Mn|Ie)", naoccA, naoccB, naoccA, navirB);
    W->sort(2314, X, 1.0, 0.0);
    W->write(psio_, PSIF_DFOCC_AMPS);
//W->print();
/*
SharedTensor2d A = std::make_shared<Tensor2d>("X (Mn|Ie)", naoccA, naoccB, naoccA, navirB);
A->add(W);
A->sort(2134, W, -1.0, 1.0);
A->print();
A.reset();
*/
    W.reset();
    X.reset();


} // end ccsdl_WMnIe_ABAB

void DFOCC::uccsdl_WmNiE_BABA()
{
    SharedTensor2d T, K, X, W;

    // W_mNiE = \sum_{Q} (t_im^Q + b_im^Q) * b_NE^Q   (106) 
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    K = std::make_shared<Tensor2d>("DF_BASIS_CC B (Q|ij)", nQ, naoccB, naoccB);
    K->copy(bQijB);
    K->add(T);
    T.reset();
    X = std::make_shared<Tensor2d>("X (im|NE)", naoccB, naoccB, naoccA, navirA);
    X->gemm(true, false, K, bQiaA, 1.0, 0.0);
    K.reset();
    W = std::make_shared<Tensor2d>("WL (mN|iE)", naoccB, naoccA, naoccB, navirA);
    W->sort(2314, X, 1.0, 0.0);
    W->write(psio_, PSIF_DFOCC_AMPS);
//W->print();
    W.reset();
    X.reset();
} // end ccsdl_WmNiE_BABA

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void DFOCC::uccsdl_WMNIJ_AAAA()
{
    SharedTensor2d J, W, I, X, Y, T;
    SharedTensor2d T2, Tau, T2new;
    // W_MNIJ Alpha Block
    // W_MNIJ = <MN||IJ> + P_(MN) P_(IJ) \sum_(Q) t_IM^Q b_JN^Q + \sum_(EF) Tau(IJ,EF) * <MN|EF> 
    // W_MNIJ = <MN||IJ> 
    J = std::make_shared<Tensor2d>("J (IM|JN)", naoccA, naoccA, naoccA, naoccA);
    J->gemm(true, false, bQijA, bQijA, 1.0, 0.0);
    W = std::make_shared<Tensor2d>("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    W->sort(1324, J, 1.0, 0.0);
    W->sort(1342, J, -1.0, 1.0);
    J.reset();
    // W_MNIJ += P_(MN) P_(IJ) \sum_(Q) t_IM^Q b_JN^Q
    // X(IM,JN) = \sum(Q) t(Q,IM) b(Q,JN)
    // W_MNIJ += P_(MN) * P_(IJ) X(IM,JN)
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (IM|JN)", naoccA, naoccA, naoccA, naoccA);
    X->gemm(true, false, T, bQijA, 1.0, 0.0);
    T.reset();
    //W->P_ijab(X);
    W->sort(2413, X, 1.0, 1.0);
    W->sort(4213, X, -1.0, 1.0);
    W->sort(2431, X, -1.0, 1.0);
    W->sort(4231, X, 1.0, 1.0);
    X.reset();
    //W_MNIJ += \sum_(EF) Tau(IJ,EF) * <MN|EF> 
    J = std::make_shared<Tensor2d>("J (ME|NF)", naoccA, navirA, naoccA, navirA);
    J->gemm(true, false, bQiaA, bQiaA, 1.0, 0.0);
    I = std::make_shared<Tensor2d>("I <MN|EF>", naoccA, naoccA, navirA, navirA);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    uccsd_tau_amps(naoccA, naoccA, navirA, navirA, Tau, T2, t1A, t1A);    
    T2.reset();
    W->gemm(false, true, I, Tau, 1.0, 1.0);
    I.reset();
    Tau.reset();
    W->write(psio_, PSIF_DFOCC_AMPS);
    W.reset();
} // end ccsdl_WMNIJ()

void DFOCC::uccsdl_Wmnij_BBBB()
{
    SharedTensor2d J, W, I, X, Y, T;
    SharedTensor2d T2, Tau, T2new;

    // W_mnij Beta Block
    // W_mnij = <mn||ij> + P_(mn) P_(ij) \sum_(Q) t_im^Q b_jn^Q + \sum_(ef) Tau(ij,ef) * <mn|ef> 
    // W_mnij = <mn||ij> 
    J = std::make_shared<Tensor2d>("J (im|jn)", naoccB, naoccB, naoccB, naoccB);
    J->gemm(true, false, bQijB, bQijB, 1.0, 0.0);
    W = std::make_shared<Tensor2d>("W <mn|ij>", naoccB, naoccB, naoccB, naoccB);
    W->sort(1324, J, 1.0, 0.0);
    W->sort(1342, J, -1.0, 1.0);
    J.reset();
    // W_mnij += P_(mn) P_(ij) \sum_(Q) t_im^Q b_jn^Q
    // X(im,jn) = \sum(Q) t(Q,im) b(Q,jn)
    // W_mnij += P_(mn) * P_(ij) X(im,jn)
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (im|jn)", naoccB, naoccB, naoccB, naoccB);
    X->gemm(true, false, T, bQijB, 1.0, 0.0);
    T.reset();
    //W->P_ijab(X);
    W->sort(2413, X, 1.0, 1.0);
    W->sort(4213, X, -1.0, 1.0);
    W->sort(2431, X, -1.0, 1.0);
    W->sort(4231, X, 1.0, 1.0);
    X.reset();
    //W_mnij += \sum_(ef) Tau(ij,ef) * <mn|ef> 
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <ij|ab>", naoccB, naoccB, navirB, navirB);
    uccsd_tau_amps(naoccB, naoccB, navirB, navirB, Tau, T2, t1B, t1B);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (me|nf)", naoccB, navirB, naoccB, navirB);
    J->gemm(true, false, bQiaB, bQiaB, 1.0, 0.0);
    I = std::make_shared<Tensor2d>("I <mn|ef>", naoccB, naoccB, navirB, navirB);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();
    W->gemm(false, true, I, Tau, 1.0, 1.0);
    W->write(psio_, PSIF_DFOCC_AMPS);
    W.reset();
    I.reset();
    Tau.reset();
} // end ccsd_W_mnij()

void DFOCC::uccsdl_WMnIj_ABAB()
{
    SharedTensor2d J, W, I, X, Y, T;
    SharedTensor2d T2, Tau, T2new;

    // W_MnIj Alpha-Beta Block
    // W_MnIj = <Mn|Ij> + \sum_(Q) t_IM^Q b_jn^Q + \sum_(Ef) Tau(Ij,Ef) * <Mn|Ef> 
    // W_MnIj = <Mn|Ij> 
    W = std::make_shared<Tensor2d>("W <Mn|Ij>", naoccA, naoccB, naoccA, naoccB);
    J = std::make_shared<Tensor2d>("J (MI|nj)", naoccA, naoccA, naoccB, naoccB);
    J->gemm(true, false, bQijA, bQijB, 1.0, 0.0);
    W->sort(1324, J, 1.0, 0.0);
    J.reset();
    // X(IM,jn) = \sum(Q) t(Q,IM) b(Q,jn)
    // W_MnIj += X(IM,jn)
    X = std::make_shared<Tensor2d>("X (IM|jn)", naoccA, naoccA, naoccB, naoccB);
    T = std::make_shared<Tensor2d>("T1 (Q|IJ)", nQ, naoccA, naoccA);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, T, bQijB, 1.0, 0.0);
    T.reset();
    W->sort(2413, X, 1.0, 1.0);
    X.reset();
    // X(IM,jn) = \sum(Q) t(Q,jn) b(Q,IM)
    // W_MnIj += X(IM,jn)
    X = std::make_shared<Tensor2d>("X (IM|jn)", naoccA, naoccA, naoccB, naoccB);
    T = std::make_shared<Tensor2d>("T1 (Q|ij)", nQ, naoccB, naoccB);
    T->read(psio_, PSIF_DFOCC_AMPS);
    X->gemm(true, false, bQijA, T, 1.0, 0.0);
    T.reset();
    W->sort(2413, X, 1.0, 1.0);
    X.reset();
    //W_MnIj += \sum_(Ef) Tau(Ij,Ef) * <Mn|Ef> 
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (ME|nf)", naoccA, navirA, naoccB, navirB);
    J->gemm(true, false, bQiaA, bQiaB, 1.0, 0.0);
    I = std::make_shared<Tensor2d>("I <Mn|Ef>", naoccA, naoccB, navirA, navirB);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();
    W->gemm(false, true, I, Tau, 1.0, 1.0);
    W->write(psio_, PSIF_DFOCC_AMPS);
    W.reset();
    I.reset();
    Tau.reset();
} // end ccsd_W_MnIj()

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void DFOCC::uccsdl_WMBIJ_AAAA()
{
    SharedTensor2d G, W, WL, X, Y, Z, Z_, T, T2, Tau, J, K, I;

    // W_MBIJ = <MB||IJ> - \sum_{E} t_IJ^BE Ft_ME - \sum_{N} t_N^B W_MNIJ + \sum_{E,F} Tau_IJ^EF <MB|EF> + P_(IJ) \sum_{E} t_I^E Z_MBEJ 
    // + P_(IJ) \sum_{N,E} t(JN,BE) <MN||IE> + P_(IJ) \sum_{n,e} t(Jn,Be) <Mn|Ie>                                                         (109) 

    // W_MBIJ (1) = <MB||IJ>
    J = std::make_shared<Tensor2d>("J (MI|JB)", naoccA, naoccA, naoccA, navirA);
    J->gemm(true, false, bQijA, bQiaA, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (MI|JB)", naoccA, naoccA, naoccA, navirA);
    K->sort(1324, J, -1.0, 0.0);
    K->axpy(J, 1.0);
    J.reset();
    WL = std::make_shared<Tensor2d>("WL (MB|IJ)", naoccA, navirA, naoccA, naoccA);
    WL->sort(1423, K, 1.0, 0.0);
    K.reset();
    // W_MBIJ (2) -= \sum_{E} t_IJ^BE Ft_ME   (NOTE THAT: Ft_IA = F_IA)
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X <IJ|BM>", naoccA, naoccA, navirA, naoccA);
    X->contract(false, true, naoccA * naoccA * navirA, naoccA, navirA, T2, FiaA, 1.0, 0.0);
    T2.reset();
    WL->sort(4312, X, -1.0, 1.0);
    X.reset();
    // W_MBIJ (3) -= \sum_{N} t_N^B W_MNIJ
    uccsdl_WMNIJ_AAAA();
    Y = std::make_shared<Tensor2d>("W <MN|IJ>", naoccA, naoccA, naoccA, naoccA);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    W = std::make_shared<Tensor2d>("W <IJ|MN>", naoccA, naoccA, naoccA, naoccA);
    W->trans(Y);
    Y.reset();
    X = std::make_shared<Tensor2d>("X <IJ|MB>", naoccA, naoccA, naoccA, navirA);
    X->contract(false, false, naoccA * naoccA * naoccA, navirA, naoccA, W, t1A, 1.0, 0.0);
    W.reset();
    WL->sort(3412, X, -1.0, 1.0);
    X.reset();
    // W_MBIJ (4) += \sum_{E,F} Tau_IJ^EF <MB|EF>
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <IJ|AB>", naoccA, naoccA, navirA, navirA);
    uccsd_tau_amps(naoccA, naoccA, navirA, navirA, Tau, T2, t1A, t1A);
    T2.reset();
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (ME|BF)", naoccA, navirA, navirA, navirA);
    J->gemm(true, false, bQiaA, bQabA, 1.0, 0.0);
    I = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <MB|EF>", naoccA, navirA, navirA, navirA);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();
    WL->gemm(false, true, I, Tau, 1.0, 1.0);
    I.reset();
    Tau.reset();
    // W_MBIJ (5) += P_(IJ) \sum_{E} t_I^E Z_MBEJ 
    uccsdl_ZMBEJ_AAAA();
    Y = std::make_shared<Tensor2d>("Z (ME|JB)", naoccA, navirA, naoccA, navirA);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Z = std::make_shared<Tensor2d>("Z (JB|ME)", naoccA, navirA, naoccA, navirA);
    Z->trans(Y);
    Y.reset();
    X = std::make_shared<Tensor2d>("X (JB|MI)", naoccA, navirA, naoccA, naoccA);
    X->contract(false, true, naoccA * navirA * naoccA, naoccA, navirA, Z, t1A, 1.0, 0.0);
    Z.reset();
    WL->sort(3241, X, 1.0, 1.0);
    WL->sort(3214, X, -1.0, 1.0);
    X.reset();
    // W_MBIJ (6) += P_(IJ) \sum_{N,E} t(JN,BE) <MN||IE>  
    // <MN||IE> = <MN|IE> - <MN|EI>
    // K_ (MI,NE) = (MI|NE) - (ME|NI) -------->  (ME|NI) = (NI|ME) sort:3214  
    J = std::make_shared<Tensor2d>("J (MI|NE)", naoccA, naoccA, naoccA, navirA);
    J->gemm(true, false, bQijA, bQiaA, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (MI|NE)", naoccA, naoccA, naoccA, navirA);
    K->sort(3214, J, -1.0, 0.0);
    K->axpy(J, 1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T (JB|NE)", naoccA, navirA, naoccA, navirA);
    T->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    Y = std::make_shared<Tensor2d>("Y (MI|JB)", naoccA, naoccA, naoccA, navirA);
    Y->gemm(false, true, K, T, 1.0, 0.0);
    K.reset();
    T.reset();
    WL->sort(1423, Y, 1.0, 1.0);
    WL->sort(1432, Y, -1.0, 1.0);
    Y.reset();
    // W_MBIJ (7) += P_(IJ) \sum_{n,e} t(Jn,Be) <Mn|Ie> 
    J = std::make_shared<Tensor2d>("J (MI|ne)", naoccA, naoccA, naoccB, navirB);
    J->gemm(true, false, bQijA, bQiaB, 1.0, 0.0);
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T (ne|JB)", naoccB, navirB, naoccA, navirA);
    T->sort(2413, T2, 1.0, 0.0);
    T2.reset(); 
    X = std::make_shared<Tensor2d>("X (MI|JB)", naoccA, naoccA, naoccA, navirA);
    X->gemm(false, false, J, T, 1.0, 0.0);
    J.reset();
    T.reset();
    WL->sort(1423, X, 1.0, 1.0);
    WL->sort(1432, X, -1.0, 1.0);
    WL->write(psio_, PSIF_DFOCC_AMPS);
//WL->print();
    X.reset();
    WL.reset();

} // end ccsdl_WMBIJ_AAAA


void DFOCC::uccsdl_Wmbij_BBBB()
{    
    SharedTensor2d G, W, WL, X, Y, Z, Z_, T, T2, Tau, J, K, I;

    // W_mbij = <mb||ij> - \sum_{e} t_ij^be Ft_me - \sum_{n} t_n^b W_mnij + \sum_{e,f} Tau_ij^ef <mb|ef> + P_(ij) \sum_{e} t_i^e Z_mbej 
    // + P_(ij) \sum_{n,e} t(jn,be) <mn||ie> + P_(ij) \sum_{N,E} t(Nj,Eb) <Nm|Ei>                                                         (110) 

    // W_mbij (1) +=  <mb||ij> 
    J = std::make_shared<Tensor2d>("J (mi|jb)", naoccB, naoccB, naoccB, navirB);
    J->gemm(true, false, bQijB, bQiaB, 1.0, 0.0);
    //K = std::make_shared<Tensor2d>("K (mj|ib)", naoccB, naoccB, naoccB, navirB);
    //K->sort(1324, J, -1.0, 0.0);
    //K->axpy(J, 1.0);
    WL = std::make_shared<Tensor2d>("WL (mb|ij)", naoccB, navirB, naoccB, naoccB);
    WL->sort(1432, J, -1.0, 0.0);
    WL->sort(1423, J, 1.0, 1.0);
    J.reset();
    //K.reset();
    // W_mbij (2) -= \sum_{e} t_ij^be Ft_me
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X <ij|bm>", naoccB, naoccB, navirB, naoccB);
    X->contract(false, true, naoccB * naoccB * navirB, naoccB, navirB, T2, FiaB, 1.0, 0.0);
    T2.reset();
    WL->sort(4312, X, -1.0, 1.0);
    X.reset();
    // W_mbij (3) -= \sum_{n} t_n^b W_mnij
    uccsdl_Wmnij_BBBB();
    Y = std::make_shared<Tensor2d>("W <mn|ij>", naoccB, naoccB, naoccB, naoccB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    W = std::make_shared<Tensor2d>("W <ij|mn>", naoccB, naoccB, naoccB, naoccB);
    W->trans(Y);
    Y.reset();
    X = std::make_shared<Tensor2d>("X <ij|mb>", naoccB, naoccB, naoccB, navirB);
    X->contract(false, false, naoccB * naoccB * naoccB, navirB, naoccB, W, t1B, 1.0, 0.0);
    W.reset();
    WL->sort(3412, X, -1.0, 1.0);
    X.reset();
    // W_mbij (4) += \sum_{e,f} Tau_ij^ef <mb|ef>
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    Tau = std::make_shared<Tensor2d>("Tau <ij|ab>", naoccB, naoccB, navirB, navirB);
    uccsd_tau_amps(naoccB, naoccB, navirB, navirB, Tau, T2, t1B, t1B);
    T2.reset();
    J = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints (me|bf)", naoccB, navirB, navirB, navirB);
    J->gemm(true, false, bQiaB, bQabB, 1.0, 0.0);
    I = std::make_shared<Tensor2d>("DF_BASIS_CC MO Ints <mb|ef>", naoccB, navirB, navirB, navirB);
    I->sort(1324, J, 1.0, 0.0);
    J.reset();
    WL->gemm(false, true, I, Tau, 1.0, 1.0);
    I.reset();
    Tau.reset();
    // W_mbij (5) += P_(ij) \sum_{e} t_i^e Z_mbej
    uccsdl_Zmbej_BBBB();
    Y = std::make_shared<Tensor2d>("Z (me|jb)", naoccB, navirB, naoccB, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Z = std::make_shared<Tensor2d>("Z (jb|me)", naoccB, navirB, naoccB, navirB);
    Z->trans(Y);
    Y.reset();
    X = std::make_shared<Tensor2d>("X (jb|mi)", naoccB, navirB, naoccB, naoccB);
    X->contract(false, true, naoccB * navirB * naoccB, naoccB, navirB, Z, t1B, 1.0, 0.0);
    Z.reset();
    WL->sort(3241, X, 1.0, 1.0);
    WL->sort(3214, X, -1.0, 1.0);
    X.reset();
    // W_mbij (6) += P_(ij) \sum_{n,e} t(jn,be) <mn||ie>
    // <mn||ie> = <mn|ie> - <mn|ei>
    // K_ (mi,ne) = (mi|ne) - (me|ni) -------->  (me|ni) = (ni|me) sort:3214  
    J = std::make_shared<Tensor2d>("J (mi|ne)", naoccB, naoccB, naoccB, navirB);
    J->gemm(true, false, bQijB, bQiaB, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (mi|ne)", naoccB, naoccB, naoccB, navirB);
    K->sort(3214, J, -1.0, 0.0);
    K->axpy(J, 1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T (jb|ne)", naoccB, navirB, naoccB, navirB);
    T->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    Y = std::make_shared<Tensor2d>("Y (mi|jb)", naoccB, naoccB, naoccB, navirB);
    Y->gemm(false, true, K, T, 1.0, 0.0);
    K.reset();
    T.reset();
    WL->sort(1423, Y, 1.0, 1.0);
    WL->sort(1432, Y, -1.0, 1.0);
    Y.reset();
    // W_mbij (7) += P_(ij) \sum_{N,E} t(Nj,Eb) <Nm|Ei> 
    J = std::make_shared<Tensor2d>("J (mi|NE)", naoccB, naoccB, naoccA, navirA);
    J->gemm(true, false, bQijB, bQiaA, 1.0, 0.0);
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T (NE|jb)", naoccA, navirA, naoccB, navirB);
    T->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    X = std::make_shared<Tensor2d>("X (mi|jb)", naoccB, naoccB, naoccB, navirB);
    X->gemm(false, false, J, T, 1.0, 0.0);
    J.reset();
    T.reset();
    WL->sort(1423, X, 1.0, 1.0);
    WL->sort(1432, X, -1.0, 1.0);
    WL->write(psio_, PSIF_DFOCC_AMPS);
    X.reset();
    WL.reset();
//outfile->Printf("\tI am here.\n");
} // end ccsdl_Wmbij_BBBB 

void DFOCC::uccsdl_WMbIj_ABAB()
{
    SharedTensor2d G, W, WL, X, Y, Z, Z_, T2, Tau, J, K, T;

    // W_MbIj = <Mb|Ij> + \sum_{E} t_Ij^Eb Ft_ME - \sum_{n} t_n^b W_MnIj + \sum_{E,f} Tau_Ij^Ef <Mb|Ef> + \sum_{E} t_I^E Z_MbEj 
    // - \sum_{e} t_j^e Z_MbeI + \sum_{n,e} t_jn^be <Mn|Ie> + \sum_{N,E} t_Nj^Eb <MN||IE> - \sum_{n,E} t_In^Eb <Mn|Ej>                       (111)                                 
    // W_MbIj (1) += <Mb|Ij>
    J = std::make_shared<Tensor2d>("J (MI|jb)", naoccA, naoccA, naoccB, navirB);
    J->gemm(true, false, bQijA, bQiaB, 1.0, 0.0);
    WL = std::make_shared<Tensor2d>("WL (Mb|Ij)", naoccA, navirB, naoccA, naoccB);
    WL->sort(1423, J, 1.0, 0.0);
    J.reset();
    // W_MbIj (2) += \sum_{E} t_Ij^Eb Ft_ME 
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X <Ab|Ij>", navirA, navirB, naoccA, naoccB);
    X->trans(T2);
    T2.reset();
    WL->contract(false, false, naoccA, navirB * naoccA * naoccB, navirA, FiaA, X, 1.0, 1.0);
    X.reset();
    // W_MbIj (3) -= \sum_{n} t_n^b W_MnIj 
    uccsdl_WMnIj_ABAB();
    W = std::make_shared<Tensor2d>("W <Mn|Ij>", naoccA, naoccB, naoccA, naoccB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    Y = std::make_shared<Tensor2d>("W <Ij|Mn>", naoccA, naoccB, naoccA, naoccB);
    Y->trans(W);
    W.reset();
    X = std::make_shared<Tensor2d>("X <Ij|Mb>", naoccA, naoccB, naoccA, navirB);
    X->contract(false, false, naoccA * naoccB * naoccA, navirB, naoccB, Y, t1B, 1.0, 0.0);
    Y.reset();
    WL->sort(3412, X, -1.0, 1.0);
    X.reset();
    // W_MbIj (4) += \sum_{E,f} Tau_Ij^Ef <Mb|Ef>
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (ME|bf)", naoccA, navirA, navirB, navirB);
    J->gemm(true, false, bQiaA, bQabB, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G (Mb|Ef)", naoccA, navirB, navirA, navirB);
    G->sort(1324, J, 1.0, 0.0);
    J.reset();
    WL->gemm(false, true, G, Tau, 1.0, 1.0);
    G.reset();
    Tau.reset();
    // W_MbIj (5) += \sum_{E} t_I^E Z_MbEj 
    uccsdl_ZMbEj_ABAB();
    Y = std::make_shared<Tensor2d>("Z (ME|jb)", naoccA, navirA, naoccB, navirB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (jb|MI)", naoccB, navirB, naoccA, naoccA);
    Z = std::make_shared<Tensor2d>("Z (jb|ME)", naoccB, navirB, naoccA, navirA);
    Z->trans(Y);
    Y.reset();
    X->contract(false, true, naoccB * navirB * naoccA, naoccA, navirA, Z, t1A, 1.0, 0.0);
    Z.reset();
    WL->sort(3241, X, 1.0, 1.0);
    X.reset();
    // W_MbIj (6) -= \sum_{e} t_j^e Z_MbeI
    uccsdl_ZMbeJ_ABBA();
    Z = std::make_shared<Tensor2d>("Z (Me|Jb)", naoccA, navirB, naoccA, navirB);
    Z->read(psio_, PSIF_DFOCC_AMPS);
    Y = std::make_shared<Tensor2d>("Y (Ib|Me)", naoccA, navirB, naoccA, navirB);
    Y->trans(Z);
    Z.reset();
    X = std::make_shared<Tensor2d>("X (Ib|Mj)", naoccA, navirB, naoccA, naoccB);
    X->contract(false, true, naoccA * navirB * naoccA, naoccB, navirB, Y, t1B, 1.0, 0.0);
    Y.reset();
    WL->sort(3214, X, -1.0, 1.0);
    X.reset();
    // W_MbIj (7) += \sum_{n,e} t_jn^be <Mn|Ie>
    T2 = std::make_shared<Tensor2d>("T2 <ij|ab>", naoccB, naoccB, navirB, navirB);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T (jb|ne)", naoccB, navirB, naoccB, navirB); 
    T->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (MI|ne)", naoccA, naoccA, naoccB, navirB);
    J->gemm(true, false, bQijA, bQiaB, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (MI|jb)", naoccA, naoccA, naoccB, navirB);
    X->gemm(false, true, J, T, 1.0, 0.0);
    J.reset();
    T.reset();
    WL->sort(1423, X, 1.0, 1.0);
    X.reset(); 
    // W_MbIj (8) += \sum_{N,E} t_Nj^Eb <MN||IE> 
    // <MN||IE> = <MN|IE> - <MN|EI>
    // K_ (MI,NE) = (MI|NE) - (ME|NI) -------->  (ME|NI) = (NI|ME) sort:3214  
    J = std::make_shared<Tensor2d>("J (MI|NE)", naoccA, naoccA, naoccA, navirA);
    J->gemm(true, false, bQijA, bQiaA, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (MI|NE)", naoccA, naoccA, naoccA, navirA);
    K->sort(3214, J, -1.0, 0.0);
    K->axpy(J, 1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <NE|jb>", naoccA, navirA, naoccB, navirB);
    T->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    X= std::make_shared<Tensor2d>("X (MI|jb)", naoccA, naoccA, naoccB, navirB);
    X->gemm(false, false, K, T, 1.0, 0.0);
    K.reset();
    T.reset();
    WL->sort(1423, X, 1.0, 1.0);
    X.reset();
    // W_MbIj (9) -= \sum_{n,E} t_In^Eb <Mn|Ej> 
    J = std::make_shared<Tensor2d>("J (ME|nj)", naoccA, navirA, naoccB, naoccB);
    J->gemm(true, false, bQiaA, bQijB, 1.0, 0.0);
    Y = std::make_shared<Tensor2d>("Y (nE|Mj)", naoccB, navirA, naoccA, naoccB);
    Y->sort(3214, J, 1.0, 0.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <Ib|nE>", naoccA, navirB, naoccB, navirA);
    T->sort(1423, T2, 1.0, 0.0);
    T2.reset();
    X = std::make_shared<Tensor2d>("X (Ib|Mj)", naoccA, navirB, naoccA, naoccB);
    X->gemm(false, false, T, Y, 1.0, 0.0);
    T.reset();
    Y.reset();
    WL->sort(3214, X, -1.0, 1.0);
    X.reset();
    WL->write(psio_, PSIF_DFOCC_AMPS);
//WL->print();
    WL.reset();    
//outfile->Printf("\tI am here.\n");

} // end ccsdl_WMbIj_ABAB

void DFOCC::uccsdl_WmBiJ_BABA()
{
    SharedTensor2d G, W, WL, X, Y, Z, Z_, T2, Tau, J, K, T;

    // W_mBiJ = <Bm|Ji> + \sum_{e} t_Ji^Be Ft_me - \sum_{N} t_N^B W_mNiJ + \sum_{e,F} Tau_Ji^Fe <Bm|Fe> + \sum_{e} t_i^e Z_mBeJ 
    // - \sum_{E} t_J^E Z_mBEi + \sum_{N,E} t_JN^BE <Nm|Ei> + \sum_{n,e} t_Jn^Be <mn||ie> - \sum_{N,e} t_Ni^Be <Nm|Je>                       (112)                                 
    // W_mBiJ (1) += <Bm|Ji>
    J = std::make_shared<Tensor2d>("J (mi|JB)", naoccB, naoccB, naoccA, navirA);
    J->gemm(true, false, bQijB, bQiaA, 1.0, 0.0);
    WL = std::make_shared<Tensor2d>("WL (mB|iJ)", naoccB, navirA, naoccB, naoccA);
    WL->sort(1423, J, 1.0, 0.0);
    J.reset();
    // W_mBiJ (2) += \sum_{e} t_Ji^Be Ft_me
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X (Ji|Bm)", naoccA, naoccB, navirA, naoccB);
    X->contract(false, true, naoccA * naoccB * navirA, naoccB, navirB, T2, FiaB, 1.0, 0.0);
    T2.reset();
    WL->sort(4321, X, 1.0, 1.0);
    X.reset();
    // W_mBiJ (3) -= \sum_{N} t_N^B W_mNiJ = \sum_{N} t_N^B W_NmJi 
    uccsdl_WMnIj_ABAB();
    W = std::make_shared<Tensor2d>("W <Mn|Ij>", naoccA, naoccB, naoccA, naoccB);
    W->read(psio_, PSIF_DFOCC_AMPS);
    X = std::make_shared<Tensor2d>("X <mJ|iB>", naoccB, naoccA, naoccB, navirA);
    X->contract(true, false, naoccB * naoccA * naoccB, navirA, naoccA, W, t1A, 1.0, 0.0);
    W.reset();
    WL->sort(1432, X, -1.0, 1.0);
    X.reset();
    // W_mBiJ (4) += \sum_{e,F} Tau_Ji^Fe <Bm|Fe>
    Tau = std::make_shared<Tensor2d>("Tau <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    uccsd_tau_amps_OS(naoccA, naoccB, navirA, navirB, Tau, T2, t1A, t1B);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (BF|me)", navirA, navirA, naoccB, navirB);
    J->gemm(true, false, bQabA, bQiaB, 1.0, 0.0);
    G = std::make_shared<Tensor2d>("G (Bm|Fe)", navirA, naoccB, navirA, navirB);
    G->sort(1324, J, 1.0, 0.0);
    J.reset();
    Y = std::make_shared<Tensor2d>("Y (Bm|Ji)", navirA, naoccB, naoccA, naoccB);
    Y->gemm(false, true, G, Tau, 1.0, 0.0);
    G.reset();
    Tau.reset();
    WL->sort(2143, Y, 1.0, 1.0);
    Y.reset();
    // W_mBiJ (5) += \sum_{e} t_i^e Z_mBeJ
    uccsdl_ZmBeJ_BABA();
    Y = std::make_shared<Tensor2d>("Z (me|JB)", naoccB, navirB, naoccA, navirA);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Z = std::make_shared<Tensor2d>("Z (JB|me)", naoccA, navirA, naoccB, navirB);
    Z->trans(Y);
    Y.reset();
    X = std::make_shared<Tensor2d>("X (JB|mi)", naoccA, navirA, naoccB, naoccB);
    X->contract(false, true, naoccA * navirA * naoccB, naoccB, navirB, Z, t1B, 1.0, 0.0);
    Z.reset();
    WL->sort(3241, X, 1.0, 1.0);
    X.reset();
    // W_mBiJ (6) -= \sum_{E} t_J^E Z_mBEi 
    uccsdl_ZmBEj_BAAB();
    Z = std::make_shared<Tensor2d>("Z (mE|jB)", naoccB, navirA, naoccB, navirA);
    Z->read(psio_, PSIF_DFOCC_AMPS);
    Y = std::make_shared<Tensor2d>("Z (iB|mE)", naoccB, navirA, naoccB, navirA);
    Y->trans(Z);
    Z.reset();
    X = std::make_shared<Tensor2d>("X (iB|mJ)", naoccB, navirA, naoccB, naoccA);
    X->contract(false, true, naoccB * navirA * naoccB, naoccA, navirA, Y, t1A, 1.0, 0.0);
    Y.reset();
    WL->sort(3214, X, -1.0, 1.0);
    X.reset();
    // W_mBiJ (7) += \sum_{N,E} t_JN^BE <Nm|Ei>
    T2 = std::make_shared<Tensor2d>("T2 <IJ|AB>", naoccA, naoccA, navirA, navirA);
    T2->read_anti_symm(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <JB|NE>", naoccA, navirA, naoccA, navirA);   
    T->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    J = std::make_shared<Tensor2d>("J (NE|mi)", naoccA, navirA, naoccB, naoccB);
    J->gemm(true, false, bQiaA, bQijB, 1.0, 0.0);
    X = std::make_shared<Tensor2d>("X (JB|mi)", naoccA, navirA, naoccB, naoccB);
    X->gemm(false, false, T, J, 1.0, 0.0);  
    T.reset();
    J.reset();
    WL->sort(3241, X, 1.0, 1.0);
    X.reset();
    // W_mBiJ (8) += \sum_{n,e} t_Jn^Be <mn||ie>
    // <mn||ie> = <mn|ie> - <mn|ei>
    // K_ (mi,ne) = (mi|ne) - (me|ni) --------> (me|ni) = (ni|me) sort:3214  
    J = std::make_shared<Tensor2d>("J (mi|ne)", naoccB, naoccB, naoccB, navirB);
    J->gemm(true, false, bQijB, bQiaB, 1.0, 0.0);
    K = std::make_shared<Tensor2d>("K (mi|ne)", naoccB, naoccB, naoccB, navirB);
    K->sort(3214, J, -1.0, 0.0);
    K->axpy(J, 1.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T (JB|ne)", naoccA, navirA, naoccB, navirB);
    T->sort(1324, T2, 1.0, 0.0);
    T2.reset();
    X= std::make_shared<Tensor2d>("X (mi|JB)", naoccB, naoccB, naoccA, navirA);
    X->gemm(false, true, K, T, 1.0, 0.0);
    K.reset();
    T.reset();
    WL->sort(1423, X, 1.0, 1.0);
    X.reset();
    // W_mBiJ (9) -= \sum_{N,e} t_Ni^Be <Nm|Je>
    J = std::make_shared<Tensor2d>("J (NJ|me)", naoccA, naoccA, naoccB, navirB);
    J->gemm(true, false, bQijA, bQiaB, 1.0, 0.0);
    Y = std::make_shared<Tensor2d>("Y (Ne|mJ)", naoccA, navirB, naoccB, naoccA);
    Y->sort(1432, J, 1.0, 0.0);
    J.reset();
    T2 = std::make_shared<Tensor2d>("T2 <Ij|Ab>", naoccA, naoccB, navirA, navirB);
    T2->read(psio_, PSIF_DFOCC_AMPS);
    T = std::make_shared<Tensor2d>("T <iB|Ne>", naoccB, navirA, naoccA, navirB);
    T->sort(2314, T2, 1.0, 0.0);
    T2.reset();
    X = std::make_shared<Tensor2d>("X (iB|mJ)", naoccB, navirA, naoccB, naoccA);
    X->gemm(false, false, T, Y, 1.0, 0.0);
    T.reset();
    Y.reset();
    WL->sort(3214, X, -1.0, 1.0);
    WL->write(psio_, PSIF_DFOCC_AMPS);
    X.reset();
    WL.reset();
//outfile->Printf("\tI am here.\n");

/*
    // Control

    SharedTensor2d A, B;

    X = std::make_shared<Tensor2d>("WL (MB|IJ)", naoccA, navirA, naoccA, naoccA);
    X->read(psio_, PSIF_DFOCC_AMPS);
    X->print();
    X.reset();

    Y = std::make_shared<Tensor2d>("WL (mb|ij)", naoccB, navirB, naoccB, naoccB);
    Y->read(psio_, PSIF_DFOCC_AMPS);
    Y->print(); 
    Y.reset();

    Z = std::make_shared<Tensor2d>("WL (Mb|Ij)", naoccA, navirB, naoccA, naoccB);
    Z->read(psio_, PSIF_DFOCC_AMPS);
    B = std::make_shared<Tensor2d>("WL (mB|iJ)", naoccB, navirA, naoccB, naoccA);
    B->read(psio_, PSIF_DFOCC_AMPS);
    A = std::make_shared<Tensor2d>("A (Mb|Ij)", naoccA, navirB, naoccA, naoccB);
    A->copy(Z);
    Z.reset();
    A->sort(1243, B, -1.0, 1.0);
    B.reset();
    A->print();
    A.reset();
*/
} // end ccsdl_WmBiJ_BABA

}  // namespace dfoccwave
}  // namespace psi

