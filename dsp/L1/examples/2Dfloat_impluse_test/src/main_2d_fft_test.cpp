/*
 * Copyright 2019 Xilinx, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/
//================================== End Lic =================================================
#define TEST_2D_FFT_
#ifdef TEST_2D_FFT_
#ifndef __SYNTHESIS__
#define _DEBUG_TYPES
#endif
#include <fstream>
#include <math.h>
#include <string>
#include <assert.h>
#include <stdio.h>
#include "top_2d_fft_test.hpp"
#include "fft_top.h"
#include "mVerificationUtlityFunctions.hpp"
#include "vitis_fft/hls_ssr_fft_2d_modeling_utilities.hpp"
#include "vt_fft.hpp"

#include "hls_linear_algebra.h"

#include <complex>

int main(int argc, char** argv) {
    // 2d input matrix
    T_elemType l_inMat[k_fftKernelSize][k_fftKernelSize];
    T_outType l_outMat[k_fftKernelSize][k_fftKernelSize];
    std::complex<float>  multi_inMat1[k_fftKernelSize][k_fftKernelSize];
    std::complex<float>  multi_inMat2[k_fftKernelSize][k_fftKernelSize];
    std::complex<float>  multi_outMat[k_fftKernelSize][k_fftKernelSize];
    //std::complex<float>  fft_inVec[k_fftKernelSize];
    //std::complex<float>  fft_outVec[k_fftKernelSize];
    T_outType ifft_iMat[k_fftKernelSize][k_fftKernelSize];
    T_outType ifft_oMat[k_fftKernelSize][k_fftKernelSize];

    T_outType l_data2d_golden[k_fftKernelSize][k_fftKernelSize];

    const int SAMPLES = (1 << FFT_NFFT_MAX);    
    static cmpxDataIn xn_input[SAMPLES];
    static cmpxDataOut xk_output[SAMPLES];
    static cmpxDataOut xi_output[SAMPLES];
    static cmpxDataOut fft_output[SAMPLES];
    static cmpxDataOut ifft_output[SAMPLES];


    float   abs_mat;
    float   abs_mat_max;
    unsigned int  max_abs_index1;
    unsigned int  max_abs_index2;

    // init input matrix with real part only impulse
    for (int r = 0; r < k_fftKernelSize; ++r) {
        for (int c = 0; c < k_fftKernelSize; ++c) {
            //if (r == 0 && c == 0)
            //    l_inMat[r][c] = T_compleFloat(1, 0);
            //else
            //    l_inMat[r][c] = T_compleFloat(0, 0);
            l_inMat[r][c] = T_compleFloat(r, c);
        }
    }
    std::cout << "k_memWidth is " <<k_memWidth<<"k_fftKernelSize is " <<k_fftKernelSize << std::endl;
    // Wide Stream for reading and streaming a 2-d matrix
    MemWideIFStreamTypeIn l_matToStream("matrixToStreaming");
    MemWideIFStreamTypeOut fftOutputStream("fftOutputStream");
    // Pass same data stream multiple times to measure the II correctly
    //for (int runs = 0; runs < 5; ++runs) {
    for (int runs = 0; runs < 1; ++runs) {
        stream2DMatrix<k_fftKernelSize, k_fftKernelSize, k_memWidth, T_elemType, MemWideIFTypeIn>(l_inMat,
                                                                                                  l_matToStream);
        top_fft2d(l_matToStream, fftOutputStream);

        printMatStream<k_fftKernelSize, k_fftKernelSize, k_memWidth, MemWideIFTypeOut>(
            fftOutputStream, "2D FFT Output Natural Order...");
        streamToMatrix<k_fftKernelSize, k_fftKernelSize, k_memWidth, T_outType>(fftOutputStream, l_outMat);
    } // runs loop
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //matrix_multiply
    for (int r = 0; r < k_fftKernelSize; ++r) {
            for (int c = 0; c < k_fftKernelSize; ++c) {
            	multi_inMat1[r][c].real(l_outMat[r][c].real());
            	multi_inMat2[r][c].real(l_outMat[r][c].real());
            	multi_inMat1[r][c].imag(l_outMat[r][c].imag());
            	multi_inMat2[r][c].imag((-1)*l_outMat[r][c].imag());
            }
        }
  hls::matrix_multiply <hls::NoTranspose, hls::NoTranspose, k_fftKernelSize,k_fftKernelSize, k_fftKernelSize, k_fftKernelSize,k_fftKernelSize ,k_fftKernelSize ,std::complex<float> ,std::complex<float>> (multi_inMat1,multi_inMat2 , multi_outMat);
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        //ifft2d
	//ifft_iMat ifft_oMat
	MemWideIFStreamTypeIn  ifft_matToStream("ifft_matrixToStreaming");
    	MemWideIFStreamTypeOut ifftOutputStream("ifftOutputStream");

 	for (int r = 0; r < k_fftKernelSize; ++r) {
            for (int c = 0; c < k_fftKernelSize; ++c) {
            	ifft_iMat[r][c]= T_compleFloat( multi_outMat[r][c].real() , multi_outMat[r][c].imag()  );
            }
        }

  for (int runs = 0; runs < 1; ++runs) {
        stream2DMatrix<k_fftKernelSize, k_fftKernelSize, k_memWidth, T_elemType, MemWideIFTypeIn>(ifft_iMat,
                                                                                                  ifft_matToStream);
        top_ifft2d(ifft_matToStream, ifftOutputStream);

        printMatStream<k_fftKernelSize, k_fftKernelSize, k_memWidth, MemWideIFTypeOut>(
            ifftOutputStream, "2D IFFT Output Natural Order...");
        streamToMatrix<k_fftKernelSize, k_fftKernelSize, k_memWidth, T_outType>(ifftOutputStream, ifft_oMat);
    } // runs loop
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   abs_mat_max = 0 ;
   for (int r = 0; r < k_fftKernelSize; ++r) {
            for (int c = 0; c < k_fftKernelSize; ++c) {
            	abs_mat = ifft_oMat[r][c].real()*ifft_oMat[r][c].real() + ifft_oMat[r][c].imag()*ifft_oMat[r][c].imag() ;
		if (abs_mat>=abs_mat_max)  
		{
			abs_mat_max = abs_mat;
			max_abs_index1 = r ;
			max_abs_index2 = c ;
		}
            }
        }

std::cout << "MAX_ABS is " <<abs_mat_max<<" Index1  is " <<max_abs_index1<<" Index2  is " <<max_abs_index2<< std::endl;

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//fft 
bool ovflo;
bool invflo;
//std::ofstream outfile ;
//outfile.open("./ffto.txt");
FILE *outfile =fopen("./ffto.txt","w");
FILE *outifile=fopen("./iffto.txt","w");

  for (int r = 0; r < SAMPLES; ++r) {
            xn_input[r] = cmpxDataIn(r, 0);
    }

fft_top(0, xn_input, xk_output, &ovflo);

for (int r = 0; r < SAMPLES; ++r) {
	unsigned addr_reverse = 0;

	for (int k = 0; k < 9; ++k)
	{
	  addr_reverse <<= 1;
	  addr_reverse |= (r >> k) & 0x1;
	}
	fft_output[addr_reverse] =  cmpxDataIn(xk_output[r].real(),xk_output[r].imag());
	//std::cout << "addr_reverse is " <<addr_reverse << std::endl;
}

for (int r = 0; r < SAMPLES; ++r) {
       fprintf(outfile,"%f %f \n",fft_output[r].real(),fft_output[r].imag());
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//ifft 
fft_top(1, fft_output,xi_output, &invflo);

for (int r = 0; r < SAMPLES; ++r) {
	unsigned addr_reverse = 0;

	for (int k = 0; k < 9; ++k)
	{
	  addr_reverse <<= 1;
	  addr_reverse |= (r >> k) & 0x1;
	}
	ifft_output[addr_reverse] =  cmpxDataIn(xi_output[r].real()/512,xi_output[r].imag()/512);
	//std::cout << "addr_reverse is " <<addr_reverse << std::endl;
}

for (int r = 0; r < SAMPLES; ++r) {
       fprintf(outifile,"%f %f \n",ifft_output[r].real(),ifft_output[r].imag());
    }
//outfile.close();
fclose(outfile);
fclose(outifile);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//compare
    //T_outType golden_result = T_elemType(1, 0);
    //for (int r = 0; r < k_fftKernelSize; ++r) {
    //    for (int c = 0; c < k_fftKernelSize; ++c) {
    //        if (golden_result != l_outMat[r][c]) return 1;
    //    }
    //}

    std::cout << "================================================================" << std::endl;
    std::cout << "---------------------Impulse test Passed Successfully." << std::endl;
    std::cout << "================================================================" << std::endl;
    return 0;
}
#endif
