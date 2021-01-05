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
#include "fft_top.h"

#include "hls_linear_algebra.h"

#include <complex>

int main(int argc, char** argv) {
	/////////////////////////////////////////////////////////////////////////////////////////
	//define
	//T_elemType l_inMat[k_fftKernelSize][k_fftKernelSize];
	static cmpxDataIn l_inMat[k_fftKernelSize][k_fftKernelSize];
	static cmpxDataIn l_inMat2[k_fftKernelSize][k_fftKernelSize];
        //std::complex<float>  multi_inMat1[k_fftKernelSize][k_fftKernelSize];
    	static cmpxDataIn xinmat[k_fftKernelSize*k_fftKernelSize];
    	static cmpxDataOut xoutmat[k_fftKernelSize*k_fftKernelSize];
	FILE *dafile =fopen("./ori_data.txt","w");

	////////////////////////////////////////////////////////////////////////////////////////////////
	//init
	  // init input matrix with real part only impulse
   	 for (int r = 0; r < k_fftKernelSize; ++r) {
   	     for (int c = 0; c < k_fftKernelSize; ++c) {
   	         //if (r == 0 && c == 0)
   	         //    l_inMat[r][c] = T_compleFloat(1, 0);
   	         //else
   	         //    l_inMat[r][c] = T_compleFloat(0, 0);
   	         l_inMat[r][c]  = cmpxDataIn(r, c);
   	         l_inMat2[r][c] = cmpxDataIn(c, r);
   	     }
   	 }

	////////////////////////////////////////////////////////////////////////////////////////////////
	//map to mem 
	for (int r = 0; r < k_fftKernelSize; ++r) {
   	     for (int c = 0; c < k_fftKernelSize; ++c) {
   	         xinmat[r*k_fftKernelSize+c]=cmpxDataIn(l_inMat[r][c].real(),l_inMat[r][c].imag());
		 //std::cout << " xinmat.real: "<<xinmat[r*k_fftKernelSize+c].real()<<" xinmat.imag: "<<xinmat[r*k_fftKernelSize+c].imag()<< std::endl;
     		 fprintf(dafile,"%f %f \n",xinmat[r*k_fftKernelSize+c].real() ,xinmat[r*k_fftKernelSize+c].imag() );		 
   	     }
   	 }
    	std::cout << "---------------------FFT 2D start " << std::endl;
	fclose(dafile);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//fft_2d
	//fft_2d(0,xinmat,xoutmat);
	//fft_2d(0,xinmat);
	//mfft_2d(0,1,l_inMat);
	int  row,col ;
	mdsp(l_inMat,l_inMat2,&row,&col);

	FILE *outfile =fopen("./fft2d_out.txt","w");
	FILE *outfile2 =fopen("./fft2d2_out.txt","w");

for (int idx1=0;idx1<k_fftKernelSize;++idx1)
{
	for (int idx2=0;idx2<k_fftKernelSize;++idx2)
	{
       		//fprintf(outfile,"%f +j*%f  ",xinmat[idx1*FFT_LENGTH+idx2].real(),xinmat[idx1*FFT_LENGTH+idx2].imag());										
       		fprintf(outfile,"%f +j*%f  ",l_inMat[idx1][idx2].real(),l_inMat[idx1][idx2].imag());										
       		fprintf(outfile2,"%f +j*%f  ",l_inMat2[idx1][idx2].real(),l_inMat2[idx1][idx2].imag());										
	}
       		fprintf(outfile," \n");										
       		fprintf(outfile2," \n");										
}


fclose(outfile);
fclose(outfile2);

    std::cout << "======OUTPUT : row is  " << row <<" col is "<< col <<std::endl;
	
    std::cout << "================================================================" << std::endl;
    std::cout << "---------------------Impulse test Passed Successfully." << std::endl;
    std::cout << "================================================================" << std::endl;
    return 0;

}
#endif
