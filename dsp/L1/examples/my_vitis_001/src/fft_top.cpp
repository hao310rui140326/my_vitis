/*******************************************************************************
Vendor: Xilinx 
Associated Filename: fft_top.cpp
Purpose: Xilinx FFT IP-XACT IP in Vivado HLS
Revision History: September 26, 2013 - initial release
                                                
*******************************************************************************
#-  (c) Copyright 2011-2019 Xilinx, Inc. All rights reserved.
#-
#-  This file contains confidential and proprietary information
#-  of Xilinx, Inc. and is protected under U.S. and
#-  international copyright and other intellectual property
#-  laws.
#-
#-  DISCLAIMER
#-  This disclaimer is not a license and does not grant any
#-  rights to the materials distributed herewith. Except as
#-  otherwise provided in a valid license issued to you by
#-  Xilinx, and to the maximum extent permitted by applicable
#-  law: (1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND
#-  WITH ALL FAULTS, AND XILINX HEREBY DISCLAIMS ALL WARRANTIES
#-  AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, INCLUDING
#-  BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-
#-  INFRINGEMENT, OR FITNESS FOR ANY PARTICULAR PURPOSE; and
#-  (2) Xilinx shall not be liable (whether in contract or tort,
#-  including negligence, or under any other theory of
#-  liability) for any loss or damage of any kind or nature
#-  related to, arising under or in connection with these
#-  materials, including for any direct, or any indirect,
#-  special, incidental, or consequential loss or damage
#-  (including loss of data, profits, goodwill, or any type of
#-  loss or damage suffered as a result of any action brought
#-  by a third party) even if such damage or loss was
#-  reasonably foreseeable or Xilinx had been advised of the
#-  possibility of the same.
#-
#-  CRITICAL APPLICATIONS
#-  Xilinx products are not designed or intended to be fail-
#-  safe, or for use in any application requiring fail-safe
#-  performance, such as life-support or safety devices or
#-  systems, Class III medical devices, nuclear facilities,
#-  applications related to the deployment of airbags, or any
#-  other applications that could lead to death, personal
#-  injury, or severe property or environmental damage
#-  (individually and collectively, "Critical
#-  Applications"). Customer assumes the sole risk and
#-  liability of any use of Xilinx products in Critical
#-  Applications, subject only to applicable laws and
#-  regulations governing limitations on product liability.
#-
#-  THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS
#-  PART OF THIS FILE AT ALL TIMES. 
#- ************************************************************************


This file contains confidential and proprietary information of Xilinx, Inc. and 
is protected under U.S. and international copyright and other intellectual 
property laws.

DISCLAIMER
This disclaimer is not a license and does not grant any rights to the materials 
distributed herewith. Except as otherwise provided in a valid license issued to 
you by Xilinx, and to the maximum extent permitted by applicable law: 
(1) THESE MATERIALS ARE MADE AVAILABLE "AS IS" AND WITH ALL FAULTS, AND XILINX 
HEREBY DISCLAIMS ALL WARRANTIES AND CONDITIONS, EXPRESS, IMPLIED, OR STATUTORY, 
INCLUDING BUT NOT LIMITED TO WARRANTIES OF MERCHANTABILITY, NON-INFRINGEMENT, OR 
FITNESS FOR ANY PARTICULAR PURPOSE; and (2) Xilinx shall not be liable (whether 
in contract or tort, including negligence, or under any other theory of 
liability) for any loss or damage of any kind or nature related to, arising under 
or in connection with these materials, including for any direct, or any indirect, 
special, incidental, or consequential loss or damage (including loss of data, 
profits, goodwill, or any type of loss or damage suffered as a result of any 
action brought by a third party) even if such damage or loss was reasonably 
foreseeable or Xilinx had been advised of the possibility of the same.

CRITICAL APPLICATIONS
Xilinx products are not designed or intended to be fail-safe, or for use in any 
application requiring fail-safe performance, such as life-support or safety 
devices or systems, Class III medical devices, nuclear facilities, applications 
related to the deployment of airbags, or any other applications that could lead 
to death, personal injury, or severe property or environmental damage 
(individually and collectively, "Critical Applications"). Customer assumes the 
sole risk and liability of any use of Xilinx products in Critical Applications, 
subject only to applicable laws and regulations governing limitations on product 
liability. 

THIS COPYRIGHT NOTICE AND DISCLAIMER MUST BE RETAINED AS PART OF THIS FILE AT 
ALL TIMES.

*******************************************************************************/


#include "fft_top.h"

#include <stdio.h>

void dummy_proc_fe(
    bool direction,
    config_t* config, 
    cmpxDataIn in[FFT_LENGTH], 
    cmpxDataIn out[FFT_LENGTH])
{
    int i; 
    config->setDir(direction);
    config->setSch(0x2AB);
    for (i=0; i< FFT_LENGTH; i++)  {
        out[i] = cmpxDataIn(in[i].real(),in[i].imag());
	//std::cout <<"dummy_proc_fe.real "<<out[i].real()<<" dummy_proc_fe.imag: "<<out[i].imag()<< std::endl;								    
    }
}

void dummy_proc_be(
    status_t* status_in, 
    bool* ovflo,
    cmpxDataOut in[FFT_LENGTH], 
    cmpxDataOut out[FFT_LENGTH])
{
    int i; 
    for (i=0; i< FFT_LENGTH; i++)
        out[i] = in[i];
    *ovflo = status_in->getOvflo() & 0x1;
}

void proc_2d_fe(
    bool direction,
    config_t* config 
    )
{
    config->setDir(direction);
    config->setSch(0x2AB);
}



void fft_top(
    bool direction,
    complex<data_in_t> in[FFT_LENGTH],
    complex<data_out_t> out[FFT_LENGTH],
    bool* ovflo)
{
#pragma HLS interface ap_hs port=direction
#pragma HLS interface ap_fifo depth=1 port=ovflo
#pragma HLS interface ap_fifo depth=FFT_LENGTH port=in,out
#pragma HLS data_pack variable=in
#pragma HLS data_pack variable=out
#pragma HLS dataflow
    complex<data_in_t> xn[FFT_LENGTH];
    complex<data_out_t> xk[FFT_LENGTH];
    config_t fft_config;
    status_t fft_status;
   
    dummy_proc_fe(direction, &fft_config, in, xn);
    // FFT IP
    hls::fft<config1>(xn, xk, &fft_status, &fft_config);
    dummy_proc_be(&fft_status, ovflo, xk, out);
}



void fft_2d(
		bool direction,
		complex<data_in_t> in[FFT_LENGTH*FFT_LENGTH])//,
		//complex<data_in_t> out[FFT_LENGTH*FFT_LENGTH]
    //)
{
#pragma HLS INTERFACE ap_memory storage_type=ram_2p port=in
	#pragma HLS interface ap_hs port=direction
	//#pragma HLS interface ap_fifo depth=1 port=ovflo
	//#pragma HLS interface ap_fifo depth=FFT_LENGTH port=in
	//#pragma HLS data_pack variable=in1
	#pragma HLS dataflow
    	const int SAMPLES = (1 << FFT_NFFT_MAX);    	
	//static cmpxDataIn   fft_input[SAMPLES];
	//static cmpxDataOut  fft_output[SAMPLES];
	//static cmpxDataOut  fft_ioutput[SAMPLES];

	static cmpxDataIn   fft_input[SAMPLES];
	static cmpxDataOut  fft_output[SAMPLES];
	cmpxDataOut  fft_ioutput[SAMPLES];


	//complex<data_in_t>   fft_input[FFT_LENGTH];
	//complex<data_in_t>   fft_output[FFT_LENGTH];
	//complex<data_in_t>   fft_ioutput[FFT_LENGTH];
	
	static  int d =0;
	static  int r =0;
	static  int c =0;
	static  int t =0;
	static  int k =0;
	static  int s =0;

	config_t fft_config;
	status_t fft_status;

	bool fvflo;
	//FILE *outfile =fopen("./fft2d_out.txt","w");

	//proc_2d_fe(direction, &fft_config);

	for (d = 0; d < 2; ++d) {
	//{
		//std::cout << " ///////////////////////////////////////////  d_is  "<<d<< std::endl;				
		for (r = 0; r < FFT_LENGTH; ++r) {
		//{
			//std::cout << " ///////////////////////////////////////////  fft_start  "<<r<< std::endl;							
			//fft
			for (c = 0; c < FFT_LENGTH; ++c ){
				if (d==0) {
					fft_input[c]=cmpxDataIn(in[r+c*FFT_LENGTH].real(),in[r+c*FFT_LENGTH].imag());	
				}
				else {
					fft_input[c]=cmpxDataIn(in[c+r*FFT_LENGTH].real(),in[c+r*FFT_LENGTH].imag());	
				}
		 	//	std::cout <<"num "<<c<< "  fft_input.real: "<<fft_input[c].real()<<" fft_input.imag: "<<fft_input[c].imag()<< std::endl;
       				//fprintf(outfile,"%f %f \n",fft_input[c].real(),fft_input[c].imag());
			}
			//fclose(outfile);

			fft_top(0, fft_input,fft_output, &fvflo);
			//hls::fft<config1>(fft_input, fft_output, &fft_status, &fft_config);

			//std::cout << " ///////////////////////////////////////////  fft_end  "<<r<< std::endl;										
			//addr inverse 
			//for (t = 0; t < FFT_LENGTH; ++t ){
			//	unsigned addr_reverse = 0;
			//	for (k = 0; k < FFT_NFFT_MAX; ++k)
			//	{
			//	  addr_reverse <<= 1;
			//	  addr_reverse |= (t >> k) & 0x1;
			//	}
			//	fft_ioutput[addr_reverse] =  cmpxDataOut(fft_output[t].real(),fft_output[t].imag());
		 	//	std::cout <<"num "<<t << " fft_output.real: "<<fft_output[t].real()<<" fft_output.imag: "<<fft_output[t].imag()<< std::endl;
       			//	//fprintf(outfile,"%f %f \n",fft_output[t].real(),fft_output[t].imag());								
			//}
			//write back 
			for ( s = 0; s < FFT_LENGTH; ++s ){
				if (d==0) {
					in[r+s*FFT_LENGTH]=cmpxDataIn(fft_output[s].real(),fft_output[s].imag());	
				}
				else {
					in[s+r*FFT_LENGTH]=cmpxDataIn(fft_output[s].real(),fft_output[s].imag());	
				}
			}
		}
	}


//for (int idx1=0;idx1<FFT_LENGTH;++idx1)
//{
//	for (int idx2=0;idx2<FFT_LENGTH;++idx2)
//	{
//       		fprintf(outfile,"%f +j*%f  ",in[idx1*FFT_LENGTH+idx2].real(),in[idx1*FFT_LENGTH+idx2].imag());										
//	}
//       		fprintf(outfile," \n");										
//}
//
//
//fclose(outfile);

}




