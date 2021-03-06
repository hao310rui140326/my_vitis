/*******************************************************************************
Vendor: Xilinx 
Associated Filename: fft_top.h
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


#ifndef _H_FFT_TEST_H_
#define _H_FFT_TEST_H_

#include <hls_stream.h>

#include "ap_fixed.h"
#include "hls_fft.h"

// configurable params
//const char FFT_INPUT_WIDTH                     = 16;
//const char FFT_OUTPUT_WIDTH                    = FFT_INPUT_WIDTH;
const char FFT_CONFIG_WIDTH                    = 16;
//const char FFT_CONFIG_WIDTH                    = 8;
const char FFT_NFFT_MAX                        = 9; 
//const char FFT_NFFT_MAX                        = 4; 
const int  FFT_LENGTH                          = 1 << FFT_NFFT_MAX; 

const int k_fftKernelSize = FFT_LENGTH;


#include <complex>
using namespace std;

struct config1 : hls::ip_fft::params_t {
    static const unsigned ordering_opt = hls::ip_fft::natural_order;
    static const unsigned config_width = FFT_CONFIG_WIDTH;
    static const unsigned phase_factor_width = 25;
    static const unsigned max_nfft = FFT_NFFT_MAX;
    static const unsigned stages_block_ram = 0;
    static const unsigned arch_opt = hls::ip_fft::pipelined_streaming_io;
};

typedef hls::ip_fft::config_t<config1> config_t;
typedef hls::ip_fft::status_t<config1> status_t;

//typedef ap_fixed<FFT_INPUT_WIDTH,1> data_in_t;
//typedef ap_fixed<FFT_OUTPUT_WIDTH,FFT_OUTPUT_WIDTH-FFT_INPUT_WIDTH+1> data_out_t;

typedef float  data_in_t;
typedef float  data_out_t;

typedef std::complex<data_in_t> cmpxDataIn;
typedef std::complex<data_out_t> cmpxDataOut;

void dummy_proc_fe(
    bool direction,
    config_t* config, 
    cmpxDataIn in[FFT_LENGTH], 
    cmpxDataIn out[FFT_LENGTH]);

void dummy_proc_be(
    status_t* status_in, 
    bool* ovflo,
    cmpxDataOut in[FFT_LENGTH], 
    cmpxDataOut out[FFT_LENGTH]);

void fft_top(
    bool direction,
    cmpxDataIn in[FFT_LENGTH],
    cmpxDataOut out[FFT_LENGTH],
    bool* ovflo);

void fft_2d(
    bool direction,		
    cmpxDataIn in[FFT_LENGTH*FFT_LENGTH]);//,

void mfft_2d(
    bool direction,		
    bool conj,		
    cmpxDataIn in[FFT_LENGTH][FFT_LENGTH]);//,

void mifft_2d(
		complex<data_in_t> in[FFT_LENGTH][FFT_LENGTH],
     		//hls::stream<float> &strm_real_out,
     		//hls::stream<float> &strm_imag_out
     		hls::stream<cmpxDataIn> &strm_out
		);

void max_abs(
		complex<data_in_t> in[FFT_LENGTH][FFT_LENGTH],
		int *shift_row,
		int *shift_col
	);



void strm_max_abs(
		hls::stream<cmpxDataIn> &strm_out,
		//hls::stream<float> &strm_real_out,
     		//hls::stream<float> &strm_imag_out,
		int *shift_row,
		int *shift_col
	);

void MatrixMult(
		complex<data_in_t> in1[FFT_LENGTH][FFT_LENGTH],
		complex<data_in_t> in2[FFT_LENGTH][FFT_LENGTH],
		complex<data_in_t> out[FFT_LENGTH][FFT_LENGTH]
		);

void mdsp(
    cmpxDataIn in1[FFT_LENGTH][FFT_LENGTH],
    cmpxDataIn in2[FFT_LENGTH][FFT_LENGTH],
	cmpxDataIn matrix_mult[FFT_LENGTH][FFT_LENGTH],
	//hls::stream<float> &strm_real_out,
     	//hls::stream<float> &strm_imag_out
     	hls::stream<cmpxDataIn> &strm_out
	//int *shift_row,
	//int *shift_col
    );//,


#endif


