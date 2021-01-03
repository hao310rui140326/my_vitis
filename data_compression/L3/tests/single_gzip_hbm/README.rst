====================
GZIP HBM Application
====================

This section presents brief introduction about Gzip HBM application and step by step
guidelines to build and deployment.

Overview
--------

GZIP is an Open Source data compression library* which provides
high compression ratio compared to Limpel Ziev based data compression algorithms
(Byte Compression). It applies two levels of compression,

*  Byte Level (Limpel Ziev  LZ Based Compression Scheme)
*  Bit Level (Huffman Entropy)

Due to its high compression ratio it takes higher precedence over LZ based
compression schemes. Traditionally the CPU based solutions are limited to MB/s
speed but there is a high demand for accelerated ZLIB which provides throughput
in terms of GB/s. 

This demo is aimed at showcasing Xilinx Alveo U50 (HBM Platform) acceleration of ZLIB for both
compression and decompression. 

.. code-block:: bash

   Tested Tool: 2020.2 
   Tested XRT: 2020.2
   Tested XSA: xilinx_u50_gen3x16_xdma_201920_3 


Executable Usage
----------------

This application is present under ``L3/tests/single_gzip_hbm/`` directory. Follow build instructions to generate executable and binary.

The host executable generated is named as "**xgzip.exe**" and it is generated in ``./build`` directory.

Following is the usage of the executable:

1. To execute single file for compression 	      : ``./build/xgzip.exe -sx ./build/xclbin_<xsa_name>_<TARGET mode>/xgzip.xclbin -c <input file_name>``
2. To execute single file for decompression           : ``./build/xgzip.exe -sx ./build/xclbin_<xsa_name>_<TARGET mode>/xgzip.xclbin -d <compressed file_name>``
3. To validate single file (compress & decompress)    : ``./build/xgzip.exe -sx ./build/xclbin_<xsa_name>_<TARGET mode>/xgzip.xclbin -t <input file_name>``
4. To validate multiple files (compress & decompress) : ``./build/xgzip.exe -sx ./build/xclbin_<xsa_name>_<TARGET mode>/xgzip.xclbin -l <files.list>``

	- ``<files.list>``: Contains various file names with current path

The usage of the generated executable is as follows:

.. code-block:: bash
  Usage: ./build_dir.sw_emu.xilinx_u50_gen3x16_xdma_201920_3/xgzip.exe [Options] [Files] 

          --help,              -h        Print Help Options
          --compress,          -c        Compress
          --decompress,        -d        DeCompress
          --single_xclbin,     -sx       Single XCLBIN [Optional]
          --test,              -t        Compress Decompress
          --c_file_list,       -cfl      Compress list files
          --d_file_list,       -dfl      Decompress list files
          --file_list,         -l        List of Input Files
          --zlib,              -zlib     [0:GZIP, 1:ZLIB]             Default: [0]
          --ck,                -ck       Compress CU                  Default: [0]
          --dk,                -dk       Decompress CU                Default: [0]
          --id,                -id       Device ID                    Default: [0]
          --max_cr,            -mcr      Maximum CR                   Default: [10]

===========================================================
