#
# Copyright 2019 Xilinx, Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

source settings.tcl
set DIR_NAME "zlib_treegen"
set DESIGN_PATH "${XF_PROJ_ROOT}/L1/tests/${DIR_NAME}"
set PROJ "zlib_treegen_test.prj"
set SOLN "sol1"

if {![info exists CLKP]} {
  set CLKP 3.3
}

open_project -reset $PROJ

add_files $XF_PROJ_ROOT/common/libs/logger/logger.cpp -cflags "-I${XF_PROJ_ROOT}/common/libs/logger"
add_files $XF_PROJ_ROOT/common/libs/cmdparser/cmdlineparser.cpp -cflags "-I${XF_PROJ_ROOT}/common/libs/cmdparser -I${XF_PROJ_ROOT}/common/libs/logger"
add_files zlib_treegen_test.cpp -cflags "-I${XF_PROJ_ROOT}/L1/include/hw -I${XF_PROJ_ROOT}/L2/include -I${XF_PROJ_ROOT}/common/libs/cmdparser -I${XF_PROJ_ROOT}/common/libs/logger"

add_files -tb $XF_PROJ_ROOT/common/libs/logger/logger.cpp -cflags "-I${XF_PROJ_ROOT}/common/libs/logger"
add_files -tb $XF_PROJ_ROOT/common/libs/cmdparser/cmdlineparser.cpp -cflags "-I${XF_PROJ_ROOT}/common/libs/cmdparser -I${XF_PROJ_ROOT}/common/libs/logger"
add_files -tb zlib_treegen_test.cpp -cflags "-I${XF_PROJ_ROOT}/L1/include/hw -I${XF_PROJ_ROOT}/L2/include -I${XF_PROJ_ROOT}/common/libs/cmdparser"


# Set the top-level function
set_top treegenMM

# Create a solution
open_solution -reset $SOLN

# Define technology and clock rate
set_part $XPART
create_clock -period $CLKP

config_compile -pragma_strict_mode

if {$CSIM == 1} {
  csim_design -argv "-l ${DESIGN_PATH}/test.list -p ${DESIGN_PATH}"
}

if {$CSYNTH == 1} {
  csynth_design
}

if {$COSIM == 1} {
  cosim_design -argv "-o ${DESIGN_PATH}/sample.txt -p ${DESIGN_PATH}"
}

if {$VIVADO_SYN == 1} {
  export_design -flow syn -rtl verilog
}

if {$VIVADO_IMPL == 1} {
  export_design -flow impl -rtl verilog
}

exit
