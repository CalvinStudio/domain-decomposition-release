#!/bin/bash
#
PROJ_PATH="fd3d_proj_small/cube100/"
#mpi node
NODE_INFO="node5 slots=4"
# NODE_INFO="node5 slots=2 \nnode5 slots=2"
EXEC_FILE="rtm3d_blend_elastic_cpml_model_mpi"
#
RTM_TEST_PATH="/home/calvin/WKSP/domain-decomposition-release/src/jarvis_cpp/jarvis/.jarvis_2_app/.jarvis_3_seismology/elastic_1_rtm3d_model_mpi/rtm3d_2_model_mpi_blend_test/"
#cu file
CUDA_SRC[0]="${RTM_TEST_PATH}rtm3d_blend_elastic_model_mpi_device.cu"
#cpp file
CPP_SRC[0]="${RTM_TEST_PATH}rtm3d_blend_elastic_model_mpi_main.cpp"
#
JARVIS_SH_PATH="/home/calvin/WKSP/domain-decomposition-release/build/jarvis/"
. ${JARVIS_SH_PATH}jarvis_0_env
. ${JARVIS_SH_PATH}jarvis_1_make
. ${JARVIS_SH_PATH}jarvis_2_run
. ${JARVIS_SH_PATH}jarvis_3_del
#
MF_FILE="${PROJ_PATH}.mpi_config"
echo -e $NODE_INFO > $MF_FILE
#
jarvis_logo
#
if [ "$1" == "del" ]; then
    jarvis_del "${CUDA_SRC[*]}" "${CPP_SRC[*]}"
fi
#
if [ "$1" == "make" ]; then
    jarvis_make_cuda_cpp_mpi "${CUDA_SRC[*]}" "${CPP_SRC[*]}" "${EXEC_FILE}"
fi
#
if [ "$1" == "run" ]; then
    jarvis_run_mpi "${EXEC_FILE}" "${PROJ_PATH}" "${MF_FILE}"
fi
#
if [ "$1" == "gdb" ]; then
    gdb    ${JARVIS_BIN_PATH}${EXEC_FILE}
fi
#
if [ "$1" == "" ]; then
    jarvis_del "${CUDA_SRC[*]}" "${CPP_SRC[*]}"
    jarvis_make_cuda_cpp_mpi "${CUDA_SRC[*]}" "${CPP_SRC[*]}" "${EXEC_FILE}"
    jarvis_run_mpi "${EXEC_FILE}" "${PROJ_PATH}" "${MF_FILE}"
fi