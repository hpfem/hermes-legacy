if(NOT H2D_COMPLEX)
    return()
endif(NOT H2D_COMPLEX)
project(test-examples-schroedinger-gross-pitaevski-adapt)

add_executable(${PROJECT_NAME} main.cpp ../definitions.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME} ${HERMES2D_CPLX})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-tutorial-P06-05-gross-pitaevski-adapt ${BIN})
set_tests_properties(test-tutorial-P06-05-gross-pitaevski-adapt PROPERTIES LABELS slow)

