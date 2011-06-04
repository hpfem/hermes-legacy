project(test-tutorial-P04-06-nonlinear)

add_executable(${PROJECT_NAME} main.cpp ../definitions.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-tutorial-P04-06-nonlinear ${BIN})
set_tests_properties(test-tutorial-P04-06-nonlinear PROPERTIES LABELS slow)
