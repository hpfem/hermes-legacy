project(test-example-heat-transfer-wall-on-fire-adapt-time-only)

add_executable(${PROJECT_NAME} main.cpp ../definitions.cpp)
include (${hermes2d_SOURCE_DIR}/CMake.common)
set_common_target_properties(${PROJECT_NAME})
set(BIN ${PROJECT_BINARY_DIR}/${PROJECT_NAME})
add_test(test-example-heat-transfer-wall-on-fire-adapt-time-only ${BIN})
set_tests_properties(test-example-heat-transfer-wall-on-fire-adapt-time-only PROPERTIES LABELS slow)

