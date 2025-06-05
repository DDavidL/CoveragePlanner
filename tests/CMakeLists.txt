include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

add_executable(graph_base_test test_graph_base.cpp)

target_include_directories(graph_base_test PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(graph_base_test gtest_main)

add_test(NAME graph_base_test COMMAND graph_base_test)
