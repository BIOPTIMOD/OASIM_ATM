cmake_minimum_required(VERSION 3.0)
 
add_executable(test test.f90)
set_target_properties(test PROPERTIES SUFFIX ".x")
target_link_libraries(test PRIVATE oasim)
