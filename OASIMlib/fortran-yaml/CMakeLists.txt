cmake_minimum_required(VERSION 3.0)

add_library(yaml_oasim OBJECT yaml_oasim_types.F90 yaml_oasim_settings.F90 yaml_oasim.F90)
set_property(TARGET yaml_oasim PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET yaml_oasim PROPERTY COMPILE_OPTIONS "")

target_link_libraries(yaml_oasim PUBLIC common)
