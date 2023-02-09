cmake_minimum_required(VERSION 3.0)

add_library(yaml OBJECT yaml_types.F90 yaml_settings.F90 yaml.F90)
set_property(TARGET yaml PROPERTY POSITION_INDEPENDENT_CODE ON)
set_property(TARGET yaml PROPERTY COMPILE_OPTIONS "")

target_link_libraries(yaml PUBLIC common)