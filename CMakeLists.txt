cmake_minimum_required(VERSION 3.9)

project(FalaiseModules)

find_package(Falaise REQUIRED)

add_compile_options(-W -Wall -Werror -O2)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

root_generate_dictionary(betabeta-data LINKDEF dict/betabeta-data-linkdef.hh)
root_generate_dictionary(ttd-residual-data LINKDEF dict/ttd-residual-data-linkdef.hh)
root_generate_dictionary(udd-calorimeter-waveform-flagging-data LINKDEF dict/udd-calorimeter-waveform-flagging-data-linkdef.hh)

add_library(FalaiseModules SHARED

  bank-removal-module.cc

  betabeta-selection-module.cc
  betabeta-data

  cd-killer-module.cc
  ttd-residual-module.cc
  ttd-residual-data
)

target_link_libraries(FalaiseModules Falaise::FalaiseModule)

# for file in conf/*in ...
configure_file("conf/betabeta-selection-pipeline.conf.in" "betabeta-selection-pipeline.conf")
configure_file("conf/cd-killer-pipeline.conf.in" "cd-killer-pipeline.conf")
configure_file("conf/ttd-residual-pipeline.conf.in" "ttd-residual-pipeline.conf")
