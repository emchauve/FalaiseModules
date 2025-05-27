cmake_minimum_required(VERSION 3.9)

project(FalaiseModules)

find_package(Falaise REQUIRED)

add_compile_options(-W -Wall -Werror -O2)

add_library(FalaiseModules SHARED
  cd-killer-module.cc
  ttd-residual-module.cc
)

target_link_libraries(FalaiseModules Falaise::FalaiseModule)

configure_file("conf/cd-killer-pipeline.conf.in" "cd-killer-pipeline.conf")
configure_file("conf/ttd-residual-pipeline.conf.in" "ttd-residual-pipeline.conf")
