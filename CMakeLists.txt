cmake_minimum_required(VERSION 3.9)

project(FalaiseModules)

find_package(Falaise REQUIRED)

add_library(FalaiseModules SHARED
  cd-killer-module.cc)

target_link_libraries(FalaiseModules Falaise::FalaiseModule)

configure_file("conf/cd-killer-pipeline.conf.in" "cd-killer-pipeline.conf")

