cmake_minimum_required (VERSION 3.0.0)
project (cobyla)

add_library (cobyla cobyla.c) 

target_include_directories (cobyla PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}) 