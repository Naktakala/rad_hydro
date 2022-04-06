cmake_minimum_required(VERSION 3.2)

#================================================ Change these
SET(EXE_TARGET rad_hydro)
SET(LIB_TARGET rad_hydro_lib)
SET(CHI_MODULE_DIR RadHydro)
SET(PROJECT_NAME RadHydro)

#------------------------------------------------ DEPENDENCIES
if (NOT DEFINED CHI_TECH_DIR)
    if (NOT (DEFINED ENV{CHI_TECH_DIR}))
        message(FATAL_ERROR "***** CHI_TECH_DIR is not set *****")
    else()
        set(CHI_TECH_DIR "$ENV{CHI_TECH_DIR}")
    endif()
endif()
message(STATUS "CHI_TECH_DIR set to ${CHI_TECH_DIR}")

#================================================ Set project targets
project(${PROJECT_NAME} C CXX)

#================================================ Set cmake variables
if (NOT DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
    set(CONFIG_FILE_NAME "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/config.cmake")
endif()

#================================================ Include macros
set(TARGET ${EXE_TARGET})
include(${CHI_TECH_DIR}/ChiResources/Macros/Downstream.cmake)
include(${CHI_TECH_DIR}/ChiResources/Macros/Filter.cmake)

#================================================ Default include directories
include_directories(${PROJECT_SOURCE_DIR}/${CHI_MODULE_DIR})

#================================================ Library directories
# Link to external dependencies other than ChiTech
# Example:
# link_directories("${LBSOLVER_DIR}/chi_build")

#================================================ Define source directories
# The directory added below will need to populate
# the variable ${SOURCES}
add_subdirectory(${PROJECT_SOURCE_DIR}/${CHI_MODULE_DIR})

#================================================ Set Misc flags
# set(CMAKE_CXX_FLAGS_DEBUG  "${CMAKE_CXX_FLAGS_DEBUG} -O3")
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -O3")

#================================================ Create executable
add_library(${LIB_TARGET} ${SOURCES})
add_executable(${EXE_TARGET} "main.cc")

target_link_libraries(${LIB_TARGET} ${CHI_LIBS})
target_link_libraries(${EXE_TARGET} ${LIB_TARGET})


# |------------ Write Makefile to root directory
file(WRITE ${PROJECT_SOURCE_DIR}/Makefile "subsystem:\n" "\t$(MAKE) -C chi_build \n\n"
        "clean:\n\t$(MAKE) -C chi_build clean\n")