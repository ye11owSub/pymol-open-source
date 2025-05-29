cmake_minimum_required(VERSION 3.15)

project(${TARGET_NAME})

set(CMAKE_VERBOSE_MAKEFILE on)
set(Python_FIND_VIRTUALENV FIRST)
set(Python_FIND_STRATEGY LOCATION)

find_package(Python REQUIRED COMPONENTS Interpreter Development NumPy)

add_library(${TARGET_NAME} SHARED ${ALL_SRC})

foreach(FOLDER_NAME IN ITEMS layer0 layer1 layer2 layer3 layer4 layer5 ov)
    include(${CMAKE_CURRENT_SOURCE_DIR}/${FOLDER_NAME}/sources.cmake)
endforeach()

target_link_options(${TARGET_NAME} PRIVATE
    "$<$<PLATFORM_ID:Darwin>:LINKER:-undefined,dynamic_lookup>"
)

target_compile_options(${TARGET_NAME} PRIVATE
    ${ALL_COMP_ARGS}

    $<$<NOT:$<CXX_COMPILER_ID:MSVC>>:
        -Werror=return-type
        -Wunused-variable
        -Wno-switch
        -Wno-narrowing
        -Wno-char-subscripts
        $<$<CONFIG:Debug>:-Og>
        $<$<NOT:$<CONFIG:Debug>>:-O3>
    >
    $<$<CXX_COMPILER_ID:MSVC>:/MP>
    $<$<CXX_COMPILER_ID:MSVC>:/std:c++17>
    $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:MSVC>>:/Z7>

    # optimization currently causes a clang segfault on OS X 10.9 when
    # compiling layer2/RepCylBond.cpp
    $<$<PLATFORM_ID:Darwin>:-fno-strict-aliasing>
)

execute_process(
    COMMAND ${Python_EXECUTABLE}  -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
    OUTPUT_VARIABLE PYTHON_MODULE_SUFFIX
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

function(set_output_path lib path)
    set_target_properties(${lib} PROPERTIES
	PREFIX ""
	SUFFIX "${PYTHON_MODULE_SUFFIX}"
	OUTPUT_NAME "${lib}"
        LIBRARY_OUTPUT_DIRECTORY "${path}"
        RUNTIME_OUTPUT_DIRECTORY "${path}"
    )

    foreach(CONFIG IN ITEMS RELEASE DEBUG RELWITHDEBINFO MINSIZEREL)
        set_target_properties(${lib} PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY_${CONFIG} "${path}"
            RUNTIME_OUTPUT_DIRECTORY_${CONFIG} "${path}"
        )
    endforeach()
endfunction()

set_output_path(${TARGET_NAME} "${Python_SITEARCH}/pymol")

target_compile_features(${TARGET_NAME} PRIVATE cxx_std_17)

target_include_directories(${TARGET_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/contrib/pocketfft
    ${ALL_INC_DIR}
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
)

target_link_directories(${TARGET_NAME} PUBLIC
    ${ALL_LIB_DIR}
    $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARY_DIRS}>
)
target_link_libraries(${TARGET_NAME}
    ${ALL_LIB}
    ${ALL_EXT_LINK}
    # pyconfig.py forces linking against pythonXY.lib on MSVC
    $<$<CXX_COMPILER_ID:MSVC>:${Python_LIBRARIES}>
    $<$<AND:$<CONFIG:Debug>,$<CXX_COMPILER_ID:MSVC>>:/DEBUG>
)

target_compile_definitions(${TARGET_NAME} PUBLIC
    ${ALL_DEF}
    _PYMOL_LIBPNG
    _PYMOL_FREETYPE
    _PYMOL_NUMPY
    $<$<CXX_COMPILER_ID:MSVC>:WIN32>
    $<$<PLATFORM_ID:Darwin>:PYMOL_CURVE_VALIDATE>

    # bounds checking in STL containers
    # if DEBUG and not Windows compiler
    $<$<AND:$<CONFIG:Debug>,$<NOT:$<CXX_COMPILER_ID:MSVC>>>:_GLIBCXX_ASSERTIONS>
)

add_subdirectory(contrib/champ)
