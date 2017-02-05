# Get the directory that contains this directory
get_filename_component(PARENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

include_directories(${BAT_INCLUDE_DIRS}
                    ${ROOT_INCLUDE_DIR}
                    ${YAP_SOURCE_DIR}/include
                    ${PARENT_SOURCE_DIR}/include)

set(sources
    Fit.cxx
    FitIntegrator.cxx
    FitModel.cxx
    MassBin.cxx
    RootFileHandler.cxx
    RootFitData.cxx
    parallelism.cxx)

# Make output directory for the generated data
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)

set(programs runmiFit runMIGen)
foreach(program ${programs})
    add_executable(${program} EXCLUDE_FROM_ALL ${program}.cxx ${sources})
    target_link_libraries(${program} ${ROOT_LIBS} ${BAT_LIBRARIES} YAP)
endforeach(program)

