# Get the directory that contains this directory
get_filename_component(PARENT_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

include_directories(${BAT_INCLUDE_DIRS}
                    ${ROOT_INCLUDE_DIR}
                    ${YAP_SOURCE_DIR}/include
                    ${PARENT_SOURCE_DIR}/include)

set(common
    CachedIntegrator.cxx
    FitIntegrator.cxx
    OnTheFlyIntegrator.cxx
    FitModel.cxx
    ModelIndependentFitModel.cxx
    RealTimeParameterPlot.cxx
    RTPlotMIFitModel.cxx
    MassBin.cxx
    MassRangePartition.cxx
    parallelism.cxx
    parameter_handling.cxx
    FitResultFileIterator.cxx
    FittedFreeAmplitude.cxx
    FreedWave.cxx
    KnownMassShapeFreedWave.cxx)

set(fit Fit.cxx)
set(data_generator DataGenerator.cxx)

set(root
    RootFileHandler.cxx
    RootFitData.cxx
    TTreeHandler.cxx)


# Add the models
set(models
    model/model_independent_d3pi.cxx
    model/d3pi.cxx)

# Make output directory for the generated data
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/output)

# Make a link to the helper bash script to run the fitter.
add_custom_target(link-run-script
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${PARENT_SOURCE_DIR}/fit.sh fit.sh)

# Make a link to the script to plot the results.
add_custom_target(link-plot-script
    COMMAND ${CMAKE_COMMAND} -E create_symlink ${PARENT_SOURCE_DIR}/plot_result.plt
                                               ${CMAKE_CURRENT_BINARY_DIR}/output/plot_result.plt)

add_executable(runmiFit EXCLUDE_FROM_ALL runmiFit.cxx ${common} ${fit} ${root} ${models})
target_link_libraries(runmiFit ${ROOT_LIBARIES} ${BAT_LIBRARIES} YAP)

add_executable(runMIGen EXCLUDE_FROM_ALL runMIGen.cxx ${common} ${data_generator} ${models})
target_link_libraries(runMIGen ${ROOT_LIBARIES} ${BAT_LIBRARIES} YAP)

add_executable(root_to_txt EXCLUDE_FROM_ALL root_to_txt.cxx)
target_link_libraries(root_to_txt ${ROOT_LIBRARIES})

# Make sure the script is linked
add_dependencies(runmiFit link-run-script link-plot-script)
