/**
 *
 *    @file  RealTimePlotModelIndependentFitModel.cxx
 *   @brief  
 *
 *    @date  02/22/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "RealTimePlotModelIndependentFitModel.h"

#include <fwd/Model.h>

#include "MassRangePartition.h"
#include "ModelIndependentFitModel.h"
#include "RealTimeParameterPlot.h"

#include <cassert>
#include <complex>
#include <memory>
#include <string>
#include <vector>

RealTimePlotModelIndependentFitModel::RealTimePlotModelIndependentFitModel(std::unique_ptr<yap::Model> model,
                                                                           std::unique_ptr<const MassRangePartition> mass_partition,
                                                                           std::unique_ptr<RealTimeParameterPlot> real_time_par_plot,
                                                                           const std::string model_name) :
    ModelIndependentFitModel(std::move(model), std::move(mass_partition), model_name),
    ParameterPlot_(std::move(real_time_par_plot))
{}

RealTimePlotModelIndependentFitModel::RealTimePlotModelIndependentFitModel(ModelIndependentFitModel mimf,
        std::unique_ptr<RealTimeParameterPlot> real_time_par_plot) :
    ModelIndependentFitModel(std::move(mimf)),
    ParameterPlot_(std::move(real_time_par_plot))
{}

RealTimePlotModelIndependentFitModel::~RealTimePlotModelIndependentFitModel() = default;

void RealTimePlotModelIndependentFitModel::setParameters(const std::vector<std::complex<double>>& p) noexcept {
    ModelIndependentFitModel::setParameters(p);

//    assert(false);
    // Plot the parameters.
    ParameterPlot_->setCurrentParameters(p);
}
