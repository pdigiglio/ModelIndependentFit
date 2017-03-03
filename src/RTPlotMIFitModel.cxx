/**
 *
 *    @file  RTPlotMIFitModel.cxx
 *   @brief  
 *
 *    @date  02/22/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "RTPlotMIFitModel.h"

#include <fwd/Model.h>

#include "KnownMassShapeFreedWave.h"
#include "MassRangePartition.h"
#include "ModelIndependentFitModel.h"
#include "RealTimeParameterPlot.h"

// -- Shouldn't be needed --
#include <TCanvas.h>
#include <TGraph.h>
// -------------------------

#include <cassert>
#include <complex>
#include <memory>
#include <string>
#include <vector>

std::vector<RealTimeParameterPlot> make_plots(const std::vector<std::shared_ptr<const KnownMassShapeFreedWave>>& fws) {
    std::vector<RealTimeParameterPlot> rtpp;
    rtpp.reserve(fws.size());

    std::transform(std::begin(fws), std::end(fws), std::back_inserter(rtpp),
                   [] (const auto& fw) { return RealTimeParameterPlot(binned_mass_shape(fw->massShape(), *fw->massRangePartition(), 1), fw->binName()); });

    return rtpp;
}

std::vector<std::shared_ptr<const FreedWave>> downcast_freed_waves(const std::vector<std::shared_ptr<const KnownMassShapeFreedWave>>& fws) {

    std::vector<std::shared_ptr<const FreedWave>> f;
    f.reserve(fws.size());
    std::transform(std::begin(fws), std::end(fws), std::back_inserter(f),
                   [] (const auto& fw) { return std::static_pointer_cast<const FreedWave>(fw); });
    return f;
}

RTPlotMIFitModel::RTPlotMIFitModel(std::unique_ptr<yap::Model> model,
                                   const std::vector<std::shared_ptr<const KnownMassShapeFreedWave>>& freed_waves,
                                   const std::string model_name) :
    ModelIndependentFitModel(std::move(model), downcast_freed_waves(freed_waves), model_name),
    ParameterPlots_(std::move(make_plots(freed_waves)))
{}

RTPlotMIFitModel::~RTPlotMIFitModel() = default;

void RTPlotMIFitModel::setParameters(const std::vector<std::vector<std::complex<double>>>& p) noexcept {
    ModelIndependentFitModel::setParameters(p);

    assert(p.size() == ParameterPlots_.size());
    for (size_t i = 0; i < ParameterPlots_.size(); ++ i)
        ParameterPlots_[i].setCurrentParameters(p[i]);
}
