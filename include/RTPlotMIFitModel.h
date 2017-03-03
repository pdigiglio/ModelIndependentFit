/**
 *
 *    @file  RTPlotMIFitModel.h
 *   @brief  
 *
 *    @date  02/21/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  REAL_TIME_PLOT_MODEL_INDEPENDENT_FIT_MODEL_H
#define  REAL_TIME_PLOT_MODEL_INDEPENDENT_FIT_MODEL_H

#include "fwd/KnownMassShapeFreedWave.h"
#include "fwd/RealTimeParameterPlot.h"

#include <fwd/Model.h>

#include "ModelIndependentFitModel.h"

#include <complex>
#include <memory>
#include <vector>

/// @ingroup Models
/// @brief   Real-time-plot model-independent fit model.
/// @details Handle a `yap::Model` for a model-independent fit while showing a real-time plot
///          of the fit parameters (useful for debugging).
class RTPlotMIFitModel final : public ModelIndependentFitModel {
public:

    /// @brief Constructor.
    /// @param model       The YAP model of the decay.
    /// @param freed_waves The freed waves (with known mass shapes) in the model.
    /// @param model_name  The (optional) name of the model.
    explicit RTPlotMIFitModel(std::unique_ptr<yap::Model> model,
                              const std::vector<std::shared_ptr<const KnownMassShapeFreedWave>>& freed_waves,
                              const std::string model_name = "");

    /// @brief _Default_ destructor.
    ~RTPlotMIFitModel() final;

    /// @brief Sets the FreeAmplitude's in the model and plots them.
    /// @param p The BAT parameters to convert into YAP parameters.
    void setParameters(const std::vector<std::vector<std::complex<double>>>& p) noexcept final;

private:

    /// The vector of plots of the fit parameters (one entry per freed wave).
    std::vector<RealTimeParameterPlot> ParameterPlots_;
};

#endif
