/**
 *
 *    @file  RealTimePlotModelIndependentFitModel.h
 *   @brief  
 *
 *    @date  02/21/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  REAL_TIME_PLOT_MODEL_INDEPENDENT_FIT_MODEL_H
#define  REAL_TIME_PLOT_MODEL_INDEPENDENT_FIT_MODEL_H

#include "fwd/MassRangePartition.h"
#include "fwd/RealTimeParameterPlot.h"

#include "ModelIndependentFitModel.h"

#include <complex>
#include <memory>

/// @ingroup Models
/// @brief Class to handle a `yap::Model` for a model-independent fit while showing a real-time plot of the fit parameters (useful for debugging).
class RealTimePlotModelIndependentFitModel final : public ModelIndependentFitModel {
public:

    /// @brief Constructor.
    /// @param model              The YAP model of the decay.
    /// @param mass_partition     The Dalitz-plot mass partition.
    /// @param real_time_par_plot The tool to visualize the fit in real time.
    /// @param model_name         The (optional) name of the model.
    explicit RealTimePlotModelIndependentFitModel(std::unique_ptr<yap::Model> model,
                                                  std::unique_ptr<const MassRangePartition> mass_partition,
                                                  std::unique_ptr<RealTimeParameterPlot> real_time_par_plot,
                                                  const std::string model_name = "");

    /// @brief Constuctor from a ModelIndependentFitModel and a RealTimeParameterPlot.
    /// @param mimf               The model to perform the model-independent fit.
    /// @param real_time_par_plot The visualization tool.
    explicit RealTimePlotModelIndependentFitModel(ModelIndependentFitModel mimf,
                                                  std::unique_ptr<RealTimeParameterPlot> real_time_par_plot);

    /// @brief _Default_ destructor.
    ~RealTimePlotModelIndependentFitModel() final;

    /// @brief Sets the FreeAmplitude's in the model and plots them.
    /// @param p The BAT parameters to convert into YAP parameters.
    void setParameters(const std::vector<std::complex<double>>& p) noexcept final;

private:
    /// A ROOT canvas that will be updated every time the model parameters change.
    std::unique_ptr<RealTimeParameterPlot> ParameterPlot_;
};

#endif