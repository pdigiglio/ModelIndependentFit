/**
 *
 *    @file  ModelIndependentFitModel.h
 *   @brief  
 *
 *    @date  02/14/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  MODEL_INDEPENDENT_FIT_MODEL_H
#define  MODEL_INDEPENDENT_FIT_MODEL_H

#include "fwd/Fit.h"

#include <fwd/Model.h>

#include "FitModel.h"

#include <iosfwd>
#include <vector>

/// Class to handle a YAP model that will perform a model-independent fit.
class ModelIndependentFitModel : public FitModel {
public:

    /// @brief Constructor.
    /// @param model          The YAP model of the decay.
    /// @param mass_partition The Dalitz-plot mass partition.
    /// @param model_name     The (optional) name of the model.
    explicit ModelIndependentFitModel(std::unique_ptr<yap::Model> model,
                                      const std::vector<double>& mass_partition,
                                      const std::string model_name = "") :
        FitModel(std::move(model), model_name),
        MassPartition_(mass_partition)
    {}

    /// Return the Dalitz-plot mass partition.
    const std::vector<double>& massPartition() const noexcept
    { return MassPartition_; }

private:
    /// The partition of the Dalitz-plot mass range.
    const std::vector<double> MassPartition_;
};

#endif
