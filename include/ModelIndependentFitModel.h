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

#include "fwd/MassRangePartition.h"
#include "fwd/FreedWave.h"

#include <fwd/Model.h>

#include "FitModel.h"

#include <complex>
#include <functional>
#include <memory>
#include <string>
#include <vector>

/// @ingroup Models
/// @brief Class to handle a `yap::Model` that will be used to perform a model-independent fit.
/// @attention The bin sorting relies on the mass-bin functions being named `Bin(i)`, being `i` an integer.
class ModelIndependentFitModel : public FitModel {
public:

    /// @brief Constructor.
    /// @param model          The YAP model of the decay.
    /// @param freed_waves    The freed waves in the model.
    /// @param model_name     The (optional) name of the model.
    explicit ModelIndependentFitModel(std::unique_ptr<yap::Model> model,
                                      const std::vector<std::shared_ptr<const FreedWave>>& freed_waves,
                                      const std::string model_name = "");

    /// _Default_ move constructor.
    ModelIndependentFitModel(ModelIndependentFitModel&&) noexcept            = default;
    /// _Default_ move assignment operator.
    ModelIndependentFitModel& operator=(ModelIndependentFitModel&&) noexcept = default;

    /// _Default_ destructor.
    virtual ~ModelIndependentFitModel();

    /// Return the freed-wave vector.
    const std::vector<std::shared_ptr<const FreedWave>>& freedWaves() const noexcept
    { return FreedWaves_; }

private:
    /// The freed waves in the model.
    std::vector<std::shared_ptr<const FreedWave>> FreedWaves_;

    /// Sort the bin free amplitudes consistently with the freed waves.
    void sortBinFreeAmplitudes();
};

/// @brief Check if the FreeAmplitude's are sored like their associated MassBin's.
/// @param m The model-independent fit model to check.
const bool mass_bin_sorted(const ModelIndependentFitModel& m);

#endif
