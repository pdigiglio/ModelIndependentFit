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
#include "fwd/MassBin.h"
#include <fwd/Model.h>

#include "FitModel.h"

#include <functional>
#include <memory>
#include <string>

/// @ingroup Models
/// @brief Class to handle a `yap::Model` that will be used to perform a model-independent fit.
/// @attention The bin sorting relies on the mass-bin functions being named `Bin(i)`, being `i` an integer.
class ModelIndependentFitModel : public FitModel {
public:

    /// @brief Constructor.
    /// @param model          The YAP model of the decay.
    /// @param mass_partition The Dalitz-plot mass partition.
    /// @param model_name     The (optional) name of the model.
    explicit ModelIndependentFitModel(std::unique_ptr<yap::Model> model,
                                      std::unique_ptr<const MassRangePartition> mass_partition,
                                      const std::string model_name = "");

    /// _Default_ move constructor.
    ModelIndependentFitModel(ModelIndependentFitModel&&) = default;
    /// _Default_ move assignment operator.
    ModelIndependentFitModel& operator=(ModelIndependentFitModel&&) = default;

    /// _Default_ destructor.
    virtual ~ModelIndependentFitModel();

    /// Return the Dalitz-plot mass partition.
    const std::unique_ptr<const MassRangePartition>& massRangePartition() const noexcept
    { return MassRangePartition_; }

    /// Return the low edges of the Dalitz-plot mass partition.
    const std::vector<double>& massPartition() const noexcept;

    /// Check if the MassBin's associated with the FreeAmplitude's are sorted like the mass partition.
    const bool massBinSorted() const;

    /// Sort the bin free amplitudes.
    void sortBinFreeAmplitudes();

private:
    /// The partition of the Dalitz-plot mass range.
    std::unique_ptr<const MassRangePartition> MassRangePartition_;
};

/// @brief Returns the mass shape of the bin corresponding to _fa_.
/// @param fa The free amplitude whose associated mass bin to access.
/// @param M  The model in which to look for the mass bin associated to _fa_.
/// @exception yap::Exception If the bin free amplitudes are not correctly sorted.
std::shared_ptr<const MassBin> bin_mass_shape(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                              const ModelIndependentFitModel& M);

/// @brief Return the value of the mass shape on the low edges of the partition.
/// @param mass_shape The mass shape one wants to evaluate.
/// @param fm         The model-independent fit model whose partitioning to use.
std::vector<std::complex<double>> binned_mass_shape(std::function<std::complex<double>(double)> mass_shape,
                                                    const ModelIndependentFitModel& fm);

#endif
