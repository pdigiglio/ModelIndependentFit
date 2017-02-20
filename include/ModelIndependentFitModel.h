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

#include "fwd/MassBin.h"
#include <fwd/Model.h>

#include "FitModel.h"

#include <string>
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
                                      const std::string model_name = "");

    /// Return the Dalitz-plot mass partition.
    const std::vector<double>& massPartition() const noexcept
    { return MassPartition_; }

    /// Check if the MassBin's associated with the FreeAmplitude's are sorted like the mass partition.
    const bool massBinSorted() const noexcept;

    /// Sort the bin free amplitudes.
    void sortBinFreeAmplitudes();

private:
    /// The partition of the Dalitz-plot mass range.
    const std::vector<double> MassPartition_;
};

/// @brief Returns the mass shape of the bin corresponding to _fa_.
/// @param fa The free amplitude whose associated mass bin to access.
/// @param M  The model in which to look for the mass bin associated to _fa_.
/// @exception yap::Exception If the bin free amplitudes are not correctly sorted.
std::shared_ptr<const MassBin> bin_mass_shape(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                              const ModelIndependentFitModel& M);

#endif
