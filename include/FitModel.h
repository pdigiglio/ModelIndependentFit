/**
 *
 *    @file  FitModel.h
 *   @brief  
 *
 *    @date  02/03/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  FIT_MODEL_H
#define  FIT_MODEL_H

#include <fwd/MassRange.h>
#include <fwd/FreeAmplitude.h>

#include "fwd/Fit.h"

#include <MassAxes.h>
#include <Model.h>

#include <memory>
#include <vector>

class FitModel {
public:
    /// @brief Constructor.
    /// @param model      The YAP model of the decay.
    /// @param model_name The (optional) name of the model.
    explicit FitModel(std::unique_ptr<yap::Model> model,
                      const std::string model_name = "");

    /// _Default_ destructor.
    ~FitModel();

    /// Returns the model name.
    const std::string& name() const noexcept
    { return Name_; }

    /// Returns the YAP model.
    const std::unique_ptr<yap::Model>& model() const noexcept
    { return Model_; }

    /// Return the default mass axes.
    const yap::MassAxes& massAxes() const noexcept
    { return MassAxes_; }

    /// Returns the mass ranges.
    const std::vector<yap::MassRange>& massRanges() const noexcept
    { return MassRanges_; }

    /// Access the non-fixed free amplitudes.
    yap::FreeAmplitudeVector& freeAmplitudes() noexcept
    { return FreeAmplitudes_; }

    /// Use a common radial size [GeV^-1] for all the resonances.
    static constexpr const double radialSize() noexcept
    { return 3.; }

    /// Returns the mass of the D+.
    static constexpr const double Dmass() noexcept
    { return 1.8691; }

private:

    /// The (optional) name of the model.
    const std::string Name_;

    /// The decay fit model.
    const std::unique_ptr<yap::Model> Model_;

    /// The default mass axes.
    const yap::MassAxes MassAxes_;

    /// The mass ranges for the decay.
    const std::vector<yap::MassRange> MassRanges_;

    /// @brief The non-fixed FreeAmplitude's of the model.
    /// @attention By chenging these, the values of the free
    /// amplitudes in the model will also change.
    yap::FreeAmplitudeVector FreeAmplitudes_;
};

/// @brief Helper function to create a YAP model.
/// @tparam T The formalism to use for the model.
template <typename T>
inline std::unique_ptr<yap::Model> make_model()
{ return std::make_unique<yap::Model>(std::make_unique<T>()); }

/// @brief Creates a new FitModel instance.
/// @param model_name The name of the model.
std::unique_ptr<FitModel> make_fit_model(const std::string& model_name);

/// Attempts an initial guess on the model parameters.
const std::vector<double> guess_parameters(const Fit& m);

/// @brief Prints the fit parameters in an output text file.
/// @param m The Fit model whose parameters to print.
void write_fit_result_to_file(Fit& m);

#endif
