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

/// Wrapper for a `yap::Model` with some facility member function to make the model-indepentent fit easier.
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
    const std::vector<yap::MassRange>& squaredMassRanges() const noexcept
    { return SquaredMassRanges_; }

    /// Access the non-fixed free amplitudes.
    yap::FreeAmplitudeVector& freeAmplitudes() noexcept
    { return FreeAmplitudes_; }

    /// Return the non-fixed free amplitudes.
    const yap::FreeAmplitudeVector& freeAmplitudes() const noexcept
    { return FreeAmplitudes_; }

    /// Use a common radial size [GeV^-1] for all the resonances.
    static constexpr const double radialSize() noexcept
    { return 3.; }

    /// Returns the mass of the D+.
    static constexpr const double Dmass() noexcept
    { return 1.86961; }

private:

    /// The (optional) name of the model.
    const std::string Name_;

    /// The decay fit model.
    const std::unique_ptr<yap::Model> Model_;

    /// The default mass axes.
    const yap::MassAxes MassAxes_;

    /// The squared-mass ranges for the decay.
    const std::vector<yap::MassRange> SquaredMassRanges_;

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

/// @brief Helper function to set up the model for generating MC data.
std::unique_ptr<yap::Model> d3pi();

/// @brief Creates a new FitModel instance.
/// @param model_name The name of the model.
std::unique_ptr<FitModel> make_fit_model(const std::string& model_name);

/// @brief Attempts an initial guess for the fit parameters based on the model used to generate MC data.
/// @param m (What is it for?)
const std::vector<double> guess_parameters(Fit& m);

/// Returns the index of a free amplitude in the FreeAmplitude vector of the fit model.
/// @param fa        The free amplitude to query for.
/// @param fit_model The model in which to look for the free amplitude.
/// @throw yap::exceptions::Exception If the free queried free amplitude is not in the vector.
const size_t free_amplitude_index(const std::shared_ptr<const yap::FreeAmplitude>& fa, const std::shared_ptr<const FitModel>& fit_model);

/// @brief Prints the fit parameters in an output text file.
/// @param m The Fit model whose parameters to print.
void write_fit_result_to_file(Fit& m);

#endif
