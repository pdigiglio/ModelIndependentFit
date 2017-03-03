// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__MIFIT__H
#define __BAT__MIFIT__H

#include "fwd/FitIntegrator.h"
#include "fwd/FitModel.h"
#include "fwd/RootFitData.h"

#include <fwd/FreeAmplitude.h>
#include <fwd/Model.h>

#include <BAT/BCModel.h>

#include <memory>

/// @brief The class that will perform the fit.
/// @tparam FitModelType The type of the FitModel used to perform the fit.
template <typename FitModelType = FitModel>
class Fit final : public BCModel {
public:
    /// @brief Constructor.
    /// @param fit_model  The decay model.
    /// @param fit_data   The data to fit the model to.
    /// @param integrator The model integrator.
    explicit Fit(std::shared_ptr<FitModelType> fit_model,
                 std::unique_ptr<RootFitData> fit_data,
                 std::unique_ptr<FitIntegrator> integrator);

    /// _Default_ destructor.
    ~Fit();

    // Overload LogLikelihood to implement model
    double LogLikelihood(const std::vector<double>& pars) final;

    // Overload LogAprioriProbability if not using built-in 1D priors
    // double LogAPrioriProbability(const std::vector<double> & pars);

    /// @brief Calculate real and imaginary parts of the FreeAmplitude's.
    /// @param pars The vector of the FreeAmplitude's being fitted.
    void CalculateObservables(const std::vector<double> & pars) final;

    /// @brief Fix the amplitude of the free amplitude _fa_.
    /// @param fa    The FreeAmplitude whose component one wants to fix.
    /// @param value The value to fix the FreeAmplitude's amplitude to.
    void fixAmplitude(const std::shared_ptr<const yap::FreeAmplitude>& fa, const double value);

    /// @brief Fix the phase of the free amplitude _fa_.
    /// @param fa    The FreeAmplitude whose component one wants to fix.
    /// @param value The value to fix the FreeAmplitude's phase to.
    void fixPhase(const std::shared_ptr<const yap::FreeAmplitude>& fa, const double value);

    /// Return the fit data.
    const std::unique_ptr<RootFitData>& rootFitData() const noexcept
    { return Data_; }

    /// Returns the fit model.
    const std::shared_ptr<const FitModelType> fitModel() const noexcept
    { return FitModel_; }

    /// @brief Sets the range of the parameter corresponding to the amplitude of _fa_.
    /// @param fa   The free amplitude whose parameter value one wants to fix.
    /// @param low  The lower boundary of the range.
    /// @param high The upper boundary of the range.
    void setAmplitudeRange(const std::shared_ptr<const yap::FreeAmplitude>& fa, double low, double high);

    /// @brief Sets the range of the parameter corresponding to the phase of _fa_.
    /// @param fa   The free amplitude whose parameter value one wants to fix.
    /// @param low  The lower boundary of the range.
    /// @param high The upper boundary of the range.
    void setPhaseRange(const std::shared_ptr<const yap::FreeAmplitude>& fa, double low, double high);

private:

    /// The fit model for the decay.
    const std::shared_ptr<FitModelType> FitModel_;

    /// The data to fit the model to.
    const std::unique_ptr<RootFitData> Data_;

    /// The model integrator.
    const std::unique_ptr<FitIntegrator> Integrator_;

    /// Helper function to get the FitModel free amplitudes.
    const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& freeAmplitudes() const noexcept;
};

#endif
