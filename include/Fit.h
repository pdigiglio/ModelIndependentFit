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
#include <string>

/// @brief The class that will perform the fit.
class Fit final : public BCModel {
public:
    /// @brief Constructor.
    /// @param fit_model  The decay model.
    /// @param fit_data   The data to fit the model to.
    /// @param integrator The model integrator.
    explicit Fit(std::shared_ptr<FitModel> fit_model,
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

    /// Returns the name of the log file.
    const std::string logFileName() const noexcept;

    /// Access the YAP model.
    const std::unique_ptr<yap::Model>& model() const noexcept;

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
    const std::shared_ptr<FitModel> FitModel_;

    /// The data to fit the model to.
    const std::unique_ptr<RootFitData> Data_;

    /// The model integrator.
    const std::unique_ptr<FitIntegrator> Integrator_;

    /// @brief Sets the free amplitudes to the values of _p_.
    /// @details This is the function where one can interpret the parameters
    /// as they want through the mathematical function one uses to assign them
    /// to the free amplitudes.
    /// @param p The new value of the free amplitudes.
    void setParameters(const std::vector<double>& p);
};

/// @brief Helper function to create a Fit class.
/// @param file_path  The directory of the input ROOT file.
/// @param file_name  The name of the input ROOT file.
/// @param model_name The name of the model.
std::unique_ptr<Fit> create_fit(const char* file_path, const char* file_name, const char* model_name = "");

#endif
