// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "Fit.h"

#include "parameter_handling.h"
#include "model/model_independent_d3pi.h" // to create the fit

#include "CachedIntegrator.h"
#include "FitModel.h"
#include "FittedFreeAmplitude.h"
#include "FreedWave.h"
#include "ModelIndependentFitModel.h"
#include "OnTheFlyIntegrator.h"
#include "RTPlotMIFitModel.h"
#include "RootFileHandler.h"
#include "RootFitData.h"

#include <DecayChannel.h>
#include <FreeAmplitude.h>
#include <Spin.h>
#include <SpinAmplitude.h>

#include <cassert>
#include <iterator>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

// #include <BAT/BCMath.h>

// Helper function to return the index of the parameter corresponding to the amplitude 
// of the queried free amplitude.
inline const size_t amplitude_parameter_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                              const std::shared_ptr<const FitModel>& fit_model) {
    return 2 * free_amplitude_index(fa, fit_model);
}

// Helper function to return the index of the parameter corresponding to the phase
// of the queried free amplitude.
inline const size_t phase_parameter_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                          const std::shared_ptr<const FitModel>& fit_model) {
    return 2 * free_amplitude_index(fa, fit_model) + 1;
}

// Helper function to return the index of the observable corresponding to the real part
// of the queried free amplitude.
inline const size_t real_observable_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                          const std::shared_ptr<const FitModel>& fit_model) {
    return 2 * free_amplitude_index(fa, fit_model);
}

// Helper function to return the index of the observable corresponding to the imaginary part
// of the queried free amplitude.
inline const size_t imag_observable_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                          const std::shared_ptr<const FitModel>& fit_model) {
    return 2 * free_amplitude_index(fa, fit_model) + 1;
}

// Helper function to get a name for a FreeAmplitude
const std::string free_amplitude_name(const yap::FreeAmplitude& fa) noexcept {
    return yap::to_string(*fa.decayChannel())
        + " L = " + std::to_string(fa.spinAmplitude()->L())
        + " S = " + yap::spin_to_string(fa.spinAmplitude()->twoS());
}

template <typename FMT>
Fit<FMT>::Fit(std::shared_ptr<FMT> fit_model,
              std::unique_ptr<RootFitData> fit_data,
              std::unique_ptr<FitIntegrator> integrator):
    FitModel_(fit_model),
    Data_(std::move(fit_data)),
    Integrator_(std::move(integrator))
{
    // Compile-time check for type compatibility.
    static_assert(std::is_base_of<FitModel, FMT>::value,
                  "FitModelType doesn't inherit from FitModel");

    // Loop over the non-fixed free amplitudes and add them to the fit parameters.
    for (const auto& fa : freeAmplitudes()) {
        const auto fa_name = free_amplitude_name(*fa);
#ifndef NDEBUG
        std::cout << "Add parameter for decay mode " << fa_name << std::endl;
#endif

        const double upper_amp_range = 20. * abs(fa->value());
        AddParameter("amp(" + fa_name + ")", 0., upper_amp_range);
        GetParameters().Back().SetPriorConstant();

        constexpr double upper_phase_range = 180.;
        AddParameter("d_phi(" + fa_name + ")", 0, upper_phase_range);
        GetParameters().Back().SetPriorConstant();
        
        // Add observables to get the real and imaginary parts of the amplitudes
        const double obs_range = upper_amp_range;
        AddObservable("Re(" + fa_name + ")", -obs_range, +obs_range);
        AddObservable("Im(" + fa_name + ")", -obs_range, +obs_range);
    }
}

template <typename FMT>
Fit<FMT>::~Fit() = default;

// ---------------------------------------------------------
template <typename FMT>
double Fit<FMT>::LogLikelihood(const std::vector<double>& pars)
{

    // Non-const access to the fit model.
    FitModel_->setParameters(fit_to_yap_parameters(pars, fitModel()));

    // Evaluate the integral of the model (with the new parameters).
    Integrator_->integrate();

    // Sums up the components of the integral matrix and takes its log.
    const auto model_integral = log(integral(Integrator_->modelIntegral()).value());
    assert(std::isfinite(model_integral));

    // Evaluates the log-likelihood.
    const auto L = yap::sum_of_log_intensity(*fitModel()->model(), Data_->partitions(), model_integral);
#ifndef NDEBUG
    std::cout << "L = " << L << " log(I) = " << model_integral << std::endl;
#endif
    fitModel()->model()->setParameterFlagsToUnchanged();

    return L;
}

// ---------------------------------------------------------
// double Fit::LogAPrioriProbability(const std::vector<double>& pars)
// {
//     // return the log of the prior probability p(pars)
//     // If you use built-in priors, leave this function commented out.
// }

// ---------------------------------------------------------
template <typename FMT>
void Fit<FMT>::CalculateObservables(const std::vector<double>& p)
{
    std::for_each(std::begin(freeAmplitudes()), std::end(freeAmplitudes()),
        [&](const auto& fa) {
            this->GetObservable(real_observable_index(fa, this->fitModel())) = std::real(fa->value());
            this->GetObservable(imag_observable_index(fa, this->fitModel())) = std::imag(fa->value());
        });
}

// ---------------------------------------------------------
template <typename FMT>
void Fit<FMT>::fixAmplitude(const std::shared_ptr<const yap::FreeAmplitude>& fa, const double value) {
    const auto i = amplitude_parameter_index(fa, fitModel());
    GetParameter(i).Fix(value);
}

// ---------------------------------------------------------
template <typename FMT>
void Fit<FMT>::fixPhase(const std::shared_ptr<const yap::FreeAmplitude>& fa, const double value) {
    const auto i = phase_parameter_index(fa, fitModel());
    GetParameter(i).Fix(value);
}

// ---------------------------------------------------------
template <typename FMT>
void Fit<FMT>::setAmplitudeRange(const std::shared_ptr<const yap::FreeAmplitude>& fa, double low, double high) {
    const auto i = amplitude_parameter_index(fa, fitModel());
    GetParameter(i).SetLimits(low, high);
}

// ---------------------------------------------------------
template <typename FMT>
void Fit<FMT>::setPhaseRange(const std::shared_ptr<const yap::FreeAmplitude>& fa, double low, double high) {
    const auto i = phase_parameter_index(fa, fitModel());
    GetParameter(i).SetLimits(low, high);
}

// ---------------------------------------------------------
template <typename FMT>
const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& Fit<FMT>::freeAmplitudes() const noexcept {
    return fitModel()->freeAmplitudes();
}

// Explicitly instantiate the classes
//template class Fit<FitModel>;
//template class Fit<ModelIndependentFitModel>;
template class Fit<RTPlotMIFitModel>;
