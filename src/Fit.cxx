// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "Fit.h"

#include "model/model_independent_d3pi.h"

#include "CachedIntegrator.h"
#include "FitModel.h"
#include "ModelIndependentFitModel.h"
#include "OnTheFlyIntegrator.h"
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

Fit::Fit(std::shared_ptr<FitModel> fit_model,
         std::unique_ptr<RootFitData> fit_data,
         std::unique_ptr<FitIntegrator> integrator):
    FitModel_(fit_model),
    Data_(std::move(fit_data)),
    Integrator_(std::move(integrator))
{
    // Loop over the non-fixed free amplitudes and add them to the fit parameters.
    for (const auto& fa : FitModel_->freeAmplitudes()) {
        const auto fa_name = free_amplitude_name(*fa);
#ifndef NDEBUG
        std::cout << "Add parameter for decay mode " << fa_name << std::endl;
#endif

        const double upper_amp_range = 20. * abs(fa->value());
        AddParameter("amp(" + fa_name + ")", 0., upper_amp_range);
        GetParameters().Back().SetPriorConstant();

        constexpr double upper_phase_range = 180.;
        AddParameter("d_phi(" + fa_name + ")", 0., upper_phase_range);
        GetParameters().Back().SetPriorConstant();
        
        // Add observables to get the real and imaginary parts of the amplitudes
        const double obs_range = upper_amp_range;
        AddObservable("Re(" + fa_name + ")", -obs_range, +obs_range);
        AddObservable("Im(" + fa_name + ")", -obs_range, +obs_range);
    }
}

const std::string Fit::logFileName() const noexcept
{ return Data_->path() + FitModel_->name() + "_log.txt"; }

const std::unique_ptr<yap::Model>& Fit::model() const noexcept
{ return FitModel_->model(); }

Fit::~Fit() = default;

// ---------------------------------------------------------
double Fit::LogLikelihood(const std::vector<double>& pars)
{
    // Sets the new values for the model free amplitudes and
    // evaluates the model integral matrix.
    setParameters(pars);

    // Sums up the components of the integral matrix and takes its log.
    const auto model_integral = log(integral(Integrator_->modelIntegral()).value());
    assert(std::isfinite(model_integral));

    // Evaluates the log-likelihood.
    const auto L = yap::sum_of_log_intensity(*FitModel_->model(), Data_->partitions(), model_integral);
#ifndef NDEBUG
    std::cout << "L = " << L << " log(I) = " << model_integral << std::endl;
#endif
    FitModel_->model()->setParameterFlagsToUnchanged();

    return L;
}

// ---------------------------------------------------------
void Fit::setParameters(const std::vector<double>& p) {
    assert(p.size() == 2 * FitModel_->freeAmplitudes().size());

    // Update the free-amplitude values.
    double cumulative_phase = 0.;
    for (auto& fa : FitModel_->freeAmplitudes()) {
        const auto amplitude = p[amplitude_parameter_index(fa, FitModel_)];

        // Add the phase difference to the cumulative phase.
        cumulative_phase += p[phase_parameter_index(fa, FitModel_)];
        *fa = std::polar<double>(amplitude, yap::rad<double>(cumulative_phase));
    }

    // Evaluate the integral of the model (with the new parameters).
    Integrator_->integrate();
}

// ---------------------------------------------------------
// double Fit::LogAPrioriProbability(const std::vector<double>& pars)
// {
//     // return the log of the prior probability p(pars)
//     // If you use built-in priors, leave this function commented out.
// }

// ---------------------------------------------------------
 void Fit::CalculateObservables(const std::vector<double>& p)
 {
     for (const auto& fa : FitModel_->freeAmplitudes()) {
         // Set the real value.
         const auto r = real_observable_index(fa, FitModel_);
         GetObservable(r) = std::real(fa->value());

         // Set the imaginary value.
         const auto i = imag_observable_index(fa, FitModel_);
         GetObservable(i) = std::imag(fa->value());
     }

#ifndef NDEBUG
     // Make sure that the free amplitudes correspond to the parametes.
     double cumulative_phase = 0.;
     for (size_t i = 0; i < FitModel_->freeAmplitudes().size(); ++ i) {
         const auto A = std::polar<double>(p[2 * i], yap::rad<double>(p[2 * i + 1] + cumulative_phase));
         cumulative_phase += p[2 * i + 1];

         assert(FitModel_->freeAmplitudes()[i]->value() == A);
     }
#endif
}

// ---------------------------------------------------------
void Fit::fixAmplitude(const std::shared_ptr<const yap::FreeAmplitude>& fa, const double value) {
    const auto i = amplitude_parameter_index(fa, FitModel_);
    GetParameter(i).Fix(value);
}

// ---------------------------------------------------------
void Fit::fixPhase(const std::shared_ptr<const yap::FreeAmplitude>& fa, const double value) {
    const auto i = phase_parameter_index(fa, FitModel_);
    GetParameter(i).Fix(value);
}

// ---------------------------------------------------------
void Fit::setAmplitudeRange(const std::shared_ptr<const yap::FreeAmplitude>& fa, double low, double high) {
    const auto i = amplitude_parameter_index(fa, FitModel_);
    GetParameter(i).SetLimits(low, high);
}

// ---------------------------------------------------------
void Fit::setPhaseRange(const std::shared_ptr<const yap::FreeAmplitude>& fa, double low, double high) {
    const auto i = phase_parameter_index(fa, FitModel_);
    GetParameter(i).SetLimits(low, high);
}

// ---------------------------------------------------------
//template <typename Integrator,
//          typename = std::enable_if_t<std::is_base_of<FitIntegrator, Integrator>::value>>
std::unique_ptr<Fit> create_fit(const char* file_path, const char* file_name, const char* model_name) {
    using Integrator = CachedIntegrator;

    // Create a model to fit the data.
    const std::shared_ptr<FitModel> fit_model(binned_d3pi(model_name));
//    assert(fit_model->freeAmplitudes().size() == 22);

    // Create the integrator.
    constexpr unsigned integration_points = 2e4;
    auto integrator(std::make_unique<Integrator>(std::static_pointer_cast<const FitModel>(fit_model), integration_points));

    auto root_file_handler(std::make_unique<RootFileHandler>(file_path, file_name));
    auto root_fit_data(std::make_unique<RootFitData>(std::move(root_file_handler), fit_model));

    // Create the BAT model for fitting the data.
    auto fit(std::make_unique<Fit>(fit_model, std::move(root_fit_data), std::move(integrator)));

//    // Fix amplitudes in the fit
//    {
//        // Get the non-fixed free amplitudes.
//        const auto fas = fit_model->freeAmplitudes();
//        // Iterator to the last (valid) one.
//        const auto fa  = std::prev(std::end(fas), 1);
//
//        // Fix the last amplitude.
//        fit->fixAmplitude(*fa, 1);
//        // Fix the first phase.
//        fit->fixPhase(*begin(fas), 0);
//    }

    return fit;
}
