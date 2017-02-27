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
#include "RealTimePlotModelIndependentFitModel.h"
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
    for (const auto& fa : freeAmplitudes()) {
        const auto fa_name = free_amplitude_name(*fa);
#ifndef NDEBUG
        std::cout << "Add parameter for decay mode " << fa_name << std::endl;
#endif

        const double upper_amp_range = 20. * abs(fa->value());
        AddParameter("amp(" + fa_name + ")", 0., upper_amp_range);
        GetParameters().Back().SetPriorConstant();

        constexpr double upper_phase_range = 180.;
        AddParameter("d_phi(" + fa_name + ")", -upper_phase_range, upper_phase_range);
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
    std::static_pointer_cast<RealTimePlotModelIndependentFitModel>(FitModel_)->setParameters(fit_to_yap_parameters(pars));
    // Evaluate the integral of the model (with the new parameters).
    Integrator_->integrate();

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
     {
         // Make sure that the free amplitudes correspond to the parametes.
         const auto yp = fit_to_yap_parameters(p);
         assert(std::all_of(std::begin(FitModel_->freeAmplitudes()), std::end(FitModel_->freeAmplitudes()),
                            [&](const auto& fa) { return fa->value() == yp[free_amplitude_index(fa, this->fitModel())]; }));
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
const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& Fit::freeAmplitudes() const noexcept {
    return fitModel()->freeAmplitudes();
}

// ---------------------------------------------------------
std::vector<std::complex<double>> fit_to_yap_parameters(const std::vector<double>& p) noexcept {
    // Check that the parameters are couple of real numbers.
    assert(p.size() % 2 == 0);
    
    // The vector of YAP parameters.
    std::vector<std::complex<double>> yap_params;
    yap_params.reserve(p.size() / 2);

    // Convert the parameters.
    double cumulative_phase = 0.;
    for (size_t i = 0; i < p.size() / 2; ++ i) {
        cumulative_phase += p[2 * i + 1];
        yap_params.emplace_back(std::polar<double>(p[2 * i], yap::rad(cumulative_phase)));
    }

    return yap_params;
}

// ---------------------------------------------------------
std::vector<double> yap_to_fit_parameters(const std::vector<std::complex<double>>& p) noexcept {
    std::vector<double> fit_params;
    fit_params.reserve(2 * p.size());

    // TODO don't hard-code this!!
    double cumulative_phase = yap::deg(std::arg(*std::begin(p)));
    for (const auto& yp : p) {
        fit_params.push_back(std::abs(yp));
        fit_params.push_back(yap::deg(std::arg(yp)) - cumulative_phase);

        cumulative_phase += fit_params.back();
    }

    return fit_params;
}

// ---------------------------------------------------------
std::vector<double> yap_to_fit_parameters(const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& fas) noexcept {
    std::vector<std::complex<double>> p;
    p.reserve(fas.size());
    std::transform(std::begin(fas), std::end(fas), std::back_inserter(p),
                   [] (const auto& fa) { return fa->value(); });

    return yap_to_fit_parameters(p);
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

    // Fix amplitudes in the fit
    {
        // Get the non-fixed free amplitudes.
        const auto fas = fit_model->freeAmplitudes();
//        // Iterator to the last (valid) one.
//        const auto fa  = std::prev(std::end(fas), 1);
//
//        // Fix the last amplitude.
//        fit->fixAmplitude(*fa, 1);

        const auto guess_params = yap_to_fit_parameters(fas);
        // Fix the first phase.
        fit->fixPhase(*begin(fas), guess_params[1]);
    }

    return fit;
}
