// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "Fit.h"

#include "FitIntegrator.h"
#include "FitModel.h"
#include "RootFileHandler.h"
#include "RootFitData.h"

//#include "model-independent/d3pi_one_resonance.h"
#include "Exceptions.h"
#include "FreeAmplitude.h"
#include "Model.h"

#include <AmplitudeBasis.h>
#include <Attributes.h>
#include <BreitWigner.h>
#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayTree.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <ImportanceSampler.h>
#include <MassAxes.h>
#include <MassRange.h>
#include <MathUtilities.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <PDL.h>
#include <PHSP.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <SpinAmplitude.h>
#include <SpinAmplitudeCache.h>
#include <VariableStatus.h>
#include <ZemachFormalism.h>
#include <logging.h>
#include <BreitWigner.h>
#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <PDL.h>
#include <SpinAmplitude.h>
#include <VariableStatus.h>

#include <cassert>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>

// #include <BAT/BCMath.h>

// Helper function to get a name for a FreeAmplitude
const std::string free_amplitude_name(const yap::FreeAmplitude& fa) noexcept {
    return yap::to_string(*fa.decayChannel())
        + " L = " + std::to_string(fa.spinAmplitude()->L())
        + " S = " + yap::spin_to_string(fa.spinAmplitude()->twoS());
}

// ---------------------------------------------------------
//Fit::Fit(const std::string& name,
//             std::unique_ptr<yap::Model> m,
//             const unsigned integration_points,
//             const unsigned integration_threads,
//             const unsigned integration_batch_size) :
//    BCModel(name),
//    Model_(std::move(m)),
//    FitData_(model()->createDataSet()),
//    FitPartitions_(1, &FitData_), // 1 partition
//    Axes_(model()->massAxes()), // Use the default axes.
//    NumberOfIntegrationPoints_(integration_points),
//    NumberOfIntegrationThreads_(integration_threads),
//    IntegrationBatchSize_(integration_batch_size),
//    Integral_(*model()),
//    IntegralData_(model()->createDataSet()),
//    IntegralDataPartitions_(1, &IntegralData_)

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

        const double upper_amp_range = 3. * abs(fa->value());
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
    FitModel_->model()->setParameterFlagsToUnchanged();

    return L;
}

void Fit::setParameters(const std::vector<double>& p) {
    assert(p.size() == 2 * FitModel_->freeAmplitudes().size());

    // Updates the values of the free amplitudes.
    double cumulative_phase = 0.;
    for (size_t i = 0; i < FitModel_->freeAmplitudes().size(); ++ i) {
        *FitModel_->freeAmplitudes()[i] = std::polar<double>(p[2 * i], yap::rad<double>(p[2 * i + 1] + cumulative_phase));
        cumulative_phase += p[2 * i + 1];
    }

    // Evaluate the integral of the model (with the new parameters).
    Integrator_->integrate();
}

//void Fit::integrate() {
//    assert(IntegrationPointGenerator_);
//    yap::ImportanceSampler::calculate(Integral_,
//                                      integrationPointGenerator(),
//                                      numberOfIntegrationPoints(),
//                                      integrationBatchSize(),
//                                      numberOfIntegrationThreads());
//}

// ---------------------------------------------------------
// double Fit::LogAPrioriProbability(const std::vector<double>& pars)
// {
//     // return the log of the prior probability p(pars)
//     // If you use built-in priors, leave this function commented out.
// }

// Helper function to find the index of a FreeAmplitude in the parameter vector.
// TODO make FreeAmplitude constant
const size_t free_amplitude_idx(const std::shared_ptr<yap::FreeAmplitude>& fa,
                                const yap::FreeAmplitudeVector& fas) {
    using namespace std;
    const auto it = find(begin(fas), end(fas), fa);

    if (it == end(fas))
        throw yap::exceptions::Exception("FreeAmplitude not found", "free_amplitude_idx");

    return 2 * distance(begin(fas), it);
}

const size_t free_amplitude_idx(const std::shared_ptr<yap::FreeAmplitude>& fa,
                                const std::shared_ptr<FitModel>& fit_model) {
    return free_amplitude_idx(fa, fit_model->freeAmplitudes());
}

// ---------------------------------------------------------
 void Fit::CalculateObservables(const std::vector<double>& p)
 {
     // Use the values of the free amplitudes that I've already set
     // in the setParameters() function, instead of evaluating them again.
     for (size_t i = 0; i < FitModel_->freeAmplitudes().size(); ++ i) {
         const auto A = FitModel_->freeAmplitudes()[i]->value();

         // TODO
         // for( fa : fas ) {
         //    i = distance(fas.begin(), fa);
         // }
         GetObservables()[i * 2]     = std::real(A);
         GetObservables()[i * 2 + 1] = std::imag(A);
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


void Fit::fixComponent(const std::shared_ptr<yap::FreeAmplitude>& fa, const unsigned char j, const double value) {
    if (j > 1)
        throw yap::exceptions::Exception("FreeAmplitude only has two components", "Fit::fixComponent");

    const auto i = free_amplitude_idx(fa, FitModel_);
    GetParameter(i + j).Fix(value);
}

void Fit::setRanges(const std::shared_ptr<yap::FreeAmplitude>& fa, double first_low, double first_high, double second_low, double second_high) {
    const auto i = free_amplitude_idx(fa, FitModel_);
    GetParameter(i + 0).SetLimits(first_low,  first_high);
    GetParameter(i + 1).SetLimits(second_low, second_high);
}
