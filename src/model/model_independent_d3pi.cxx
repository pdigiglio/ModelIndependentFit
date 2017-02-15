/**
 *
 *    @file  model_independent_d3pi.cxx
 *   @brief  
 *
 *    @date  02/14/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "model/model_independent_d3pi.h"

#include "model/d3pi.h"

#include "Fit.h"
#include "FitResultFileIterator.h"
#include "FittedFreeAmplitude.h"
#include "MassBin.h"
#include "ModelIndependentFitModel.h"

#include <Attributes.h>
#include <BreitWigner.h>
#include <DecayTree.h> // for to_string
#include <DecayingParticle.h>
#include <Exceptions.h>
#include <FinalStateParticle.h>
#include <FreeAmplitude.h>
#include <MassAxes.h>
#include <MassRange.h>
#include <Model.h>
#include <PDL.h>
#include <Parameter.h>
#include <ZemachFormalism.h>

#include <memory>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>

// Helper function.
// Splits a range _[min, max)_ into a partition with _bins_ number of bins.
std::vector<double> partition_range(double min, const double max, const unsigned int bins) {
    assert(bins > 1);
    std::vector<double> partition;
    partition.reserve(bins + 1);

    assert(min < max);
    const double step = (max - min) / bins;
    for (unsigned int i = 0; i < bins; ++i) {
        partition.push_back(min);
        min += step;
    }

    // insert last element
    partition.push_back(max);

    return partition;
}

// Helper function to create a partition with bins that have different widths.
const std::vector<double> partition_mass_axis(double low_range, double up_range) {
    assert(low_range < up_range);

    // partition
    constexpr unsigned bins1 = 3;
    auto p1 = partition_range(low_range, 0.823231, bins1);

    constexpr unsigned bins2 = 4;
    auto p2 = partition_range(0.823231, 1.18596, bins2);

    constexpr unsigned bins3 = 3;
    auto p3 = partition_range(1.18596, up_range, bins3);

    // number of bins
    constexpr unsigned bins = bins1 + bins2 + bins3;
    std::vector<double> p;
    p.reserve(bins);
    p.insert(p.end(), std::make_move_iterator(p1.begin()), std::make_move_iterator(p1.end()));
    p.insert(p.end(), std::make_move_iterator(p2.begin() + 1), std::make_move_iterator(p2.end()));
    p.insert(p.end(), std::make_move_iterator(p3.begin() + 1), std::make_move_iterator(p3.end()));

    // It should be sorted, but just to make sure...
    std::sort(p.begin(), p.end());

#ifndef NDEBUG
    for (int i = 0; i < static_cast<int>(p.size()) - 1; ++ i)
        std::cout << "bin(" << std::setw(2) << std::right << i << ") = ["
            << std::setw(8) << std::left << p[i] << " , " << p[i+1] << ")" << std::endl;
#endif

    return p;
}

// Helper function to create a binned fit model.
std::unique_ptr<ModelIndependentFitModel> create_d3pi_binned(const std::string& model_name) {
    using namespace std;
    using namespace yap;

    // Make a new model.
    auto M = make_model<ZemachFormalism>();

    auto F = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");

    // Initial state particle.
    auto D = DecayingParticle::create(F["D+"], FitModel::radialSize());

    // Final state particles.
    auto piPlus  = FinalStateParticle::create(F[ 211]);
    auto piMinus = FinalStateParticle::create(F[-211]);
    M->setFinalState(piPlus, piMinus, piPlus);

    const auto p = partition_mass_axis(std::abs(piPlus->mass() + piMinus->mass()),
                                       std::abs(F["D+"].mass() - piPlus->mass()));

    // The number of bins in the partition.
    const auto bins = p.size() - 1;

    // Add the mass bins to the fit model.
    for (unsigned i = 0; i < bins; ++i) {
        const auto resonance_name = "Bin(" + to_string(i) + ")";
        auto r = DecayingParticle::create(resonance_name, QuantumNumbers(0, 0), FitModel::radialSize(), make_shared<MassBin>(p[i], p[i+1]));
        assert(r != nullptr);

        r->addWeakDecay(piMinus, piPlus);
        D->addWeakDecay(r, piPlus);
    }

    {
        // Add the rho0 resonance.
        auto rho0 = F["rho0"];
        auto r = DecayingParticle::create(rho0, FitModel::radialSize(), std::make_shared<BreitWigner>(rho0));
        assert(r != nullptr);

        r->addWeakDecay(piMinus, piPlus);
        D->addWeakDecay(r, piPlus);

        // ---------------------------------
        // Add the initials state
        M->addInitialState(D);
        // Lock the model.
        M->lock();
        // ---------------------------------

        // Set this free amplitude as fixed.
        free_amplitude(*D, to(r))->variableStatus() = VariableStatus::fixed;
        // Fix the amplitude to its default value, which is the complex (1, 0).
    }

    // Check that each bin has an associated free amplitude
    assert(bins == free_amplitudes(*M, from(D), is_not_fixed()).size());

    return std::make_unique<ModelIndependentFitModel>(std::move(M), p, model_name);
}

const std::shared_ptr<const yap::DecayingParticle> f_0(const std::unique_ptr<yap::Model>& m) {
#ifndef NDEBUG
    std::cout << " > Fetching f_0 from model" << std::endl;
#endif
    auto f_0 = yap::particle(*m, yap::is_named("f_0"));
    assert(f_0 != nullptr);
    return std::static_pointer_cast<const yap::DecayingParticle>(f_0);
}

// Helper function to get the f0's resonance mass shape.
std::function<std::complex<double>(double)> f0_mass_shape() {
    const auto data_model = d3pi();
    assert(data_model->locked());

    const auto f0_BW = std::static_pointer_cast<const yap::BreitWigner>(f_0(data_model)->massShape());
    const double f0_mass  = f0_BW->mass()->value();
    const double f0_width = f0_BW->width()->value();
    return [=](double s) { return 1. / std::complex<double>(f0_mass * f0_mass - s, - f0_mass * f0_width); };
}

std::unique_ptr<ModelIndependentFitModel> binned_d3pi(const std::string model_name) {
    using namespace std;
    using namespace yap;

    // Create a model.
    auto M = create_d3pi_binned(model_name);
    const auto p = M->massPartition();

    // -----------------------------------------------------------------
    // Use the model that generated the MC data to set the initial
    // values for the free amplitudes of the bins.
    {
        // Get the mass shape of the f0 according to the MC-data model.
        auto bw = f0_mass_shape();

        // Get the D+ particle from the model.
        const auto D = std::static_pointer_cast<DecayingParticle>(particle(*M->model(), is_named("D+")));

        int i = 0;
        for (const auto& fa : free_amplitudes(*M->model(), from(D), is_not_fixed())) {
            // Evaluate BW on the low edge of the bin.
            const auto bin_low_edge = p[i];
            const auto bw_value = bw(bin_low_edge * bin_low_edge);
#ifndef NDEBUG
            std::cout << "[" << p[i] << p[i + 1] << ") BW("
                << p[i] * p[i] << ") = " << bw_value << std::endl;
#endif
            *fa = polar(abs(bw_value), arg(bw_value));
            ++i;
        }
    }
    // -----------------------------------------------------------------

    return M;
}

std::unique_ptr<ModelIndependentFitModel> binned_d3pi_from_file(const std::string& file_name,
                                                                const std::string model_name) {
    // Make a binned model
    auto M = create_d3pi_binned(model_name);

    const auto D = std::static_pointer_cast<yap::DecayingParticle>(yap::particle(*M->model(), yap::is_named("D+")));
    const auto fas = M->freeAmplitudes(); //yap::free_amplitudes(*M, yap::from(D), yap::is_not_fixed());

    // Check that the number of lines and the number of free amplitudes match.
    {
        std::ifstream par_fit(file_name, std::ios::in);
        if (!par_fit)
            throw yap::exceptions::Exception("Can't open " + file_name, "d3pi_binned");

        const auto lines = count_lines(par_fit);
        if (fas.size() != lines)
            throw yap::exceptions::Exception("The number of FreeAmplitudes (" + std::to_string(fas.size()) +
                    ") doesn't correspond to the file lines (" + std::to_string(lines) + ")",
                    "d3pi_binned");
    }

    // No need to check for exceptions as I already did it.
    std::ifstream par_fit(file_name, std::ios::in);
    FitResultFileIterator it(par_fit);
    for (const auto& fa : fas) {
        // Assign the input value.
        assert(it != FitResultFileIterator::end());
        // Temporary copy of the fitted free amplitude.
        const auto ffa = *it;

#ifndef NDEBUG
        std::cout << ffa << std::endl;
#endif
        // Fix the amplitude to the fitted value.
        *fa = std::polar<double>(ffa.amplitude, ffa.cumulative_phase);
        fa->variableStatus() = yap::VariableStatus::fixed;

        // Read the next line.
        ++ it;
    }

    throw;

    return M;
}

const std::vector<double> guess_parameters(Fit& m) {
    // DANGEROUS! It's not said that fitModel is a ModelIndependentFitModel!
    const auto p = std::static_pointer_cast<const ModelIndependentFitModel>(m.fitModel())->massPartition();

    const auto bins = p.size() - 1;
    assert(2 * bins == m.GetParameters().Size());

    // Get the mass shape of the f0 according to the MC-data model.
    auto bw = f0_mass_shape();

    // The fixed value of the second bin's amplitude.
    const auto fixed_amplitude_idx = m.GetParameters().Size() - 2;
    const auto normalization = 1.; // / std::abs(bw(p[fixed_amplitude_idx/2] * p[fixed_amplitude_idx/2]));

    std::vector<double> guess(m.GetParameters().Size(), 0);
    unsigned fixed_parameters = 0;

    // The initial value is such that the first phase difference is 0.
    const double phase_shift = yap::deg<double>(std::arg(bw(p[0] * p[0])));
    double cumulative_phase  = 0.; //phase_shift;
    std::ofstream par_guess("output/par_guess.txt", std::ios::out); 
    par_guess << FittedFreeAmplitude::header() << std::endl;
    for (size_t i = 0; i < bins; ++ i) {

        // Evaluate the value of the BW on the lower edge of the bin.
        const auto BW_value = bw(p[i] * p[i]);
        const auto BW_magnitude = std::abs(BW_value);
        const auto BW_ph_diff   = yap::deg<double>(std::arg(BW_value)) - cumulative_phase;

        // Store magnitude and phase difference in a vector (useful for indeces).
        const std::vector<double> guess_pars({normalization * BW_magnitude, BW_ph_diff});

        // Wrap the values into a FittedFreeAmplitude to easily stream them to the file.
        par_guess << FittedFreeAmplitude(p[i], guess_pars[0], 0., guess_pars[1], 0., guess_pars[1] + cumulative_phase, 0.) << std::endl;

        for (size_t j = 0; j < 2; ++ j) {
            const auto idx = 2 * i + j;

            const auto& parameter = m.GetParameter(idx);
            //            if (idx == 1 ) { // || idx == fixed_amplitude_idx) {
            //                assert(parameter.Fixed());
            //                guess[idx] = parameter.GetFixedValue();
            //                assert(guess[idx] == guess_pars[j]);
            //                ++ fixed_parameters;
            //            } else {
            //                assert(!parameter.Fixed());
            guess[idx] = guess_pars[j];
            //            }
        }

        cumulative_phase += BW_ph_diff;
    }

//    assert(fixed_parameters == 2);
    return guess;
}

void write_fit_result_to_file(Fit& m) {
    // DANGEROUS! It's not said that fitModel is a ModelIndependentFitModel!
    const auto p = std::static_pointer_cast<const ModelIndependentFitModel>(m.fitModel())->massPartition();

    // Chech that there are 4 real parameters for each bin.
    // They're 2 for the observables and 2 for the fit parameters.
    const auto bins = p.size() - 1;
    assert(4 * bins == m.GetBestFitParameters().size());

    // -- Write file -----------------------------------------------------
    std::ofstream par_fit("output/par_fit.txt", std::ios::out);
    par_fit << FittedFreeAmplitude::header() << std::endl;

    double cumulative_phase          = 0.;
    double cumulative_phase_variance = 0.;
    for (size_t i = 0; i < bins; ++ i) {
        const auto amp_i       = m.GetBestFitParameters()[2 * i];
        const auto delta_amp_i = m.GetBestFitParameterErrors()[2 * i];
        const auto phi_i       = m.GetBestFitParameters()[2 * i + 1];
        const auto delta_phi_i = m.GetBestFitParameterErrors()[2 * i + 1];

        cumulative_phase          += phi_i;
        cumulative_phase_variance += delta_phi_i * delta_phi_i;

        // Construct a temporary FittedFreeAmplidude to easily stream it.
        par_fit << FittedFreeAmplitude(p[i], amp_i, delta_amp_i, phi_i, delta_phi_i, cumulative_phase, std::sqrt(cumulative_phase_variance)) << std::endl;
    }
    // -------------------------------------------------------------------
}
