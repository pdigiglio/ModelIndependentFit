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
#include "MassRangePartition.h"
#include "ModelIndependentFitModel.h"
#include "RealTimePlotModelIndependentFitModel.h"
#include "RealTimeParameterPlot.h"

#include <DecayChannel.h>

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
#include <utility>
#include <vector>

// Helper function to create a partition with bins that have different widths.
std::unique_ptr<const MassRangePartition> partition_mass_axis(double low_range, double up_range) {
    assert(low_range < up_range);

    auto mrp = std::make_unique<const MassRangePartition>(std::vector<MassRangePartition::region>(
                { MassRangePartition::region(low_range, 2),
                  MassRangePartition::region( .8,  14),
                  MassRangePartition::region(1.05,  4),
                  MassRangePartition::region(1.4,  10),
                  MassRangePartition::region(1.6,   2),
                  MassRangePartition::region(up_range, 0) }));

#ifndef NDEBUG
    for (size_t i = 0; i < mrp->numberOfBins(); ++ i) {
        std::cout << "bin(" << std::setw(2) << std::right << i << ") = ["
                            << std::setw(8) << std::left << mrp->massPartition()[i]
                  << " , "  << std::setw(8) << std::left << mrp->massPartition()[i+1]
                  << ")" << std::endl;
    }
#endif

    return mrp;
}

// Helper function to create a binned fit model.
std::shared_ptr<ModelIndependentFitModel> create_d3pi_binned(const std::string& model_name) {
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

    auto p = partition_mass_axis(std::abs(piPlus->mass() + piMinus->mass()),
                                 std::abs(F["D+"].mass() - piPlus->mass()));

    // Add the mass bins to the fit model.
    for (unsigned i = 0; i < p->numberOfBins(); ++i) {
        const auto resonance_name = "Bin(" + to_string(i) + ")";
        auto r = DecayingParticle::create(resonance_name, QuantumNumbers(0, 0), FitModel::radialSize(), make_shared<MassBin>(p->massPartition()[i], p->massPartition()[i+1]));
        assert(r != nullptr);

        r->addWeakDecay(piMinus, piPlus);
        D->addWeakDecay(r, piPlus);
    }

    {
//        // Add the rho0 resonance.
//        auto rho0 = F["rho0"];
//        auto r = DecayingParticle::create(rho0, FitModel::radialSize(), std::make_shared<BreitWigner>(rho0));
//
//        r->addWeakDecay(piMinus, piPlus);
//        D->addWeakDecay(r, piPlus);

        // ---------------------------------
        // Add the initials state
        M->addInitialState(D);
        // Lock the model.
        M->lock();
        // ---------------------------------

//        // Set this free amplitude as fixed.
//        free_amplitude(*D, to(r))->variableStatus() = VariableStatus::fixed;
//        // Fix the amplitude to its default value, which is the complex (1, 0).
    }

    // Check that each bin has an associated free amplitude
    assert(p->numberOfBins() == free_amplitudes(*M, from(D), is_not_fixed()).size());

    const std::vector<double> guess(2 * p->numberOfBins(), 0);
    return std::make_unique<ModelIndependentFitModel>(std::move(M), std::move(p),  model_name);
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

const std::shared_ptr<const yap::DecayingParticle> f_0_1500(const std::unique_ptr<yap::Model>& m) {
#ifndef NDEBUG
    std::cout << " > Fetching f0(1500) from model" << std::endl;
#endif
    auto f_0 = yap::particle(*m, yap::is_named("f_0(1500)"));
    assert(f_0 != nullptr);
    return std::static_pointer_cast<const yap::DecayingParticle>(f_0);
}

// Helper function to get the f0(1500)'s resonance mass shape.
std::function<std::complex<double>(double)> f0_1500_mass_shape() {
    const auto data_model = d3pi();
    assert(data_model->locked());

    const auto f0_BW = std::static_pointer_cast<const yap::BreitWigner>(f_0_1500(data_model)->massShape());
    const double f0_mass  = f0_BW->mass()->value();
    const double f0_width = f0_BW->width()->value();
    return [=](double s) { return 1. / std::complex<double>(f0_mass * f0_mass - s, - f0_mass * f0_width); };
}

// Helper function to return the mass shape of the S wave.
std::function<std::complex<double>(double)> S_wave_mass_shape() {
    const auto f0_ms      = f0_mass_shape();
    const auto f0_1500_ms = f0_1500_mass_shape();

    return [=](double s) { return f0_ms(s) + f0_1500_ms(s); };
}

std::shared_ptr<RealTimePlotModelIndependentFitModel> binned_d3pi(const std::string model_name) {
    auto M = create_d3pi_binned(model_name);
    assert(M->model()->locked());
    assert(M->massBinSorted());

    const auto step_mass_shape = binned_mass_shape(S_wave_mass_shape(), *M);
    M->setParameters(step_mass_shape);
    auto rtpp = std::make_unique<RealTimeParameterPlot>(bins(0., *M->massRangePartition()), M->freeAmplitudes());

    return std::make_shared<RealTimePlotModelIndependentFitModel>(std::move(*M), std::move(rtpp));
}

std::shared_ptr<ModelIndependentFitModel> binned_d3pi_from_file(const std::string& file_name,
                                                                const std::string model_name) {
    // Make a binned model
    auto M = create_d3pi_binned(model_name);
    assert(M->model()->locked());
    assert(M->massBinSorted());

    // Get the non-fixed free amplitudes.
    const auto fas = M->freeAmplitudes();

    // Vector of fitted free amplitudes.
    std::vector<FittedFreeAmplitude> ffa_vector;
    ffa_vector.reserve(fas.size());

    {
        // Check that the file can be opened.
        std::ifstream par_fit(file_name, std::ios::in);
        if (!par_fit)
            throw yap::exceptions::Exception("Can't open " + file_name, "binned_d3pi_from_file");

        // Load the entries in the vector.
        std::copy(FitResultFileIterator(par_fit), FitResultFileIterator::end(),
                  std::back_inserter(ffa_vector));
    }

    // Check that the number of lines and the number of free amplitudes match.
    const auto lines = ffa_vector.size();
    if (fas.size() != lines) {
        const std::string error_string = "The number of FreeAmplitudes (" + std::to_string(fas.size()) +
                                         ") doesn't correspond to the file lines (" + std::to_string(lines) + ")";

        throw yap::exceptions::Exception(error_string, "binned_d3pi_from_file");
    }

#ifndef NDEBUG
    {
        // Check if the bin low edge of the mass bin and the one of the fitted free
        // amplitude match within a certain tollerance.
        const auto ble_matches = [] (const auto& mass_bin, const auto& ffa) {
            return std::abs(mass_bin->lowerEdge()->value() - ffa.bin_low_edge) < 1e-5;
        };

        // Check if the bin low edges of the file match the ones in the mass bins.
        assert(std::all_of(std::begin(fas), std::end(fas), [&](const auto& fa)
               { return ble_matches(bin_mass_shape(fa, *M), ffa_vector[free_amplitude_index(fa, *M)]); }));
    }
#endif

    // Get the complex-parameter vector from the input file.
    std::vector<std::complex<double>> p;
    p.reserve(fas.size());

    std::transform(std::begin(ffa_vector), std::end(ffa_vector), std::back_inserter(p),
                   [](const auto& ffa) { return std::polar(ffa.amplitude, yap::rad(ffa.cumulative_phase)); });

    M->fixParameters(p);

    return M;
}

std::vector<double> guess_parameters(Fit& m) {
    const auto& fas = m.fitModel()->freeAmplitudes();
    // DANGEROUS: this may not be castable to a ModelIndependentFitModel!!
    const auto fit_model = std::static_pointer_cast<const ModelIndependentFitModel>(m.fitModel());

#ifndef NDEBUG
    // Check if the initial values of the free amplitudes match the guessed initial parameters.
    {
        // Make the mass shape a step function.
        const auto step_mass_shape = binned_mass_shape(S_wave_mass_shape(), *fit_model);
        assert(std::equal(std::begin(fas), std::end(fas), std::begin(step_mass_shape),
                          [](const auto& fa, const auto& sms) { return fa->value() == sms; }));
    }
#endif

    // -- Open output file ------------------------------------------
    std::ofstream par_guess("output/par_guess.txt", std::ios::out);
    par_guess << FittedFreeAmplitude::header() << std::endl;
    // --------------------------------------------------------------

    const auto initial_fit_guess = yap_to_fit_parameters(fas);
    assert(initial_fit_guess.size() == 2 * fas.size());

    // Take into account a possible phase shifting.
    const double phase_shift = yap::deg(std::arg(fas[0]->value())) - initial_fit_guess[1];

    // Cumulative phase: sum of the phase differences.
    double cp = 0.;
    for (size_t i = 0; i < fas.size(); ++ i) {
        assert(std::abs(fas.at(i)->value()) == initial_fit_guess.at(2 * i));

        // Bin low edge.
        const auto ble = fit_model->massPartition().at(i);
        // Free-amplitude amplitude.
        const auto a   = initial_fit_guess.at(2 * i);
        // Free-amplitude phase difference.
        const auto pd  = initial_fit_guess.at(2 * i + 1);
        cp += pd;

        // Make sure that the phases match
        assert(std::abs((cp + phase_shift) - yap::deg(std::arg(fas.at(i)->value()))) < 1e-7);

        par_guess << FittedFreeAmplitude(ble, a, 0, pd, 0, cp, 0) << std::endl;
    }

#ifndef NDEBUG
    constexpr std::array<size_t, 1> fixed_parameters = { 1 };
    for (size_t i = 0; i < initial_fit_guess.size(); ++ i) {
        if (std::find(std::begin(fixed_parameters), std::end(fixed_parameters), i) != std::end(fixed_parameters)) {
            // Make sure it's fixed.
            assert(m.GetParameter(i).Fixed());
            std::cout << m.GetParameter(i).GetFixedValue() << " " <<  initial_fit_guess[i] << std::endl;
            assert(m.GetParameter(i).GetFixedValue() == initial_fit_guess[i]);
        } else {
            assert(!m.GetParameter(i).Fixed());
        }
    }
#endif

    return initial_fit_guess;
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

        // Construct a temporary FittedFreeAmplitude to easily stream it.
        par_fit << FittedFreeAmplitude(p[i], amp_i, delta_amp_i, phi_i, delta_phi_i, cumulative_phase, std::sqrt(cumulative_phase_variance)) << std::endl;
    }

    // XXX This would break the generation of the model from the fit results.
//    par_fit << FittedFreeAmplitude::footer() << std::endl;
    // -------------------------------------------------------------------
}
