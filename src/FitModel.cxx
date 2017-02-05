/**
 *
 *    @file  FitModel.cxx
 *   @brief  
 *
 *    @date  02/03/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "FitModel.h"

#include "MassBin.h"
#include "Fit.h"

#include <Attributes.h>
#include <Exceptions.h>
#include <BreitWigner.h>
#include <FreeAmplitude.h>
#include <MassAxes.h>
#include <MassRange.h>
#include <Model.h>
#include <Parameter.h>
#include <ZemachFormalism.h>
#include <PDL.h>
#include <FinalStateParticle.h>
#include <DecayingParticle.h>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iterator>
#include <memory>
#include <vector>

FitModel::FitModel(std::unique_ptr<yap::Model> m) :
    Model_(std::move(m)),
    MassAxes_(model()->massAxes()),
    MassRanges_(yap::squared(mass_range(FitModel::Dmass(), axes(), model()->finalStateParticles()))),
    FreeAmplitudes_()
{
    // Check if the model is valid
    if (!model())
        throw yap::exceptions::Exception("Model nullptr", "FitModel::FitModel");

    // Check for consistency
    model()->lock();
    if (!model()->consistent())
        throw yap::exceptions::Exception("Model inconsistent", "FitModel::FitModel");

    // Copy the non-fixed FreeAmplitude's in the internal storage.
    const auto model_fas = yap::free_amplitudes(*model());
    std::copy_if(std::begin(model_fas), std::end(model_fas), std::back_inserter(FreeAmplitudes_),
                 [](const typename yap::FreeAmplitudeVector::value_type& fa)
                 { return fa->variableStatus() != yap::VariableStatus::fixed; } );
}

FitModel::~FitModel() = default;

// Helper function.
// Splits a range _[min, max)_ into a partition with _bins_ number of bins.
std::vector<double> partition_range(double min, const double max, const unsigned int bins)
{
    std::vector<double> partition;
    partition.reserve(bins + 1);

    const double step = (max - min) / bins;
//    std::cout << "Bin width = " << step << std::endl;
    for (unsigned int i = 0; i < bins; ++i) {
        partition.push_back(min);
        min += step;
    }
    // insert last element
    partition.push_back(max);

    return partition;
}

const std::vector<double> partition_mass_axis(double low_range, double up_range) {
    assert(low_range < up_range);

    // partition
    constexpr unsigned bins1 = 6;
    auto p1 = partition_range(low_range, 0.823231, bins1);

    constexpr unsigned bins2 = 10;
    auto p2 = partition_range(0.823231, 1.18596, bins2);

    constexpr unsigned bins3 = 6;
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

    for (int i = 0; i < static_cast<int>(p.size()) - 1; ++ i)
        std::cout << "bin(" << std::setw(2) << std::right << i << ") = ["
                  << std::setw(8) << std::left << p[i] << " , " << p[i+1] << ")" << std::endl;

//    // For debugging purpose
//    for (const auto& bin : p)
//        std::cout << bin << std::endl;

    return p;
}

std::unique_ptr<yap::Model> d3pi_binned() {
    using namespace std;
    using namespace yap;

    // Make a new model.
    auto M = make_model<ZemachFormalism>();

    auto F = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");

    // Final state particles.
    auto piPlus  = FinalStateParticle::create(F[ 211]);
    auto piMinus = FinalStateParticle::create(F[-211]);

    M->setFinalState(piPlus, piMinus, piPlus);

    // Initial state particle.
    auto D = DecayingParticle::create(F["D+"], FitModel::radialSize());

    // Partition the mass range.
    auto p = partition_mass_axis(std::abs(piPlus->mass() + piMinus->mass()),
                                 std::abs(F["D+"].mass() - piPlus->mass()));

    // The number of bins in the partition.
    const auto bins = p.size() - 1;

    // vector of resonance pointers
    for (unsigned i = 0; i < bins; ++i) {
        const auto resonance_name = "Bin(" + to_string(i) + ")";
        auto r = DecayingParticle::create(resonance_name, QuantumNumbers(0, 0), FitModel::radialSize(), make_shared<MassBin>(p[i], p[i+1]));
        r->addWeakDecay(piMinus, piPlus);
        D->addWeakDecay(r, piPlus);
    }

    {
        // Add the rho0 resonance.
        auto r = DecayingParticle::create(F["rho0"], FitModel::radialSize(), std::make_shared<BreitWigner>(F["rho0"]));
        r->addWeakDecay(piMinus, piPlus);
        D->addWeakDecay(r, piPlus);

        // Set this free amplitude as fixed.
        free_amplitude(*D, to(r))->variableStatus() = VariableStatus::fixed;
        // Fix the amplitude to its default value, which is the complex (1, 0).
    }

    M->addInitialState(D);

    // Check that each bin has an associated free amplitude
    assert(bins == free_amplitudes(*M, from(D), is_not_fixed()).size());

    constexpr double bw_M  = 1.;
    constexpr double Gamma = .2;
    auto bw = [=](double s) { return 1. / std::complex<double>(bw_M * bw_M - s, - bw_M * Gamma ); };
    int i = 0;
    for (const auto& fa : free_amplitudes(*M, from(*D), is_not_fixed())) {

        // evaluate BW on the low edge of the bin
        const auto bw_value = bw(p[i] * p[i]);
//        std::cout << "[" << p[i] << p[i + 1] << ") BW("
//                  << p[i] * p[i] << ") = " << bw_value << std::endl;

        *fa = polar(abs(bw_value), arg(bw_value));
//        std::cout << "i = " << i << std::endl;
//        std::cout << to_string(*fa) << std::endl;

        ++i;
    }

    return M;
}

const std::vector<double> guess_parameters(const Fit& m) {
    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");
    auto r = yap::DecayingParticle::create(F["f_0"], 3., std::make_shared<yap::BreitWigner>(F["f_0"]));
    auto w = std::static_pointer_cast<yap::BreitWigner>(r->massShape())->width();
    w->setValue(4. * w->value());
    // Breit-Wigner f_0 width
    const auto BW_width = w->value();
    // f_0 mass
    const auto BW_mass  = std::static_pointer_cast<yap::BreitWigner>(r->massShape())->mass()->value();
    auto bw = [=](const double s) { return 1. / std::complex<double>(BW_mass * BW_mass - s, -BW_mass * BW_width); };
    auto pi = []() { return std::acos(-1); };

    
    auto piPlus  = yap::FinalStateParticle::create(F[ 211]);
    auto piMinus = yap::FinalStateParticle::create(F[-211]);
    auto p = partition_mass_axis(std::abs(piPlus->mass() + piMinus->mass()),
                                 std::abs(F["D+"].mass() - piPlus->mass()));

    std::cout << " guess param partition : " << p.size() << std::endl;
    //
    // The fixed value of the second bin's amplitude.
    const auto fixed_amplitude_idx = m.GetParameters().Size() - 2;
    const auto normalization = 1. / std::abs(bw(p[fixed_amplitude_idx/2] * p[fixed_amplitude_idx/2]));

    assert(2 * (p.size() - 1) == m.GetParameters().Size());

    std::cout << "Parameters:      " << m.GetParameters().Size() << std::endl;
    // run mode finding; by default using Minuit
    std::vector<double> guess(m.GetParameters().Size(), 0);
    unsigned fixed_parameters = 0;

    // The initial value is such that the first phase difference is 0.
    const double phase_shift = std::arg(bw(p[0] * p[0])) * 180 / pi();
    double cumulative_phase  = phase_shift;
    std::ofstream par_guess("output/par_guess.txt", std::ios::out); 
    par_guess << "#low_mass  abs        shifted phase (deg) phase (deg) delta_phase" << std::endl;
    for (size_t i = 0; i < p.size() - 1; ++ i) {

        // Evaluate the value of the BW on the lower edge of the bin.
        const auto BW_value = bw(p[i] * p[i]);
        const auto BW_magnitude = std::abs(BW_value);
        const auto BW_ph_diff   = std::arg(BW_value) * 180 / pi() - cumulative_phase;

//        std::cout << std::setw(9) << p[i] << " "
//                  << std::setw(9) << std::real(BW_value) << " "
//                  << std::setw(9) << std::imag(BW_value) << " "
//                  << std::setw(9) << std::abs(BW_value)  << " "
//                  << std::setw(9) << std::arg(BW_value) * 180 / pi() << std::endl; 
        // Store magnitude and phase difference in a vector (useful for indeces).
        const std::vector<double> guess_pars({normalization * BW_magnitude, BW_ph_diff});

        par_guess << std::left
                  << std::setw(10) << p[i] << " "
                  << std::setw(10) << guess_pars[0] << " "
                  << std::setw(19) << guess_pars[1] + cumulative_phase - phase_shift << " "
                  << std::setw(11) << guess_pars[1] + cumulative_phase << " "
                  << std::setw(10) << guess_pars[1] << std::endl;

        for (size_t j = 0; j < 2; ++ j) {
            const auto idx = 2 * i + j;

            const auto& parameter = m.GetParameter(idx);
            if (idx == 1 || idx == fixed_amplitude_idx) {
                assert(parameter.Fixed());
                guess[idx] = parameter.GetFixedValue();
                assert(guess[idx] == guess_pars[j]);
                ++ fixed_parameters;
            } else {
                assert(!parameter.Fixed());
                guess[idx] = guess_pars[j];
            }
        }

        cumulative_phase += BW_ph_diff;
    }

    assert(fixed_parameters == 2);

//    std::cout << "Free Parameters: " << m.GetParameters().Size() - fixed_parameters << std::endl;
//
//    throw;

    return guess;
}

void write_fit_result_to_file(Fit& m) {
    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");
    auto piPlus  = yap::FinalStateParticle::create(F[ 211]);
    auto piMinus = yap::FinalStateParticle::create(F[-211]);
    auto p = partition_mass_axis(std::abs(piPlus->mass() + piMinus->mass()),
                                 std::abs(F["D+"].mass() - piPlus->mass()));

    std::ofstream par_fit("output/par_fit.txt", std::ios::out);
    par_fit << "#low_mass    amp          d_amp        phase_diff   d_phase_diff phase        d_phase" << std::endl;
    const auto bins = p.size() - 1;
    std::cout << "bins : " << bins << " parameters : " << m.GetBestFitParameters().size() << std::endl;

    // Chech that there are 2 real parameters for each bin.
    // This fails because the parameters also include the observables.
//    assert(bins == 2 * m.GetBestFitParameters().size());

    double cumulative_phase          = 0.;
    double cumulative_phase_variance = 0.;
    for (size_t i = 0; i < bins; ++ i) {
        const auto phi_i       = m.GetBestFitParameters()[2 * i + 1];
        const auto delta_phi_i = m.GetBestFitParameterErrors()[2 * i + 1];
        
        cumulative_phase          += phi_i;
        cumulative_phase_variance += delta_phi_i * delta_phi_i;

        par_fit << std::left
                << std::setw(12) << p[i] << " "
                << std::setw(12) << m.GetBestFitParameters()[2 * i] << " "
                << std::setw(12) << m.GetBestFitParameterErrors()[2 * i] << " "
                << std::setw(12) << m.GetBestFitParameters()[2 * i + 1] << " "
                << std::setw(12) << m.GetBestFitParameterErrors()[2 * i + 1] << " "
                << std::setw(12) << cumulative_phase << " "
                << std::setw(12) << std::sqrt(cumulative_phase_variance) << std::endl;
    }
}


std::unique_ptr<FitModel> make_fit_model() {
    return std::make_unique<FitModel>(d3pi_binned());
}
