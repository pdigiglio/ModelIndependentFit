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
// for to_string
#include <DecayTree.h>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iterator>
#include <memory>
#include <vector>

FitModel::FitModel(std::unique_ptr<yap::Model> m,
                   const std::string model_name) :
    Name_(model_name),
    Model_(std::move(m)),
    MassAxes_(model()->massAxes()),
    MassRanges_(yap::squared(mass_range(FitModel::Dmass(), massAxes(), model()->finalStateParticles()))),
    FreeAmplitudes_()
{
    // Check if the model is valid
    if (!model())
        throw yap::exceptions::Exception("Model nullptr", "FitModel::FitModel");

    // Check for consistency
    model()->lock();
    if (!model()->consistent())
        throw yap::exceptions::Exception("Model inconsistent", "FitModel::FitModel");

    // Check if the mass axes are not empty.
    if (massAxes().empty())
        throw yap::exceptions::Exception("Mass axes empty", "FitModel::FitModel");

    // Check if the mass ranges are not empty.
    if (massRanges().empty())
        throw yap::exceptions::Exception("Mass ranges empty", "FitModel::FitModel");

    // Copy the non-fixed FreeAmplitude's in the internal storage.
    const auto model_fas = yap::free_amplitudes(*model());
    std::copy_if(std::begin(model_fas), std::end(model_fas), std::back_inserter(FreeAmplitudes_),
                 [](const typename yap::FreeAmplitudeVector::value_type& fa)
                 { return fa->variableStatus() != yap::VariableStatus::fixed; } );
}

FitModel::~FitModel() = default;

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
    constexpr unsigned bins1 = 8;
    auto p1 = partition_range(low_range, 0.823231, bins1);

    constexpr unsigned bins2 = 14;
    auto p2 = partition_range(0.823231, 1.18596, bins2);

    constexpr unsigned bins3 = 8;
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


std::shared_ptr<const yap::FinalStateParticle> pi_plus(const std::unique_ptr<yap::Model>& m) {
#ifndef NDEBUG
    std::cout << " > Fetching pi+ from model" << std::endl;
#endif
    auto piPlus = yap::particle(*m, yap::is_named("pi+"));
    assert(piPlus != nullptr);
    return std::static_pointer_cast<const yap::FinalStateParticle>(piPlus);
}

const std::shared_ptr<const yap::FinalStateParticle> pi_minus(const std::unique_ptr<yap::Model>& m) {
#ifndef NDEBUG
    std::cout << " > Fetching pi- from model" << std::endl;
#endif
    auto piMinus = yap::particle(*m, yap::is_named("pi-"));
    assert(piMinus != nullptr);
    return std::static_pointer_cast<const yap::FinalStateParticle>(piMinus);
}

const std::shared_ptr<const yap::DecayingParticle> f_0(const std::unique_ptr<yap::Model>& m) {
#ifndef NDEBUG
    std::cout << " > Fetching f_0 from model" << std::endl;
#endif
    auto f_0 = yap::particle(*m, yap::is_named("f_0"));
    assert(f_0 != nullptr);
    return std::static_pointer_cast<const yap::DecayingParticle>(f_0);
}

// Helper function to get the partitioning of a mass range of a model.
// Note: make sure that the model has been locked before!
const std::vector<double> mass_partition_of_model(const std::unique_ptr<yap::Model>& data_model) {
    assert(data_model != nullptr);
    assert(data_model->locked());

    // Get the final state particles to determine the phase-space region.
    auto piPlus  = pi_plus(data_model);
    auto piMinus = pi_minus(data_model);

    return partition_mass_axis(std::abs(piPlus->mass() + piMinus->mass()),
                               std::abs(FitModel::Dmass() - piPlus->mass()));
}

std::unique_ptr<yap::Model> d3pi_binned() {
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

    // -----------------------------------------------------------------
    // Use the model that generated the MC data to set the initial
    // values for the free amplitudes of the bins.
    {
        const auto data_model = d3pi();

        const auto f0_BW = std::static_pointer_cast<const BreitWigner>(f_0(data_model)->massShape());
        const double f0_mass  = f0_BW->mass()->value();
        const double f0_width = f0_BW->width()->value();
        auto bw = [=](double s) {
            return 1. / std::complex<double>(f0_mass * f0_mass - s, - f0_mass * f0_width);
        };

        int i = 0;
        for (const auto& fa : free_amplitudes(*M, from(*D), is_not_fixed())) {
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

std::unique_ptr<yap::Model> d3pi()
{
    using namespace std;
    using namespace yap;

    // Create an empty model
    auto M = make_model<ZemachFormalism>();

	// Open the particle table.
    auto F = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");

    // Create the final state particles.
    auto piPlus  = FinalStateParticle::create(F[ 211]);
    auto piMinus = FinalStateParticle::create(F[-211]);

	// Set the final state of the decay.
    M->setFinalState(piPlus, piMinus, piPlus);

    // Create the initial state particle.
    auto D = DecayingParticle::create(F["D+"], FitModel::radialSize());

	// Vector of the intermediate resonance names.
	const std::vector<std::string> resonance_names({"f_0", "rho0"});
	// Create the decay tree based on the previous resonance names.
	for (const auto& rn : resonance_names) {
		// Create the resonance.
		auto r = DecayingParticle::create(F[rn], FitModel::radialSize(), make_shared<BreitWigner>(F[rn]));

		// Set the resonance daughter particles.
		r->addWeakDecay(piPlus, piMinus);

		// Broaden the f_0's width.
		if (rn == "f_0") {
			auto w = static_pointer_cast<BreitWigner>(r->massShape())->width();
			w->setValue(4. * w->value());
//			std::cout << "Width " << w->value() << std::endl;
		}

		// Add the resonance to the D decay.
		D->addWeakDecay(r, piPlus);
		std::cout << rn << " resonance added to the D decay." << std::endl;
	}

	// Set the decay initial state.
    M->addInitialState(D);
    // Lock the model before returning it.
    M->lock();

    return M;
}

const std::vector<double> guess_parameters(Fit& m) {
    // The model according to which the MC data were generated.
    const auto data_model = d3pi();

    // Get the parameters of the f0 resonance to reproduce its shape.
    const auto f0_BW    = std::static_pointer_cast<const yap::BreitWigner>(f_0(data_model)->massShape());
    const auto BW_width = f0_BW->width()->value();
    const auto BW_mass  = f0_BW->mass()->value();

    // Reproduce the mass shape of the f0 resonance.
    auto bw = [=](const double s) {
        return 1. / std::complex<double>(BW_mass * BW_mass - s, -BW_mass * BW_width);
    };

    // Get the mass-range partition.
    const auto p = mass_partition_of_model(m.model());

    auto pi = []() { return std::acos(-1); };

//    // Find the fixed parameters.
//    std::vector<int> fixed_parameter_indices;
//    for (const auto& parameter : m.GetParameters())
//        if (parameter.Fixed())
//            fixed_parameter_indices.push_back(std::distance(std::begin(m.GetParameters()), parameter));
//
//    assert(fixed_parameter_indices.size() == 2);
//    for (const auto& fpi : fixed_parameter_indices)
//        assert(!(fpi < 0));
//
//    assert(false);


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


std::unique_ptr<FitModel> make_fit_model(const std::string& model_name) {
    return std::make_unique<FitModel>(d3pi_binned(), model_name);
}

// Helper function to find the index of a FreeAmplitude in the parameter vector.
const size_t free_amplitude_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                  const yap::FreeAmplitudeVector& fas) {
    using namespace std;
    const auto it = find(begin(fas), end(fas), fa);

    if (it == end(fas))
        throw yap::exceptions::Exception("FreeAmplitude not found", "free_amplitude_idx");

    return distance(begin(fas), it);
}

const size_t free_amplitude_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                  const std::shared_ptr<const FitModel>& fit_model) {
    return free_amplitude_index(fa, fit_model->freeAmplitudes());
}
