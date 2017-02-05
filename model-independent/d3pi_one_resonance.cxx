/**
 *
 *    @file  d3pi_one_resonance.cxx
 *   @brief  
 *
 *    @date  09/30/16
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "d3pi_one_resonance.h"
#include "MassBin.h"

#include "../miFit.h"

#include <Attributes.h>
#include <AmplitudeBasis.h>
#include <BreitWigner.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MathUtilities.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <SpinAmplitudeCache.h>
#include <ZemachFormalism.h>

#include <BAT/BCGaussianPrior.h>
//#include <BAT/BCSplitGaussianPrior.h>

#include <cassert>
#include <complex>
#include <iterator>
#include <memory>
#include <thread>


template <typename T>
inline std::unique_ptr<yap::Model> make_model()
{ return std::make_unique<yap::Model>(std::make_unique<T>()); }

//std::unique_ptr<yap::Model> d3pi_binned() {
//    using namespace std;
//    using namespace yap;
//
//    // Make a new model.
//    auto M = make_model<ZemachFormalism>();
//
//    auto F = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");
//
//    // Final state particles.
//    auto piPlus  = FinalStateParticle::create(F[ 211]);
//    auto piMinus = FinalStateParticle::create(F[-211]);
//
//    M->setFinalState(piPlus, piMinus, piPlus);
//
//    // Initial state particle.
//    auto D = DecayingParticle::create(F["D+"], radial_size);
//
//    // Partition the mass range.
//    auto p = partition_mass_axis(std::abs(piPlus->mass() + piMinus->mass()),
//                                 std::abs(F["D+"].mass() - piPlus->mass()));
//
//    // The number of bins in the partition.
//    const auto bins = p.size() - 1;
//
//    // vector of resonance pointers
//    for (unsigned i = 0; i < bins; ++i) {
//        const auto resonance_name = "Bin(" + to_string(i) + ")";
//        auto r = DecayingParticle::create(resonance_name, QuantumNumbers(0, 0), radial_size, make_shared<MassBin>(p[i], p[i+1]));
//        r->addWeakDecay(piMinus, piPlus);
//        D->addWeakDecay(r, piPlus);
//    }
//
//    {
//        // Add the rho0 resonance.
//        auto r = DecayingParticle::create(F["rho0"], radial_size, std::make_shared<BreitWigner>(F["rho0"]));
//        r->addWeakDecay(piMinus, piPlus);
//        D->addWeakDecay(r, piPlus);
//
//        // Set this free amplitude as fixed.
//        free_amplitude(*D, to(r))->variableStatus() = VariableStatus::fixed;
//        // Fix the amplitude to its default value, which is the complex (1, 0).
//    }
//
//    M->addInitialState(D);
//
//    // Check that each bin has an associated free amplitude
//    assert(bins == free_amplitudes(*M, from(D), is_not_fixed()).size());
//
//    constexpr double bw_M  = 1.;
//    constexpr double Gamma = .2;
//    auto bw = [=](double s) { return 1. / std::complex<double>(bw_M * bw_M - s, - bw_M * Gamma ); };
//    int i = 0;
//    for (const auto& fa : free_amplitudes(*M, from(*D), is_not_fixed())) {
//
//        // evaluate BW on the low edge of the bin
//        const auto bw_value = bw(p[i] * p[i]);
////        std::cout << "[" << p[i] << p[i + 1] << ") BW("
////                  << p[i] * p[i] << ") = " << bw_value << std::endl;
//
//        *fa = polar(abs(bw_value), arg(bw_value));
////        std::cout << "i = " << i << std::endl;
////        std::cout << to_string(*fa) << std::endl;
//
//        ++i;
//    }
//
//    return M;
//}

//std::unique_ptr<yap::Model> d3pi_one_resonance(std::unique_ptr<yap::Model> M)
//{
//    using namespace std;
//    using namespace yap;
//
//	// Open the particle table.
//    auto F = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");
//
//    // Create the final state particles.
//    auto piPlus  = FinalStateParticle::create(F[ 211]);
//    auto piMinus = FinalStateParticle::create(F[-211]);
//
//	// Set the final state of the decay.
//    M->setFinalState(piPlus, piMinus, piPlus);
//
//    // Create the initial state particle.
//    auto D = DecayingParticle::create(F["D+"], radial_size);
//
//	// Vector of the intermediate resonance names.
//	const std::vector<std::string> resonance_names({"f_0", "rho0"});
//	// Create the decay tree based on the previous resonance names.
//	for (const auto& rn : resonance_names) {
//		// Create the resonance.
//		auto r = DecayingParticle::create(F[rn], radial_size, make_shared<BreitWigner>(F[rn]));
//
//		// Set the resonance daughter particles.
//		r->addWeakDecay(piPlus, piMinus);
//
//		// Broaden the f_0's width.
//		if (rn == "f_0") {
//			auto w = static_pointer_cast<BreitWigner>(r->massShape())->width();
//			w->setValue(4. * w->value());
////			std::cout << "Width " << w->value() << std::endl;
//		}
//
//		// Add the resonance to the D decay.
//		D->addWeakDecay(r, piPlus);
//		std::cout << rn << " resonance added to the D decay." << std::endl;
//	}
//
//	// Set the decay initial state.
//    M->addInitialState(D);
//
//    return M;
//}
//
//constexpr const double pi() noexcept
//{ return std::acos(-1); }


miFit d3pi_binned_fit(const std::string& name)
{
    using namespace std;
    using namespace yap;

    constexpr unsigned integration_points     = 2e4;
    const     unsigned integration_threads    = available_threads();
    const     unsigned integration_batch_size = integration_points / integration_threads;

    miFit m(name, d3pi_binned(), integration_points, integration_threads, integration_batch_size);


    // Lock the model in order to get the free amplitudes only!
    m.model()->lock();

    for (const auto& fa : free_amplitudes(*m.model(), is_not_fixed()))
        m.setRanges(fa, 0., 15. * abs(fa->value()), 0., 180.);

    // Fix an amplitude
    // Use the scope to have local variables
    { 
        const auto fas = free_amplitudes(*m.model(), is_not_fixed());
        const auto& fa = prev(fas.cend(), 1);
        assert(fa != std::end(fas));

        // Fix the amplitude
        m.fixComponent(*fa, 0, 1.);
        // Fix the first bin's phase to 0
        m.fixComponent(*begin(fas), 1, 0.);
    }

    return m;
}
