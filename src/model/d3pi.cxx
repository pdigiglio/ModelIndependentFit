/**
 *
 *    @file  d3pi.cxx
 *   @brief
 *
 *    @date  02/14/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "model/d3pi.h"

#include "FitModel.h"

#include <Attributes.h>
#include <BreitWigner.h>
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
//	const std::vector<std::string> resonance_names({"rho0", "f_0", "f_0(1500)"});
//	const std::vector<std::string> resonance_names({"f_0", "f_0(1500)"});
	const std::vector<std::string> resonance_names({"f_0", "f_2"});
	// Create the decay tree based on the previous resonance names.
	for (const auto& rn : resonance_names) {
		// Create the resonance.
		auto r = DecayingParticle::create(F[rn], FitModel::radialSize(), make_shared<BreitWigner>(F[rn]));

		// Set the resonance daughter particles.
		r->addWeakDecay(piPlus, piMinus);

		// Add the resonance to the D decay.
		D->addWeakDecay(r, piPlus);

#ifndef NDEBUG
		std::cout << rn << " resonance added to the D decay." << std::endl;
#endif
	}

	// Set the decay initial state.
    M->addInitialState(D);
    // Lock the model before returning it.
    M->lock();

    // Fix the free amplitudes.
    const auto fas = yap::free_amplitudes(*M, yap::from(*D), yap::is_not_fixed());
    std::for_each(std::begin(fas), std::end(fas),
                  [](const auto& fa) { fa->variableStatus() = yap::VariableStatus::fixed; });

    return M;
}
