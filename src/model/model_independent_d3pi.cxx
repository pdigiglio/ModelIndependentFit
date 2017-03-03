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
#include "FreedWave.h"
#include "KnownMassShapeFreedWave.h"
#include "MassBin.h"
#include "MassRangePartition.h"
#include "ModelIndependentFitModel.h"
#include "RTPlotMIFitModel.h"
//#include "RealTimeParameterPlot.h"

#include <DecayChannel.h>

#include <Attributes.h>
#include <BreitWigner.h>
#include <DecayTree.h> // for to_string
#include <DecayingParticle.h>
#include <Exceptions.h>
#include <Filter.h>
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
std::unique_ptr<const MassRangePartition> D_wave_partition(double low_range, double up_range) {
    assert(low_range < up_range);

    auto mrp = std::make_unique<const MassRangePartition>(std::vector<MassRangePartition::region>(
                { MassRangePartition::region(low_range, 10),
                  MassRangePartition::region( 1.05,  14),
                  MassRangePartition::region(1.40,   6),
//                  MassRangePartition::region(1.25,   2),
//                  MassRangePartition::region(1.4,   6),
//                  MassRangePartition::region(1.6,   2),
                  MassRangePartition::region(up_range, 0) }));

    return mrp;
}

// Helper function to create a partition with bins that have different widths.
std::unique_ptr<const MassRangePartition> S_wave_partition(double low_range, double up_range) {
    assert(low_range < up_range);

    auto mrp = std::make_unique<const MassRangePartition>(std::vector<MassRangePartition::region>(
                { MassRangePartition::region(low_range, 10),
                  MassRangePartition::region( .85,  14),
                  MassRangePartition::region(1.10,   6),
//                  MassRangePartition::region(1.25,   2),
//                  MassRangePartition::region(1.4,   6),
//                  MassRangePartition::region(1.6,   2),
                  MassRangePartition::region(up_range, 0) }));

    return mrp;
}

// Helper function.
// Fetches a resonance from a model.
const std::shared_ptr<const yap::DecayingParticle> fetch_resonance(const std::string resonance_name,
                                                                   const std::unique_ptr<yap::Model>& m) {
#ifndef NDEBUG
    std::cout << " > Fetching " << resonance_name << " from model" << std::endl;
#endif
    auto res = yap::particle(*m, yap::is_named(resonance_name));
    assert(res != nullptr);
    return std::static_pointer_cast<const yap::DecayingParticle>(res);
}


// Helper function.
// Returns a Breit-Wigner with the resonance parameters.
std::function<std::complex<double>(double)> bw_mass_shape(const std::shared_ptr<const yap::DecayingParticle>& r) {
    const auto bw = std::static_pointer_cast<const yap::BreitWigner>(r->massShape());
    const auto m = bw->mass()->value();
    const auto w = bw->width()->value();

    return [=](double s) { return 1. / std::complex<double>(m * m - s, - m * w); };
}

// Helper function to return the mass shape of the S wave.
template <size_t S>
std::function<std::complex<double>(double)> wave_mass_shape() {
    // Create a model according to the data model.
    const auto data_model = d3pi();
    
    std::vector<std::shared_ptr<yap::Particle>> p;
    {
        // Get the particles with spin S.
        const auto tmp = yap::particles(*data_model, yap::s_equals(2 * S));
        assert(!tmp.empty());

        // Only copy the resonances.
        std::copy_if(std::begin(tmp), std::end(tmp), std::back_inserter(p),
                     [&](const auto& r) {
                         const auto isp = data_model->initialStates();
                         const auto it = std::find(std::begin(isp), std::end(isp), r);
                         return yap::is_decaying_particle(r) && it == std::end(isp);
                     });

        assert(p.size());
    }

#ifndef NDEBUG
    switch (S) {
        case 0:
            std::cerr << "S-wave content:" << std::endl;
            break;
        case 1:
            std::cerr << "P-wave content:" << std::endl;
            break;
        case 2:
            std::cerr << "D-wave content:" << std::endl;
            break;
        case 3:
            std::cerr << "F-wave content:" << std::endl;
            break;
        default:
            std::cerr << "Content of the " << S << "-spin wave:" << std::endl;
            break;

    };
    for (const auto& r : p)
        std::cerr << " > " << yap::to_string(*r) << std::endl;
    
#endif

    std::vector<std::function<std::complex<double>(double)>> mass_shapes;
    mass_shapes.reserve(p.size());

    std::transform(std::begin(p), std::end(p), std::back_inserter(mass_shapes),
                   [] (const auto& r) { return bw_mass_shape(std::static_pointer_cast<const yap::DecayingParticle>(r)); });
    
    return [=] (double s) {
        return std::accumulate(std::begin(mass_shapes), std::end(mass_shapes), std::complex<double>(0., 0.),
                               [=] (const auto& sum, const auto& ms) { return sum + ms(s); });
    };
}

// helper function to return the mass shape of the s wave.
std::function<std::complex<double>(double)> S_wave_mass_shape() {
    return wave_mass_shape<0>();
}

// helper function to return the mass shape of the p wave.
std::function<std::complex<double>(double)> P_wave_mass_shape() {
    return wave_mass_shape<1>();
}

// helper function to return the mass shape of the d wave.
std::function<std::complex<double>(double)> D_wave_mass_shape() {
    return wave_mass_shape<2>();
}

// Helper function to create a binned fit model.
std::shared_ptr<RTPlotMIFitModel> create_d3pi_binned(const std::string& model_name) {
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

    const auto piPlusMass  = piPlus->mass();
    const auto piMinusMass = piMinus->mass();
    const auto DPlusMass   = F["D+"].mass();


    // Function to add a resonance to the model.
    const auto add_resonance_to_model = [&] (const auto& r) {
//#ifndef NDEBUG
//        std::cerr << "Adding resonance: " << yap::to_string(*r) << std::endl;
//        std::cerr << " > Address: " << r << std::endl;
//#endif
        const auto tmp_r = std::const_pointer_cast<yap::DecayingParticle>(r);
        assert(tmp_r == r);

        tmp_r->addWeakDecay(piMinus, piPlus);
        D->addWeakDecay(tmp_r, piPlus);
    };

    std::vector<std::shared_ptr<const KnownMassShapeFreedWave>> fws;
    {
        auto mrp = D_wave_partition(std::abs(piPlusMass + piMinusMass),
                                    std::abs(DPlusMass  - piPlusMass));

        const auto mass_shape = D_wave_mass_shape();
        auto fw(std::make_shared<KnownMassShapeFreedWave>(std::move(mrp), QuantumNumbers(0, 4), mass_shape, "D-Bin"));

        // Add the resonances to the model.
        assert(!fw->massBins().empty());
        std::for_each(std::begin(fw->massBins()), std::end(fw->massBins()), add_resonance_to_model);

        // Add the freed wave to the vector
        fws.emplace_back(std::static_pointer_cast<const KnownMassShapeFreedWave>(fw));
    }

    {
        auto mrp = S_wave_partition(std::abs(piPlusMass + piMinusMass),
                                    std::abs(DPlusMass  - piPlusMass));

        const auto mass_shape = S_wave_mass_shape();
        auto fw(std::make_shared<KnownMassShapeFreedWave>(std::move(mrp), QuantumNumbers(0, 0), mass_shape, "S-Bin"));

        // Add the resonances to the model.
        assert(!fw->massBins().empty());
        std::for_each(std::begin(fw->massBins()), std::end(fw->massBins()), add_resonance_to_model);

        // Add the freed wave to the vector
        fws.emplace_back(std::static_pointer_cast<const KnownMassShapeFreedWave>(fw));
    }

//    assert(S_wave_mass_shape());
    {
//        // Add the rho0 resonance.
//        auto rho0 = F["rho0"];
//        // Make sure it's in the model: this will throw if there's no rho0 in the data model.
//        fetch_resonance("rho0");
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

    return std::make_unique<RTPlotMIFitModel>(std::move(M), fws, model_name);
}

std::shared_ptr<RTPlotMIFitModel> binned_d3pi(const std::string model_name) {
    auto M = create_d3pi_binned(model_name);
    assert(M->model()->locked());
    assert(mass_bin_sorted(*M));

    std::vector<std::vector<std::complex<double>>> pars(
            { binned_mass_shape(D_wave_mass_shape(), *M->freedWaves()[0]->massRangePartition()),
              binned_mass_shape(S_wave_mass_shape(), *M->freedWaves()[1]->massRangePartition()) });

    M->setParameters(pars);

    return M;
}

// TODO Adapt this to the new parameter data-type.
// Maybe make a parameter_from_file() function
//std::shared_ptr<ModelIndependentFitModel> binned_d3pi_from_file(const std::string& file_name,
//                                                                const std::string model_name) {
//    // Make a binned model
//    auto M = create_d3pi_binned(model_name);
//    assert(M->model()->locked());
//    assert(mass_bin_sorted(*M));
//
//    // Get the non-fixed free amplitudes.
//    const auto fas = M->freeAmplitudes();
//
//    // Vector of fitted free amplitudes.
//    std::vector<FittedFreeAmplitude> ffa_vector;
//    ffa_vector.reserve(fas.size());
//
//    {
//        // Check that the file can be opened.
//        std::ifstream par_fit(file_name, std::ios::in);
//        if (!par_fit)
//            throw yap::exceptions::Exception("Can't open " + file_name, "binned_d3pi_from_file");
//
//        // Load the entries in the vector.
//        std::copy(FitResultFileIterator(par_fit), FitResultFileIterator::end(),
//                  std::back_inserter(ffa_vector));
//    }
//
//    // Check that the number of lines and the number of free amplitudes match.
//    const auto lines = ffa_vector.size();
//    if (fas.size() != lines) {
//        const std::string error_string = "The number of FreeAmplitudes (" + std::to_string(fas.size()) +
//                                         ") doesn't correspond to the file lines (" + std::to_string(lines) + ")";
//
//        throw yap::exceptions::Exception(error_string, "binned_d3pi_from_file");
//    }
//
////#ifndef NDEBUG
////    {
////        // Check if the bin low edge of the mass bin and the one of the fitted free
////        // amplitude match within a certain tollerance.
////        const auto ble_matches = [] (const auto& mass_bin, const auto& ffa) {
////            return std::abs(mass_bin->lowerEdge()->value() - ffa.bin_low_edge) < 1e-5;
////        };
////
////        // Check if the bin low edges of the file match the ones in the mass bins.
////        assert(std::all_of(std::begin(fas), std::end(fas), [&](const auto& fa)
////               { return ble_matches(bin_mass_shape(fa, *M), ffa_vector[free_amplitude_index(fa, *M)]); }));
////    }
////#endif
//
//    // Get the complex-parameter vector from the input file.
//    std::vector<std::complex<double>> p;
//    p.reserve(fas.size());
//
//    std::transform(std::begin(ffa_vector), std::end(ffa_vector), std::back_inserter(p),
//                   [](const auto& ffa) { return std::polar(ffa.amplitude, yap::rad(ffa.cumulative_phase)); });
//
//    M->fixParameters(p);
//
//    return M;
//}
