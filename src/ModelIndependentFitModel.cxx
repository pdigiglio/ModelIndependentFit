/**
 *
 *    @file  ModelIndependentFitModel.cxx
 *   @brief  
 *
 *    @date  02/14/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "ModelIndependentFitModel.h"

#include "FreedWave.h"
#include "MassBin.h"
#include "MassRangePartition.h"

#include <Attributes.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <FreeAmplitude.h>

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>
#include <string>
#include <vector>

// XXX Can this introduce any address shifting due to type cast?
std::shared_ptr<const yap::DecayingParticle> resonance(const std::shared_ptr<const yap::FreeAmplitude>& fa) {
    std::vector<std::shared_ptr<const yap::Particle>> p;

    {
        const auto dc = fa->decayChannel();
        std::copy_if(std::begin(dc->daughters()), std::end(dc->daughters()), std::back_inserter(p),
                [](const auto& d) { return yap::is_decaying_particle(d); });
    }

    // The found element should be unique.
    assert(p.size() == 1);

    return std::static_pointer_cast<const yap::DecayingParticle>(*std::begin(p));
}

// Helper function.
// Returns the index of the free amplitude, using the position of its corresponding MassBin.
// XXX This will throw if the MassBin is not found -> it will throw for a non-binned free amplitude.
size_t free_amplitude_index(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                            const std::vector<std::shared_ptr<const FreedWave>>& fws) {
    // TODO check if the mass shape can be converted into a MassBin.
    const auto r = resonance(fa);
//#ifndef NDEBUG
//    std::cerr << "Looking for FreeAmplitude: " << yap::to_string(*fa) << std::endl;
//    std::cerr << " > Corresponding to particle: " << yap::to_string(*r) << std::endl;
//#endif

    // Look for the resonance.
    size_t offset = 0;
    for (const auto& fw : fws) {
        // Cache the mass bins.
        const auto& mbs = fw->massBins();

        // Search for the resonance in the mass-bin vector.
        const auto it = std::find(std::begin(mbs), std::end(mbs), r);

        // If the element was found, return.
        if (it != std::end(mbs)) {

            // The index in the mass-bin vector.
            const auto i = std::distance(std::begin(mbs), it);

#ifndef NDEBUG
            // -- Make sure that the mass-shape boundaries match the partitioning boundaries -------
            const auto mass_bin = std::dynamic_pointer_cast<const MassBin>((*it)->massShape());
            // Make sure mass_bin is not empty.
            assert(mass_bin.use_count());

            assert(std::abs(mass_bin->lowerEdge()->value() - fw->massRangePartition()->massPartition()[i]) < 1e-8);
            assert(std::abs(mass_bin->upperEdge()->value() - fw->massRangePartition()->massPartition()[i+1]) < 1e-8);
            // -------------------------------------------------------------------------------------
#endif

            return i + offset;
        }

        // Update the offset.
        offset += mbs.size();
    }

    // XXX this will also happen if there's a non-binned free amplitude!
    throw yap::exceptions::Exception("No MassBin found for the free amplitude " + yap::to_string(*fa),
                                     "free_amplitude_index");
}

ModelIndependentFitModel::ModelIndependentFitModel(std::unique_ptr<yap::Model> model,
                                                   const std::vector<std::shared_ptr<const FreedWave>>& freed_waves,
                                                   const std::string model_name) :
    FitModel(std::move(model), model_name),
    FreedWaves_(std::move(freed_waves))
{
    sortBinFreeAmplitudes();
    assert(mass_bin_sorted(*this));
}

ModelIndependentFitModel::~ModelIndependentFitModel() = default;

void ModelIndependentFitModel::sortBinFreeAmplitudes() {
    std::sort(std::begin(accessFreeAmplitudes()), std::end(accessFreeAmplitudes()),
              [&](const auto& lhs, const auto& rhs) { return free_amplitude_index(lhs, this->freedWaves()) < free_amplitude_index(rhs, this->freedWaves()); } );
}

const bool mass_bin_sorted(const ModelIndependentFitModel& m) {
    const auto fa_compare = [&](const auto& lhs, const auto& rhs) {
        return free_amplitude_index(lhs, m.freedWaves()) < free_amplitude_index(rhs, m.freedWaves());
    };

    return std::is_sorted(std::begin(m.freeAmplitudes()), std::end(m.freeAmplitudes()), fa_compare);
}
