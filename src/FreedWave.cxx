/**
 *
 *    @file  FreedWave.cxx
 *   @brief  
 *
 *    @date  02/27/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "FreedWave.h"

#include "FitModel.h"
#include "MassBin.h"

#include <DecayingParticle.h>
#include <MassRangePartition.h>
#include <QuantumNumbers.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

std::vector<std::shared_ptr<const yap::DecayingParticle>> create_bin_resonances(const MassRangePartition& mrp,
                                                                                const yap::QuantumNumbers& qn,
                                                                                const std::string& bin_name) {
    std::vector<std::shared_ptr<const yap::DecayingParticle>> r;
    r.reserve(mrp.numberOfBins());

    {
        const auto& p = mrp.massPartition();

        // The function to create the bin resonance.
        const auto create_bin = [&](const double low_edge, const double up_edge) {
            const auto name = bin_name + "(" + std::to_string(r.size()) + ")";
            return yap::DecayingParticle::create(name, qn, FitModel::radialSize(),
                                                 std::make_shared<MassBin>(low_edge, up_edge));
        };

        std::transform(std::begin(p), std::end(p) - 1, std::begin(p) + 1, std::back_inserter(r), create_bin);
    }

    assert(mrp.numberOfBins() == r.size());
    return r;
}

FreedWave::FreedWave(std::unique_ptr<const MassRangePartition> mrp,
                     const yap::QuantumNumbers& qn,
                     const std::string bin_name) :
    QuantumNumbers_(std::make_unique<const yap::QuantumNumbers>(qn)),
    MassRangePartition_(std::move(mrp)),
    BinName_(std::move(bin_name)),
    MassBins_(create_bin_resonances(*massRangePartition(), *quantumNumbers(), binName()))
{}

FreedWave::~FreedWave() = default;
