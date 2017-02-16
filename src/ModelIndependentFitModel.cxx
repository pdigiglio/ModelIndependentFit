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

#include "MassBin.h"

#include <Attributes.h>
#include <DecayingParticle.h>
#include <FreeAmplitude.h>

#include <cassert>
#include <memory>
#include <string>
#include <vector>

ModelIndependentFitModel::ModelIndependentFitModel(std::unique_ptr<yap::Model> model,
                                                   const std::vector<double>& mass_partition,
                                                   const std::string model_name) :
    FitModel(std::move(model), model_name),
    MassPartition_(mass_partition)
{
    sortBinFreeAmplitudes();
}

// Helper function to retrieve the bin number from the free amplitude name.
const size_t bin_number(const std::shared_ptr<const yap::FreeAmplitude>& fa) {
    const auto fa_name = yap::to_string(*fa);
    assert(fa_name.find("Bin(") != std::string::npos);

    std::string fa_bin_number;
    for (const auto& c : fa_name) {

        // If the character is a number, push it back.
        if (std::isdigit(c))
            fa_bin_number.push_back(c);
        else
            if (!fa_bin_number.empty())
                break;
    }

    // A number should be found.
    assert(!fa_bin_number.empty());

//    std::cout << "Name: " << fa_name << std::endl;
//    std::cout << "Bin:  " << std::stoul(fa_bin_number) << std::endl;

    return std::stoul(fa_bin_number);
}

const bool ModelIndependentFitModel::massBinSorted() const noexcept {
    // Check that the free amplitudes are as many as the bins.
    assert(massPartition().size() == freeAmplitudes().size() + 1);

    for (const auto& fa : freeAmplitudes()) {
        // Get the position of the free amplitude in the free-amplitude vector.
        const auto i = free_amplitude_index(fa, *this);

        // Check that the position in the vector corresponds to the bin number.
//        assert(i == bin_number(fa));
        if (i != bin_number(fa))
            return false;

        const auto bin_name = "Bin(" + std::to_string(i) + ")";
        const auto bin = std::static_pointer_cast<const yap::DecayingParticle>(yap::particle(*model(), yap::is_named(bin_name)));
        const auto mass_bin = std::static_pointer_cast<const MassBin>(bin->massShape());

//        assert(mass_bin->lowerEdge()->value() == massPartition()[i]);
        if (mass_bin->lowerEdge()->value() != massPartition()[i])
            return false;

//        assert(mass_bin->upperEdge()->value() == massPartition()[i+1]);
        if (mass_bin->upperEdge()->value() != massPartition()[i+1])
            return false;
    }

    return true;
}

void ModelIndependentFitModel::sortBinFreeAmplitudes() {
    std::sort(std::begin(freeAmplitudes()), std::end(freeAmplitudes()),
              [](const auto& lhs, const auto& rhs) { return bin_number(lhs) < bin_number(rhs); } );
}
