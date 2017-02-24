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
#include "MassRangePartition.h"

#include <Attributes.h>
#include <DecayingParticle.h>
#include <FreeAmplitude.h>

#include <cassert>
#include <functional>
#include <memory>
#include <string>
#include <vector>

ModelIndependentFitModel::ModelIndependentFitModel(std::unique_ptr<yap::Model> model,
                                                   std::unique_ptr<const MassRangePartition> mass_partition,
                                                   const std::string model_name) :
    FitModel(std::move(model), model_name),
    MassRangePartition_(std::move(mass_partition))
{
    sortBinFreeAmplitudes();
}

ModelIndependentFitModel::~ModelIndependentFitModel() = default;

const std::vector<double>& ModelIndependentFitModel::massPartition() const noexcept
{ return massRangePartition()->massPartition(); }

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

const bool ModelIndependentFitModel::massBinSorted() const {
    // Check that the free amplitudes are as many as the bins.
    assert(massRangePartition()->massPartition().size() == freeAmplitudes().size() + 1);

    for (const auto& fa : freeAmplitudes()) {
        // Get the position of the free amplitude in the free-amplitude vector.
        const auto i = free_amplitude_index(fa, *this);

        // Get the mass shape of the free amplitude.
        const auto mass_bin = bin_mass_shape(fa, *this);

//        assert(mass_bin->lowerEdge()->value() == massPartition()[i]);
        if (mass_bin->lowerEdge()->value() != massRangePartition()->massPartition()[i])
            return false;

//        assert(mass_bin->upperEdge()->value() == massPartition()[i+1]);
        if (mass_bin->upperEdge()->value() != massRangePartition()->massPartition()[i+1])
            return false;
    }

    return true;
}

void ModelIndependentFitModel::sortBinFreeAmplitudes() {
    std::sort(std::begin(accessFreeAmplitudes()), std::end(accessFreeAmplitudes()),
              [](const auto& lhs, const auto& rhs) { return bin_number(lhs) < bin_number(rhs); } );
}

std::shared_ptr<const MassBin> bin_mass_shape(const std::shared_ptr<const yap::FreeAmplitude>& fa,
                                              const ModelIndependentFitModel& M) {

    // Get the position of the free amplitude in the free-amplitude vector.
    const auto i = free_amplitude_index(fa, M);

    // Check that the position in the vector corresponds to the bin number.
    if (i != bin_number(fa))
        throw yap::exceptions::Exception("Bin free amplitudes are not sorted", "bin_mass_shape");

    const auto bin_name = "Bin(" + std::to_string(i) + ")";
    const auto bin = std::static_pointer_cast<const yap::DecayingParticle>(yap::particle(*M.model(), yap::is_named(bin_name)));

    return std::static_pointer_cast<const MassBin>(bin->massShape());
}

std::vector<std::complex<double>> binned_mass_shape(std::function<std::complex<double>(double)> mass_shape,
                                                    const ModelIndependentFitModel& fm) {
    std::vector<std::complex<double>> binned_ms;
    binned_ms.reserve(fm.massRangePartition()->massPartition().size());

    // Evaluate the mass shape on the low edges of the bins (exclude the last one).
    std::transform(std::begin(fm.massRangePartition()->massPartition()), --std::end(fm.massRangePartition()->massPartition()), std::back_inserter(binned_ms),
                   [&] (const auto& ble) { return mass_shape(ble * ble); });

    assert(binned_ms.size() == fm.freeAmplitudes().size());
    return binned_ms;
}
