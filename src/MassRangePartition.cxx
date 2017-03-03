/**
 *
 *    @file  MassRangePartition.cxx
 *   @brief  
 *
 *    @date  02/23/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "MassRangePartition.h"

#include <Exceptions.h>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

inline unsigned count_bins(const std::vector<MassRangePartition::region>& regions) noexcept {
    return std::accumulate(std::begin(regions), std::end(regions), 0,
                           [](unsigned a, const auto& b) { return a + b.second; });
}

// Helper function.
// Splits a range _[min, max)_ into a partition with _bins_ number of bins.
std::vector<double> partition_range(double min, const double max, const unsigned int bins) {
    assert(bins > 1);

    // The bin size.
    assert(min < max);
    const double step = (max - min) / bins;

    std::vector<double> p;
    p.reserve(bins + 1);
    std::generate_n(std::back_inserter(p), bins + 1, [&]() { const auto m = min; min += step; return m; });
    assert(p.size() == bins + 1);

    return p;
}

std::vector<double> partition_mass_range(std::vector<MassRangePartition::region> regions) {
    assert(regions.size() > 1);

    // The partition vector.
    std::vector<double> p;
    p.reserve(count_bins(regions) + 1);

    // Make the partition.
    unsigned i;
    for (i = 0; i < regions.size() - 1; ++ i) {
        const auto bins = regions.at(i).second;
        assert(bins != 0);

        const auto low_range = regions.at(i).first;
        const auto up_range  = regions.at(i + 1).first;
        assert(low_range < up_range);

        // Partition the current range.
        auto rp = partition_range(low_range, up_range, bins);
        assert(low_range  == *std::begin(rp));
        assert(std::abs(up_range - *(std::end(rp) - 1)) < 1e-8);
        assert(bins + 1   == rp.size());

        // Remove the last element (which will be the first one in the next iteration).
        rp.pop_back();

        // Move the vector into the partition.
        std::copy(std::make_move_iterator(std::begin(rp)), std::make_move_iterator(std::end(rp)),
                  std::back_inserter(p));
    }

    // Push back the upper boundary of the partition.
    assert(regions[i].second == 0);
    p.push_back(regions[i].first);

    return p;
}

MassRangePartition::MassRangePartition(const std::vector<region>& r) noexcept :
    Regions_(r),
    NumberOfBins_(count_bins(regions())),
    MassPartition_(partition_mass_range(regions()))
{
    // Check if it's sorted.
    assert(std::is_sorted(std::begin(massPartition()), std::end(massPartition())));
    assert(numberOfBins() + 1 == massPartition().size());

#ifndef NDEBUG
    for (size_t i = 0; i < numberOfBins(); ++ i) {
        std::cout << "bin(" << std::setw(2) << std::right << i << ") = ["
                            << std::setw(8) << std::left << massPartition()[i]
                  << " , "  << std::setw(8) << std::left << massPartition()[i+1]
                  << ")" << std::endl;
    }
#endif
}

std::vector<double> bins(const double x, const MassRangePartition& mrp) {
    if (x < 0. or x > 1.)
        throw yap::exceptions::Exception("The parameter x must be in [0,1]", "MassRangePartition::bins");

    std::vector<double> p;
    p.reserve(mrp.numberOfBins());

    {
        const auto& mp = mrp.massPartition();
        std::transform(std::begin(mp), std::end(mp) - 1, std::begin(mp) + 1, std::back_inserter(p),
                       [&](const auto& a, const auto& b) { return a + x * (b - a); });
    }

    return p;
}

std::pair<std::vector<std::complex<double>>, std::vector<double>>
binned_mass_shape(std::function<std::complex<double>(double)> mass_shape,
                  const MassRangePartition& mrp,
                  size_t samples_per_bin) {

    std::vector<std::complex<double>> binned_ms;
    binned_ms.reserve(samples_per_bin * mrp.numberOfBins());

    std::vector<double> sample_points;
    sample_points.reserve(samples_per_bin * mrp.numberOfBins());

    for (size_t i = 0; i < mrp.numberOfBins(); ++ i) {
        // Bin low edge.
        const auto ble = mrp.massPartition()[i];
        // Sampling step.
        const auto ss  = (mrp.massPartition()[i + 1] - ble) / samples_per_bin;

        for (size_t s = 0; s < samples_per_bin; ++ s) {
            sample_points.emplace_back(ble + s * ss);
            binned_ms.emplace_back(mass_shape(sample_points.back() * sample_points.back()));
        }
    }

    return std::make_pair(binned_ms, sample_points);
}
