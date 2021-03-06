/**
 *
 *    @file  MassRangePartition.h
 *   @brief  
 *
 *    @date  02/23/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  MASS_RANGE_PARTITION_H
#define  MASS_RANGE_PARTITION_H

#include <complex>
#include <functional>
#include <utility>
#include <vector>

/// @brief Class to partition a mass range.
class MassRangePartition {
public:
    /// @brief The region type.
    /// @details The regions are a pair of (_initial point_, _number of bins_). The upper boundary
    /// of the partition has 0 bins so that (_point_, _0_) will be considered the upper boundary.
    /// No other region is allowed to have 0 bins.
    using region = std::pair<double, unsigned>;

    /// @brief Constructor.
    /// @param regions The regions in which the mass range will be divided.
    MassRangePartition(const std::vector<region>& regions) noexcept;

    /// Return the low bins of the mass partition.
    const std::vector<double>& massPartition() const noexcept
    { return MassPartition_; }

    /// Return the mass-range regions.
    const std::vector<region> regions() const noexcept
    { return Regions_; }

    /// Return the number of bins.
    unsigned numberOfBins() const noexcept
    { return NumberOfBins_; }

private:
    /// The regions in which the mass range has been divided.
    std::vector<region> Regions_;

    /// The number of bins.
    unsigned NumberOfBins_;

    /// The partition of the Dalitz-plot mass range.
    std::vector<double> MassPartition_;
};

/// @brief Returns a vector of point evaluated according the formula \f$p_i = b_i + x(b_{i+1} - b_i)\f$, being \f$b_i\f$ the _i_-th bin low edge.
/// @param x   The displacement (relative to the bin width) from the bin low edge.
/// @param mrp The mass-range partition whose binning to use.
/// @details If, for example, _x = 0.5_ the bin centers will be returned.
/// @throw yap::Exception If _x_ is not in the range [_0_, _1_].
std::vector<double> bins(const double x, const MassRangePartition& mrp);


/// @brief Returns a vector with the bin centers. 
/// @param mrp The mass-range partition whose binning to use.
inline const std::vector<double> bin_centers(const MassRangePartition& mrp) noexcept
{ return bins(.5, mrp); }

/// @brief Return the values of the mass shape on samples evaluated from the mass partition.
/// @param mass_shape           The mass shape one wants to evaluate.
/// @param mass_range_partition The mass partition to sample the mass shape.
/// @param samples_per_bin      The number of samples per bin.
/// @return A pair of the mass-shape values and the corresponding mass samples.
std::pair<std::vector<std::complex<double>>, std::vector<double>>
binned_mass_shape(std::function<std::complex<double>(double)> mass_shape,
                  const MassRangePartition& mass_range_partition,
                  size_t samples_per_bin);

/// @brief Return the values of the mass shape on the low edges of the partition.
/// @param mass_shape The mass shape one wants to evaluate.
/// @param mass_range_partition The mass partition to sample the mass shape.
inline std::vector<std::complex<double>> binned_mass_shape(std::function<std::complex<double>(double)> mass_shape,
                                                           const MassRangePartition& mass_range_partition)
{ return binned_mass_shape(mass_shape, mass_range_partition, 1).first; }

#endif
