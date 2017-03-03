/**
 *
 *    @file  KnownMassShapeFreedWave.h
 *   @brief  
 *
 *    @date  02/28/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  KNOWN_MASS_SHAPE_FREED_WAVE_H
#define  KNOWN_MASS_SHAPE_FREED_WAVE_H

#include "FreedWave.h"

#include <complex>
#include <functional>

/// @brief Freed wave whose initial mass shape is known (i.e. because the wave content is known).
class KnownMassShapeFreedWave final : public FreedWave {
public:

    /// @brief Constructor.
    /// @param mrp      The mass-range partitioning.
    /// @param qn       The resonance quantum numbers.
    /// @param ms       The freed-wave expected mass shape.
    /// @param bin_name The (optional) name of the resonance.
    explicit KnownMassShapeFreedWave(std::unique_ptr<const MassRangePartition> mrp,
                                     const yap::QuantumNumbers& qn,
                                     std::function<std::complex<double>(double)> ms,
                                     const std::string bin_name = "Bin");

    /// _Default_ destructor.
    ~KnownMassShapeFreedWave() final;

    /// Return the freed-wave mass shape.
    const std::function<std::complex<double>(double)> massShape() const noexcept
    { return MassShape_; }

private:

    /// The expected mass shape of the freed wave.
    std::function<std::complex<double>(double)> MassShape_;
};

#endif
