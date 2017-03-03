/**
 *
 *    @file  FreedWave.h
 *   @brief  
 *
 *    @date  02/26/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  FREED_WAVE_H
#define  FREED_WAVE_H

#include "fwd/MassRangePartition.h"

#include <fwd/DecayingParticle.h>
#include <fwd/QuantumNumbers.h>

#include <memory>
#include <string>
#include <vector>

class FreedWave {
public:

    /// @brief Constructor.
    /// @param mrp      The mass-range partitioning.
    /// @param qn       The resonance quantum numbers.
    /// @param bin_name The (optional) name of the resonance.
    explicit FreedWave(std::unique_ptr<const MassRangePartition> mrp,
                       const yap::QuantumNumbers& qn,
                       const std::string bin_name = "Bin");

    /// _Deleted_ copy constructor.
    FreedWave(const FreedWave&)            = delete;
    /// _Deleted_ copy assignment operator.
    FreedWave& operator=(const FreedWave&) = delete;

    /// _Default_ move constructor.
    FreedWave(FreedWave&&)            = default;
    /// _Default_ move assignment operator.
    FreedWave& operator=(FreedWave&&) = default;

    /// _Default_ destructor.
    virtual ~FreedWave();

    /// Return the quantum numbers of the wave.
    const std::unique_ptr<const yap::QuantumNumbers>& quantumNumbers() const noexcept
    { return QuantumNumbers_; }

    /// Return the mass-range partition.
    const std::unique_ptr<const MassRangePartition>& massRangePartition() const noexcept
    { return MassRangePartition_; }

    /// Return the name of the mass bins.
    const std::string& binName() const noexcept
    { return BinName_; }

    /// Return the bins of the freed wave.
    const std::vector<std::shared_ptr<const yap::DecayingParticle>>& massBins() const noexcept
    { return MassBins_; }

private:

    /// The wave quantum numbers.
    std::unique_ptr<const yap::QuantumNumbers> QuantumNumbers_;

    /// The partitioning of the mass range.
    std::unique_ptr<const MassRangePartition> MassRangePartition_;

    /// The name of the bins.
    std::string BinName_;

    /// The mass bins.
    std::vector<std::shared_ptr<const yap::DecayingParticle>> MassBins_;
};

#endif
