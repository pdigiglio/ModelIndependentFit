/**
 *
 *    @file  KnownMassShapeFreedWave.cxx
 *   @brief  
 *
 *    @date  02/28/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "KnownMassShapeFreedWave.h"

#include "MassRangePartition.h"

#include <QuantumNumbers.h>

#include <cassert>
#include <complex>
#include <functional>
#include <memory>
#include <string>

KnownMassShapeFreedWave::KnownMassShapeFreedWave(std::unique_ptr<const MassRangePartition> mrp,
                                                 const yap::QuantumNumbers& qn,
                                                 std::function<std::complex<double>(double)> ms,
                                                 const std::string bin_name) :
    FreedWave(std::move(mrp), qn, bin_name),
    MassShape_(ms)
{
    assert(MassShape_);
}

KnownMassShapeFreedWave::~KnownMassShapeFreedWave() = default;
