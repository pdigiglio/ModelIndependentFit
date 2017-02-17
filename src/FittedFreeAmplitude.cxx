/**
 *
 *    @file  FittedFreeAmplitude.cxx
 *   @brief  
 *
 *    @date  02/13/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "FittedFreeAmplitude.h"

#include <iomanip>
#include <ostream>

std::string FittedFreeAmplitude::header() noexcept {
    std::ostringstream s;

    s << std::left
      << std::setw(16) << "#bin_low_edge "
      << std::setw(16) << "amplitude "
      << std::setw(16) << "d_amplitude "
      << std::setw(16) << "phase_diff "
      << std::setw(16) << "d_phase_diff "
      << std::setw(16) << "phase "
      << std::setw(16) << "d_phase ";

    return s.str();
}

std::string FittedFreeAmplitude::footer(double upper_edge) noexcept {
    std::ostringstream s;
    
    s << std::left
      << std::setw(15) << upper_edge << " "
      << std::setw(15) << 0 << " "
      << std::setw(15) << 0 << " "
      << std::setw(15) << 0 << " "
      << std::setw(15) << 0 << " "
      << std::setw(15) << 0 << " "
      << std::setw(15) << 0;

    return s.str();
}

std::ostream& operator<<(std::ostream& s, const FittedFreeAmplitude& ffa) {
    s << std::left
      << std::setw(15) << ffa.bin_low_edge << " "
      << std::setw(15) << ffa.amplitude << " "
      << std::setw(15) << ffa.amplitude_uncertainty << " "
      << std::setw(15) << ffa.phase_difference << " "
      << std::setw(15) << ffa.phase_difference_uncertainty << " "
      << std::setw(15) << ffa.cumulative_phase << " "
      << std::setw(15) << ffa.cumulative_phase_uncertainty;
    return s;
}
