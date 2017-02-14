/**
 *
 *    @file  FittedFreeAmplitude.h
 *   @brief  
 *
 *    @date  02/13/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  FITTED_FREE_AMPLITUDE_H
#define  FITTED_FREE_AMPLITUDE_H

#include <iosfwd>

/// Structure containing the values of a line of the file storing the fitted parameters.
struct FittedFreeAmplitude {
    double bin_low_edge;
    double amplitude;
    double amplitude_uncertainty;
    double phase_difference;
    double phase_difference_uncertainty;
    double cumulative_phase;
    double cumulative_phase_uncertainty;

    constexpr explicit
    FittedFreeAmplitude(double ble, double a,
                        double au,  double pd,
                        double pdu, double cp,
                        double cpu) noexcept :
        bin_low_edge(ble),
        amplitude(a),
        amplitude_uncertainty(au),
        phase_difference(pd),
        phase_difference_uncertainty(pdu),
        cumulative_phase(cp),
        cumulative_phase_uncertainty(cpu)
    {}

    /// (const) Arrow operator (needed to use `FitResultFileIterator::operator->()`)
    constexpr const FittedFreeAmplitude* operator->() const noexcept
    { return this; }

    /// Returns a header for the output fit-result file.
    static std::string header() noexcept;
};

/// Output stream iterator for the FittedFreeAmplitude.
std::ostream& operator<<(std::ostream& s, const FittedFreeAmplitude& ffa);

#endif
