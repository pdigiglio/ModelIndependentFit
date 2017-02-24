/**
 *
 *    @file  RealTimeParameterPlot.h
 *   @brief  
 *
 *    @date  02/21/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  REAL_TIME_PARAMETER_PLOT_H
#define  REAL_TIME_PARAMETER_PLOT_H
 
class TCanvas;
class TGraph;

#include <fwd/FreeAmplitude.h>

#include <complex>
#include <memory>
#include <vector>

/// Class for real-time visualization of the fit parameters against the guessed ones.
class RealTimeParameterPlot {
public:
    /// @brief Constructor.
    /// @param bin_low_edges The low edges of the mass bins.
    /// @param init_fas      The initial guess of the free amplitudes.
    RealTimeParameterPlot(const std::vector<double>& bin_low_edges,
                          const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& init_fas) noexcept;

    /// _Default_ destructor.
    ~RealTimeParameterPlot();

    /// _Default_ move constructor.
    RealTimeParameterPlot(RealTimeParameterPlot&&) = default;

    /// @brief Updates the graphs with the current patameters.
    /// @param p The parameter vector.
    void setCurrentParameters(const std::vector<std::complex<double>>& p) const noexcept;

private:
    /// The canvas for the parameter plots.
    std::unique_ptr<TCanvas> AmplitudeCanvas_;

    /// Plot the amplitudes of the fit parameters.
    std::unique_ptr<TGraph> FitAmplitude_;
    /// Plot the guessed amplitudes.
    std::unique_ptr<TGraph> GuessAmplitude_;

    std::unique_ptr<TCanvas> PhaseCanvas_;
    /// Plot the phases of the fit parameters.
    std::unique_ptr<TGraph> FitPhase_;
    /// Plot the guessed phases.
    std::unique_ptr<TGraph> GuessPhase_;
};

#endif
