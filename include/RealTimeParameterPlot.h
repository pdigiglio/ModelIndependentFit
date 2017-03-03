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
#include <string>
#include <utility>
#include <vector>

/// Class for real-time visualization of the fit parameters against the guessed ones.
class RealTimeParameterPlot {
public:
    /// @brief A pair of a complex vector of the values of a function, \f$f = (f(m_1), \dots, f(m_n))\f$,
    /// and the points on which the function values were sampled, \f$m = (m_1, \dots,m_n)\f$.
    using binned_complex_function = std::pair<std::vector<std::complex<double>>, std::vector<double>>;

    /// @brief Constructor.
    /// @param step_mass_shape The initial mass shape to initialize the plots to.
    /// @param wave_name       The wave name that will be prepended to the canvas name to make it unique.
    explicit RealTimeParameterPlot(binned_complex_function step_mass_shape,
                                   const std::string& wave_name) noexcept;

    /// _Default_ destructor.
    ~RealTimeParameterPlot();

    /// _Default_ move constructor.
    RealTimeParameterPlot(RealTimeParameterPlot&&) = default;

    /// @brief Updates the graphs with the current patameters.
    /// @param p The parameter vector.
    void setCurrentParameters(const std::vector<std::complex<double>>& p) const noexcept;

private:
    /// The unique name for the amplitude TCanvas.
    std::string AmplitudeCanvasName_;
    /// The canvas for the amplitude plot.
    std::unique_ptr<TCanvas> AmplitudeCanvas_;
    /// Plot the amplitudes of the fit parameters.
    std::unique_ptr<TGraph> FitAmplitude_;
    /// Plot the guessed amplitudes.
    std::unique_ptr<TGraph> GuessAmplitude_;

    /// The unique name for the amplitude TCanvas.
    std::string PhaseCanvasName_;
    /// The canvas for the phase plot.
    std::unique_ptr<TCanvas> PhaseCanvas_;
    /// Plot the phases of the fit parameters.
    std::unique_ptr<TGraph> FitPhase_;
    /// Plot the guessed phases.
    std::unique_ptr<TGraph> GuessPhase_;
};

#endif
