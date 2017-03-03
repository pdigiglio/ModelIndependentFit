/**
 *
 *    @file  model_independent_d3pi.h
 *   @brief  
 * @details
 *
 *    @date  02/14/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  MODEL_INDEPENDENT_D3PI_INC
#define  MODEL_INDEPENDENT_D3PI_INC

#include "fwd/ModelIndependentFitModel.h"
#include "fwd/RTPlotMIFitModel.h"

#include <complex>
#include <functional>
#include <iosfwd>
#include <memory>
#include <vector>

/// @brief Returns a binned model for a \f$D^+ \to \pi^+\pi^-\pi^+\f$ decay.
/// @param model_name The (optional) name of the model.
std::shared_ptr<RTPlotMIFitModel> binned_d3pi(const std::string model_name = "");

/// @brief Returns a binned model for a \f$D^+ \to \pi^+\pi^-\pi^+\f$ decay.
/// @param file_name  The name of the file from which to take the bin amplitudes.
/// @param model_name The (optional) name of the model.
std::shared_ptr<ModelIndependentFitModel> binned_d3pi_from_file(const std::string& file_name,
                                                                const std::string model_name = "");


/// Returns the _S_ wave of the \f$D^+\to\pi^+\pi^-\pi^+\f$ that was used to generate the MC data.
std::function<std::complex<double>(double)> S_wave_mass_shape();

/// Returns the _P_ wave of the \f$D^+\to\pi^+\pi^-\pi^+\f$ that was used to generate the MC data.
std::function<std::complex<double>(double)> P_wave_mass_shape();

/// Returns the _D_ wave of the \f$D^+\to\pi^+\pi^-\pi^+\f$ that was used to generate the MC data.
std::function<std::complex<double>(double)> D_wave_mass_shape();

#endif
