/**
 *
 *    @file  parameter_handling.h
 *   @brief  Contains helper functions to convert YAP <-> fit parameters.
 *
 *    @date  03/01/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  PARAMETER_HANDLING_INC
#define  PARAMETER_HANDLING_INC

#include "fwd/ModelIndependentFitModel.h"

#include <fwd/FreeAmplitude.h>

#include <complex>
#include <vector>
#include <memory>

/// @brief Converts the fit parameters into FreeAmplitude's.
/// @details This is the function that actually interprets the fit parameters, as
/// they will have a different meaning depending on the function that is applied
/// to them during the conversion.
/// @param p The BAT fit parameters.
/// @param m The model according to which the parameters will be interpreted.
std::vector<std::vector<std::complex<double>>> fit_to_yap_parameters(const std::vector<double>& p,
                                                                     const std::shared_ptr<const ModelIndependentFitModel>& m) noexcept;

/// @brief Converts the complex free-amplitude parameters into fit amplitudes.
/// @attention This should be the inverse of `fit_to_yap_parameters`.
/// @param p The YAP parameters.
std::vector<double> yap_to_fit_parameters(const std::vector<std::complex<double>>& p) noexcept;

/// @brief Converts the complex free-amplitude parameters into fit amplitudes.
/// @attention This should be the inverse of `fit_to_yap_parameters`.
/// @param p The YAP parameters.
std::vector<double> yap_to_fit_parameters(const std::vector<std::vector<std::complex<double>>>& p) noexcept;

/// @brief Converts the complex free-amplitude parameters into fit amplitudes.
/// @param fas The vector of free amplitudes.
/// @param m The model according to which the parameters will be interpreted.
std::vector<double> yap_to_fit_parameters(const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& fas,
                                          const std::shared_ptr<const ModelIndependentFitModel>& m) noexcept;

/// @brief Helper function to check a vector of FreeAmplitude's agains some parameters.
/// @param fas The vector of FreeAmplitude's to check.
/// @param p   The parameters to check the FreeAmplitude's against.
bool check_fas_against_parameters(const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& fas,
                                  const std::vector<std::vector<std::complex<double>>>& p) noexcept;
#endif
