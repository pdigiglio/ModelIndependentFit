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

#include "fwd/Fit.h"
#include "fwd/ModelIndependentFitModel.h"

#include <iosfwd>
#include <memory>
#include <vector>

/// @brief Returns a binned model for a \f$D^+ \to \pi^+\pi^-\pi^+\f$ decay.
/// @param model_name The (optional) name of the model.
std::unique_ptr<ModelIndependentFitModel> binned_d3pi(const std::string model_name = "");

/// @brief Returns a binned model for a \f$D^+ \to \pi^+\pi^-\pi^+\f$ decay.
/// @param file_name  The name of the file from which to take the bin amplitudes.
/// @param model_name The (optional) name of the model.
std::unique_ptr<ModelIndependentFitModel> binned_d3pi_from_file(const std::string& file_name,
                                                                const std::string model_name = "");

/// @brief Attempts an initial guess for the fit parameters based on the model used to generate MC data.
/// @param m (What is it for?)
const std::vector<double> guess_parameters(Fit& m);

/// @brief Prints the fit parameters in an output text file.
/// @param m The Fit model whose parameters to print.
void write_fit_result_to_file(Fit& m);

#endif
