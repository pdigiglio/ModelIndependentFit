/**
 *
 *    @file  fit_utilities.h
 *   @brief  
 *
 *    @date  03/01/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  FIT_UTILITIES_H
#define  FIT_UTILITIES_H

#include "Fit.h"

#include "parameter_handling.h"
#include "model/model_independent_d3pi.h"

#include "MassRangePartition.h"
#include "FreedWave.h"
#include "RTPlotMIFitModel.h"
#include "FittedFreeAmplitude.h"
#include "parameter_handling.h"
#include "RootFileHandler.h"
#include "RootFitData.h"
#include "CachedIntegrator.h"

#include <FreeAmplitude.h>

#include <cassert>
#include <memory>
#include <vector>

/// @brief Returns the name of the log file.
/// @param fit The Fit class whose log file name to return.
/// @tparam FitModelType The fit-model type.
template <typename FitModelType>
std::string log_file_name(const Fit<FitModelType>& fit) noexcept
{ return fit.rootFitData()->path() + fit.fitModel()->name() + "_log.txt"; }

/// @brief Helper function to create a Fit class.
/// @param file_path  The directory of the input ROOT file.
/// @param file_name  The name of the input ROOT file.
/// @param model_name The name of the model.
/// @tparam FitModelType The fit-model type.
template <typename FitModelType>
std::unique_ptr<Fit<FitModelType>> create_fit(const char* file_path, const char* file_name, const char* model_name) {
    using Integrator = CachedIntegrator;

    // Create a model to fit the data.
    const std::shared_ptr<FitModelType> fit_model(binned_d3pi(model_name));

    // Create the integrator.
    constexpr unsigned integration_points = 2e4;
    auto integrator(std::make_unique<Integrator>(std::static_pointer_cast<const FitModelType>(fit_model), integration_points));

    auto root_file_handler(std::make_unique<RootFileHandler>(file_path, file_name, model_name));
    auto root_fit_data(std::make_unique<RootFitData>(std::move(root_file_handler), fit_model));

    // Create the BAT model for fitting the data.
    auto fit(std::make_unique<Fit<FitModelType>>(fit_model, std::move(root_fit_data), std::move(integrator)));

    // Fix amplitudes in the fit
    {
        // Get the non-fixed free amplitudes.
        const auto fas = fit_model->freeAmplitudes();
//        // Iterator to the last (valid) one.
//        const auto fa  = std::prev(std::end(fas), 1);
//
//        // Fix the last amplitude.
//        fit->fixAmplitude(*fa, 1);

        const auto guess_params = yap_to_fit_parameters(fas, fit->fitModel());
        // Fix the first phase.
        fit->fixPhase(*begin(fas), guess_params[1]);
    }

    return fit;
}

/// @brief Double-check the parameters of the fit model in _m_ and print the initial (guessed) values to a file.
/// @return A vector of parameter for the `Fit` class.
/// @param m The `Fit` class for which to generate the initial guess.
/// @tparam FitModelType The fit-model type.
template <typename FitModelType>
std::vector<double> guess_parameters(Fit<FitModelType>& m) {
    // Check the type of the fit model.
    static_assert(std::is_base_of<ModelIndependentFitModel, FitModelType>::value,
                  "FitModelType doesn't inherit from FitModelType");
    const auto fit_model = std::static_pointer_cast<const ModelIndependentFitModel>(m.fitModel());

    // Make sure that the parameters of the model match these.
    const std::vector<std::vector<std::complex<double>>> yp({
            binned_mass_shape(D_wave_mass_shape(), *fit_model->freedWaves()[0]->massRangePartition()),
            binned_mass_shape(S_wave_mass_shape(), *fit_model->freedWaves()[1]->massRangePartition()) });
    assert(yp.size() == fit_model->freedWaves().size());
    assert(check_fas_against_parameters(fit_model->freeAmplitudes(), yp));

    // Write the initial parameters of the model to file.
    std::ofstream par_guess("output/par_guess.txt", std::ios::out);
    for (size_t j = 0; j < yp.size(); ++ j) {
        const auto fp = yap_to_fit_parameters(yp[j]);
        assert(2 * yp[j].size() == fp.size());

        par_guess << "# " << fit_model->freedWaves()[j]->binName() << std::endl;
        par_guess << FittedFreeAmplitude::header() << std::endl;
        
        // Cumulative phase.
        double cp = 0.;
        for (size_t i = 0; i < yp[j].size(); ++ i) {
            // Bin low edge.
            const auto ble = fit_model->freedWaves()[j]->massRangePartition()->massPartition().at(i);
            // Free-amplitude amplitude.
            const auto a   = fp.at(2 * i);
            // Free-amplitude phase difference.
            const auto pd  = fp.at(2 * i + 1);
            cp += pd;

            par_guess << FittedFreeAmplitude(ble, a, 0, pd, 0, cp, 0) << std::endl;
        }

        // Leave an empty line between the plots.
        par_guess << std::endl;
    }

    return yap_to_fit_parameters(yp);
}

/// @brief Prints the fit parameters in an output text file.
/// @param m The Fit model whose parameters to print.
/// @tparam FitModelType The fit-model type.
template <typename FitModelType>
void write_fit_result_to_file(Fit<FitModelType>& m) {
    // Check the type of the fit model.
    static_assert(std::is_base_of<ModelIndependentFitModel, FitModelType>::value,
                  "FitModelType doesn't inherit from FitModelType");
    const auto fit_model = std::static_pointer_cast<const ModelIndependentFitModel>(m.fitModel());

    // Check that #parameters == 4 * #bins (i.e. 4 * #FreeAmplitude's).
    // There are 2 parameters + 2 observables per bin.
    assert(m.GetBestFitParameters().size() == 
           std::accumulate(std::begin(fit_model->freedWaves()), std::end(fit_model->freedWaves()), size_t(0),
               [](const size_t s, const auto& fw) { return s + 4 * fw->massRangePartition()->numberOfBins(); }));

    // Print the wave's best parameters to the output file.
    std::ofstream par_fit("output/par_fit.txt", std::ios::out);

    size_t offset = 0;
    for (const auto& fw : fit_model->freedWaves()) {
        par_fit << "# " << fw->binName() << std::endl;
        par_fit << FittedFreeAmplitude::header() << std::endl;

        // Cache the mass partitioning.
        const auto& p = fw->massRangePartition()->massPartition();

        // Cumulative phase.
        double cp  = 0.;
        // Cumulative-phase variance.
        double cpv = 0.;
        for (size_t j = 0; j < fw->massRangePartition()->numberOfBins(); ++ j) {
            // Shift the index to get the proper parameters.
            const auto i = j + offset;

            const auto amp   = m.GetBestFitParameters()[2 * i];
            const auto d_amp = m.GetBestFitParameterErrors()[2 * i];
            const auto phi   = m.GetBestFitParameters()[2 * i + 1];
            const auto d_phi = m.GetBestFitParameterErrors()[2 * i + 1];

            cp  += phi;
            cpv += d_phi * d_phi;

            // Construct a temporary FittedFreeAmplitude to easily stream it.
            par_fit << FittedFreeAmplitude(p[j], amp, d_amp, phi, d_phi, cp, std::sqrt(cpv)) << std::endl;
        }

        // Separate the plots with an empty line
        par_fit << std::endl;

        // Update the offset.
        offset += fw->massRangePartition()->numberOfBins();
    }

    // Check that all the parameters (w/o the observables) have been printed.
    assert(4 * offset == m.GetBestFitParameters().size());
}

#endif
