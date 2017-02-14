/**
 *
 *    @file  DataGenerator.cxx
 *   @brief  
 *
 *    @date  02/10/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "DataGenerator.h"

#include "FitModel.h"

#include <CalculationStatus.h>
#include <DataSet.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <VariableStatus.h>

#include <cassert>
#include <iostream>
#include <iterator>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

DataGenerator::DataGenerator(std::shared_ptr<const FitModel> fit_model,
                             const std::string model_name) :
    BCModel(model_name),
    FitModel_(fit_model),
    Data_(GetNChains(), FitModel_->model()->createDataSet(1)) // XXX why?
{
//    // Check that there are no non-fixed free amplitudes.
//    std::cout << "Non-fixed free amplitudes: " << FitModel_->freeAmplitudes().size() << std::endl;
//    assert(FitModel_->freeAmplitudes().empty());

    const auto& m2r = FitModel_->squaredMassRanges();
    const std::string mass_unit = "[GeV/c^2]";
    for (size_t i = 0; i < FitModel_->massAxes().size(); ++ i) {
        const auto axis = FitModel_->massAxes()[i];

        // Evaluate the name of the axis.
        const auto axis_label = "m2_" + yap::indices_string(*axis).substr(1, 2);

        // Add a new BAT parameter for this mass axis.
        AddParameter(axis_label, m2r[i][0], m2r[i][1], axis_label, mass_unit);

#ifndef NDEBUG
        std::cout << "Added parameter: " << axis_label
                  << " with range: [" << m2r[i][0] 
                  << ", " << m2r[i][1] << "]" << std::endl;
#endif
    }
}

double DataGenerator::LogLikelihood(const std::vector<double>& p) {
    const auto c = GetCurrentChain();
    // XXX why?
    return yap::sum_of_log_intensity(*FitModel_->model(), Data_[c]);
}

double DataGenerator::LogAPrioriProbability(const std::vector<double>& p) {
    // Calculate the four momenta of the final-state particles.
    const auto P = yap::calculate_four_momenta(FitModel::Dmass(), *FitModel_->model(), FitModel_->massAxes(), p);

    // If it's outside the phase space, log(prob) = log(0) = - infinity
    if (P.empty())
        return -std::numeric_limits<double>::infinity();

    // XXX why is that?
    const unsigned c = GetCurrentChain();
    Data_[c].setAll(yap::VariableStatus::changed);
    Data_[c].setAll(yap::CalculationStatus::uncalculated);
    FitModel_->model()->setFinalStateMomenta(Data_[c][0], P, Data_[c]);

    return 0; // log(prob) = log(1)
}
