/**
 *
 *    @file  DataGenerator.h
 *   @brief  
 *
 *    @date  02/10/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  DATA_GENERATOR_H
#define  DATA_GENERATOR_H

#include "fwd/FitModel.h"
#include <fwd/DataSet.h>

#include <BAT/BCModel.h>

#include <memory>
#include <vector>

class DataGenerator final : public BCModel {
public:
    /// @brief Constructor.
    /// @param fit_model The model according to which to generate data.
    explicit DataGenerator(std::shared_ptr<const FitModel> fit_model);

    /// @brief The logarithm of the likelihood function.
    /// @param p The (unused) model parameters.
    double LogLikelihood(const std::vector<double>& p) final;

    /// @brief The logarithm of the priors.
    /// @param p 
    double LogAPrioriProbability(const std::vector<double>& p) final;

private:
    /// The model according to which to generate data.
    const std::shared_ptr<const FitModel> FitModel_;

    /// The data sets.
    std::vector<yap::DataSet> Data_;
};

#endif
