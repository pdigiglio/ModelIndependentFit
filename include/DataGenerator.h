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

/// Class to generate Monte-Carlo data according to a FitModel.
class DataGenerator final : public BCModel {
public:
    /// @brief Constructor.
    /// @param fit_model  The model according to which to generate data.
    explicit DataGenerator(std::shared_ptr<const FitModel> fit_model);

    /// _Default_ destructor.
    ~DataGenerator() final;

    /// _Deleted_ copy constructor.
    DataGenerator(const DataGenerator&)            = delete;
    /// _Deleted_ copy assignment operator.
    DataGenerator& operator=(const DataGenerator&) = delete;

    /// _Deleted_ move constructor.
    DataGenerator(DataGenerator&&)            = delete;
    /// _Deleted_ move assignment operator.
    DataGenerator& operator=(DataGenerator&&) = delete;

    /// @brief The logarithm of the likelihood function.
    /// @param p The (unused) suqared masses of the final-state particles.
    double LogLikelihood(const std::vector<double>& p) final;

    /// @brief The logarithm of the priors.
    /// @details It's _0 = log(1)_ if _p_ lies in the phase space, _log(0)_ otherwise.
    /// @param p The squared masses of the final-state particles.
    double LogAPrioriProbability(const std::vector<double>& p) final;

private:
    /// The model according to which to generate data.
    const std::shared_ptr<const FitModel> FitModel_;

    /// The data sets.
    std::vector<yap::DataSet> Data_;
};

#endif
