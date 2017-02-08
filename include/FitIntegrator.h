/**
 *
 *    @file  FitIntegrator.h
 *   @brief  
 *
 *    @date  02/03/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  INTEGRATOR_H
#define  INTEGRATOR_H

#include "fwd/FitModel.h"

#include <fwd/FourVector.h>

#include <DataSet.h>
#include <ModelIntegral.h>

#include <functional>
#include <random>

/// Sample-integrator of the intensity over the model.
class FitIntegrator {
public: 
    /// Function for generating new points for integration.
    using Generator = std::function<std::vector<yap::FourVector<double> >()>;

    /// The type of the random-number engine.
    using RandomEngine = std::mt19937;

    /// @brief Constructor.
    explicit FitIntegrator(std::shared_ptr<const FitModel> fit_model,
                           const unsigned integration_points) noexcept;

    /// Returns the point generator (const version).
    const Generator& pointGenerator() const noexcept
    { return PointGenerator_; }

    /// Access the point generator.
    Generator& pointGenerator() noexcept
    { return PointGenerator_; }

    /// Returns the number of integration points.
    const unsigned numberOfPoints() const noexcept
    { return NumberOfPoints_; }

    const unsigned batchSize() const noexcept
    { return BatchSize_; }

    /// Returns the number of threads used for the integration.
    const unsigned numberOfThreads() const noexcept
    { return NumberOfThreads_; }

    /// Returns the matrix of the model component integrals.
    const yap::ModelIntegral& modelIntegral() const noexcept
    { return ModelIntegral_; }

    /// Evaluates the integral of the model.
    void integrate();

private:

    /// @brief The fit model.
    /// @details It's shared with the `miFit`, which will change
    /// the free amplitudes of the model.
    const std::shared_ptr<const FitModel> FitModel_;

    /// The random number engine.
    RandomEngine RandomEngine_;

    /// The number of integration points (?)
    const unsigned NumberOfPoints_;
    /// Number of threads used for the integration.
    const unsigned NumberOfThreads_;
    /// (?)
    const unsigned BatchSize_;

    yap::DataSet             Data_;
    yap::DataPartitionVector DataPartitions_;
    yap::ModelIntegral       ModelIntegral_;

    /// Point generator for the integration.
    Generator PointGenerator_;
};

#endif
