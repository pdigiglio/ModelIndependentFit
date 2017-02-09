/**
 *
 *    @file  OnTheFlyIntegrator.h
 *   @brief  
 *
 *    @date  02/09/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  ON_THE_FLY_INTEGRATOR_H
#define  ON_THE_FLY_INTEGRATOR_H

#include "fwd/FitModel.h"

#include <fwd/FourVector.h>

#include "FitIntegrator.h"

#include <ModelIntegral.h>

#include <functional>
#include <random>

/// Sample-integrator of the intensity over the model.
class OnTheFlyIntegrator : public FitIntegrator {
public:
    /// Function for generating new points for integration.
    using Generator = std::function<std::vector<yap::FourVector<double> >()>;

    /// The type of the random-number engine.
    using RandomEngine = std::mt19937;

    /// @brief Constructor.
    /// @param fit_model The model according to which one wants to generate integration data.
    /// @param integration_points The number of points to generate.
    explicit OnTheFlyIntegrator(std::shared_ptr<const FitModel> fit_model,
                                const unsigned integration_points) noexcept;

    /// Returns the point generator (const version).
    const Generator& pointGenerator() const noexcept
    { return PointGenerator_; }

    /// @brief Returns the batch size.
    /// @todo Meaning of this?
    const unsigned batchSize() const noexcept
    { return BatchSize_; }

    /// Evaluates the integral of the model by generating data points on the fly.
    void integrate() final;

private:

    /// The random number engine.
    RandomEngine RandomEngine_;

    /// @brief The batch size.
    /// @todo Meaning of this?
    const unsigned BatchSize_;

    /// Point generator for the integration.
    Generator PointGenerator_;

    /// Access the point generator.
    Generator& pointGenerator() noexcept
    { return PointGenerator_; }
};

#endif
