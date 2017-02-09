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

#ifndef  FIT_INTEGRATOR_H
#define  FIT_INTEGRATOR_H

#include "fwd/FitModel.h"

#include <fwd/FourVector.h>

#include <ModelIntegral.h>

#include <functional>
#include <random>

/// Abstract base class for sample-integrators of the intensity over the model.
class FitIntegrator {
public: 

    /// @brief Constructor.
    /// @param fit_model The model according to which one wants to generate integration data.
    /// @param integration_points The number of points to generate.
    explicit FitIntegrator(std::shared_ptr<const FitModel> fit_model,
                           const unsigned integration_points) noexcept;

    /// _Default_ destructor.
    virtual ~FitIntegrator() = default;

    /// Returns the number of integration points.
    const unsigned numberOfPoints() const noexcept
    { return NumberOfPoints_; }

    /// Returns the number of threads used for the integration.
    const unsigned numberOfThreads() const noexcept
    { return NumberOfThreads_; }

    /// Returns the matrix of the model component integrals.
    const yap::ModelIntegral& modelIntegral() const noexcept
    { return ModelIntegral_; }

    /// Returns the fit model.
    const std::shared_ptr<const FitModel>& fitModel() const noexcept
    { return FitModel_; }

    /// Evaluates the integral of the model.
    virtual void integrate() = 0;

private:

    /// @brief The fit model.
    /// @details It's shared with the `miFit`, which will change
    /// the free amplitudes of the model.
    const std::shared_ptr<const FitModel> FitModel_;

    /// The number of integration points.
    const unsigned NumberOfPoints_;

    /// Number of threads used for the integration.
    const unsigned NumberOfThreads_;

protected:

    /// The integral of the model.
    yap::ModelIntegral ModelIntegral_;
};

#endif
