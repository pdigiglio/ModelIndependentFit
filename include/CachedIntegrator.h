/**
 *
 *    @file  CachedIntegrator.h
 *   @brief  
 *
 *    @date  02/09/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  CACHED_INTEGRATOR_H
#define  CACHED_INTEGRATOR_H

#include "fwd/FitModel.h"

#include <fwd/DataPartition.h>
#include <fwd/DataSet.h>

#include "FitIntegrator.h"

/// @brief Sample-integrator of the intensity over the model.
/// @details Its integration data set is evaluated once at construction time and cached.
class CachedIntegrator : public FitIntegrator {
public:
    /// @brief Constructor.
    /// @param fit_model The model according to which one wants to generate integration data.
    /// @param integration_points The number of points to generate.
    explicit CachedIntegrator(std::shared_ptr<const FitModel> fit_model,
                              const unsigned integration_points) noexcept;

    /// _Default_ destructor.
    ~CachedIntegrator() final;

    /// Evaluates the integral of the model using the cached data points.
    void integrate() final;

private:
    /// Monte-Carlo generated integration data.
    const std::unique_ptr<yap::DataSet> Data_;
    yap::DataPartitionVector Partitions_;
};

#endif
