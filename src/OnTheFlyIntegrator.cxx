/**
 *
 *    @file  OnTheFlyIntegrator.cxx
 *   @brief  
 *
 *    @date  02/09/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "OnTheFlyIntegrator.h"

#include "FitModel.h"
#include "parallelism.h"

#include <FourVector.h>
#include <ImportanceSampler.h>
#include <PHSP.h>

#include <memory>
#include <random>
#include <type_traits>

OnTheFlyIntegrator::OnTheFlyIntegrator(std::shared_ptr<const FitModel> fit_model,
                             const unsigned integration_points) noexcept :
    FitIntegrator(fit_model, integration_points),
    RandomEngine_(0),
    BatchSize_(numberOfPoints() / numberOfThreads()),
    PointGenerator_([=]() { return yap::phsp(*fitModel()->model(), FitModel::Dmass(), fitModel()->massAxes(), fitModel()->squaredMassRanges(), RandomEngine_, std::numeric_limits<unsigned>::max()); })
{} 

void OnTheFlyIntegrator::integrate() {
    yap::ImportanceSampler::calculate(ModelIntegral_, pointGenerator(), numberOfPoints(), batchSize(), numberOfThreads());
}
