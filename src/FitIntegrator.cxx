/**
 *
 *    @file  FitIntegrator.cxx
 *   @brief  
 *
 *    @date  02/03/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "FitIntegrator.h"

#include "FitModel.h"
#include "parallelism.h"

#include <Model.h>
#include <ModelIntegral.h>

#include <cassert>
#include <memory>

FitIntegrator::FitIntegrator(std::shared_ptr<const FitModel> fit_model,
                             const unsigned integration_points) noexcept :
    FitModel_(fit_model),
    NumberOfPoints_(integration_points),
    NumberOfThreads_(available_threads()),
    ModelIntegral_(*FitModel_->model())
{
    assert(FitModel_);
}
