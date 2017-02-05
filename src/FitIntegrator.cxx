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

#include <DataPartition.h>
#include <DataSet.h>
#include <FourVector.h>
#include <ImportanceSampler.h>
#include <PHSP.h>

#include <memory>
#include <random>
#include <type_traits>

FitIntegrator::FitIntegrator(std::shared_ptr<const FitModel> fit_model,
                             const unsigned integration_points) noexcept :
    FitModel_(fit_model),
    RandomEngine_(0),
    NumberOfPoints_(integration_points),
    NumberOfThreads_(available_threads()),
    BatchSize_(NumberOfPoints_ / NumberOfThreads_),
    Data_(FitModel_->model()->createDataSet()),
    DataPartitions_(NumberOfThreads_, &Data_),
    ModelIntegral_(*fit_model->model())
{
//    /// This doesn't work if the integrator is instanciated as constant
//    const_cast<FitIntegrator*>(this)->pointGenerator() = std::bind(yap::phsp<RandomEngine>, std::cref(*FitModel_->model()), FitModel::Dmass(), FitModel_->massRanges(), RandomEngine_, std::numeric_limits<unsigned>::max());
    PointGenerator_ = [&](){
        return yap::phsp(*FitModel_->model(), FitModel::Dmass(), FitModel_->axes(), FitModel_->massRanges(), RandomEngine_, std::numeric_limits<unsigned>::max());
    };

}

void FitIntegrator::integrate() {
//    assert(PointGenerator_);
    yap::ImportanceSampler::calculate(ModelIntegral_, pointGenerator(), numberOfPoints(), batchSize(), numberOfThreads());
}
