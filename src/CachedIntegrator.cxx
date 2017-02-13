/**
 *
 *    @file  CachedIntegrator.cxx
 *   @brief  
 *
 *    @date  02/09/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "CachedIntegrator.h"

#include "FitModel.h"

#include <Model.h>
#include <ModelIntegral.h>

#include <DataPartition.h>
#include <DataSet.h>
#include <FourVector.h>
#include <ImportanceSampler.h>
#include <PHSP.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>
#include <random>
#include <type_traits>


yap::DataSet generate_data(const std::shared_ptr<const FitModel>& fm,
                           const unsigned integration_points) {
    std::mt19937 g(0);
    const auto max_attempts = std::numeric_limits<unsigned>::max();
    auto point_generator = [&]() {
        return yap::phsp(*fm->model(), FitModel::Dmass(), fm->massAxes(), fm->squaredMassRanges(), g, max_attempts);
    };

    // Create empty data.
    auto data_set = fm->model()->createDataSet();

    // Fill the data set with MC-generated data.
    std::generate_n(std::back_inserter(data_set), integration_points, point_generator);

    // Check if the data set is not empty.
    assert(!data_set.empty());
    return data_set;
}

CachedIntegrator::~CachedIntegrator() = default;

CachedIntegrator::CachedIntegrator(std::shared_ptr<const FitModel> fit_model,
                                   const unsigned integration_points) noexcept :
    FitIntegrator(fit_model, integration_points),
    Data_(std::make_unique<yap::DataSet>(generate_data(fitModel(), numberOfPoints()))),
    Partitions_(yap::DataPartitionBlock::create(*Data_, numberOfThreads()))
{
    assert(Data_ != nullptr);
    assert(!Data_->empty());
}

void CachedIntegrator::integrate() {
    yap::ImportanceSampler::calculate(ModelIntegral_, Partitions_);
}
