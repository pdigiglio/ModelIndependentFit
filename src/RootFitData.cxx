/**
 *
 *    @file  RootFitData.cxx
 *   @brief  
 *
 *    @date  02/04/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "RootFitData.h"

#include "FitModel.h"
#include "RootFileHandler.h"
#include "TTreeHandler.h"
#include "parallelism.h"

#include <FourMomenta.h>
#include <FourVector.h>
#include <DataPartition.h>
#include <DataSet.h>
#include <logging.h>

#include <TTree.h>

#include <cassert>
#include <iterator>
#include <memory>
#include <vector>

// TODO make lambda
inline void set_address(const yap::MassAxes::value_type& a,
                        std::vector<double>& m2,
                        const std::unique_ptr<TTreeHandler>& t_mcmc)
{
    m2.emplace_back(0);
    t_mcmc->setBranchAddress(indices_string(*a, "m2_", ""), &m2.back());
}

// Helper function to load data from the ROOT file to the data set.
yap::DataSet load_data( std::shared_ptr<const FitModel> fit_model,
                       const std::unique_ptr<TTreeHandler>& t_mcmc,
                       int N /* = -1*/,
                       int lag  /*= -1*/)
{
    // Make sure that the mass azes are not empty.
    const auto A = fit_model->massAxes();
    assert(!A.empty());

    // Set the branch addresses.
    std::vector<double> m2;
    m2.reserve(A.size());

    // Link the ROOT branches to m2's entries.
    std::for_each(std::begin(A), std::end(A),
                  [&](const yap::MassAxes::value_type& a) { return set_address(a, m2, t_mcmc); });

    // Check that all the mass axes have been loaded.
    assert(m2.size() == A.size());

    // Load data.
    auto Phase     = t_mcmc->setBranchAddress<int>("Phase");
    auto Iteration = t_mcmc->setBranchAddress<unsigned>("Iteration");
    auto Chain     = t_mcmc->setBranchAddress<unsigned>("Chain");

    const auto n_entries = t_mcmc->entries();

    // If N is negative, try to load all the data.
    N   = ((N < 0)   ? n_entries : N);
    // Evaluate the lag.
    lag = ((lag < 0) ? (n_entries / N) : lag);
    lag = std::max(lag, 1);

    // Create an empty data set.
    auto data_set = fit_model->model()->createDataSet();
    // Store the current data size.
    const size_t old_size = data_set.size();

    // Read input data and make sure that they're in the pase space.
    int n_attempted = 0;
    for (long long n = 0; n < n_entries and n_attempted < N; ++n) {
        t_mcmc->getEntry(n);

        if (*Phase <= 0)
            continue;

        if (*Iteration % lag != 0)
            continue;

        // if (fabs(m2[0] - 1.35 * 1.35) > 0.1 or m2[1] > 1.55 or m2[1] < 0.58)
        //     continue;

        ++n_attempted;

        auto P = calculate_four_momenta(FitModel::Dmass(), *fit_model->model(), A, m2);
        if (P.empty()) {
            std::cout << "Point is out of phase space!" << std::endl;

#ifndef NDEBUG
            // Print the vector
            std::cout << "( ";
            std::copy(std::begin(m2), std::end(m2), std::ostream_iterator<double>(std::cout, " "));
            std::cout << ")" << std::endl << std::endl;
#endif

            continue;
        }

        data_set.push_back(P);
    }

    // Check that some data were loaded.
    assert(!data_set.empty());

    LOG(INFO) << "Loaded " << data_set.size() - old_size << " data points ("
              << ((data_set.size() - old_size) * data_set[0].bytes() * 1.e-6) << " MB)"
              << " from a tree of size " << n_entries << ", with a lag of " << lag;

    if (int(data_set.size() - old_size) < N)
        LOG(WARNING) << "could not load as many data points as requested. Reduce the lag (or set it to -1 to automatically determine the lag).";

    return data_set;
}

RootFitData::RootFitData(std::unique_ptr<RootFileHandler> rfh,
                         std::shared_ptr<const FitModel> fit_model) :
    RootFile_(std::move(rfh)),
    FitModel_(fit_model),
    Data_(load_data(FitModel_, RootFile_->mcmcTree(), 1e4, -1)),
    Partitions_(yap::DataPartitionBlock::create(Data_, available_threads()))
{
    assert(RootFile_ != nullptr);
    assert(fit_model != nullptr);
    assert(Data_.size() == 1e4);
}

const std::string RootFitData::path() const noexcept
{ return RootFile_->path(); }
