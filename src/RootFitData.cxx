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
#include "parallelism.h"

#include <FourMomenta.h>
#include <FourVector.h>
#include <DataPartition.h>
#include <DataSet.h>
#include <logging.h>

#include <iterator>
#include <memory>
#include <vector>

// TODO make lambda
inline void set_address(const yap::MassAxes::value_type& a, std::vector<double>& m2, TTree& t_mcmc)
{
    m2.push_back(0);
    t_mcmc.SetBranchAddress(indices_string(*a, "m2_", "").data(), &m2.back());
}

// TODO return a yap::DataSet
size_t load_data(yap::DataSet& data, const yap::Model& M, const yap::MassAxes& A,
                 double initial_mass, TTree& t_mcmc, int N /* = -1*/, int lag  /*= -1*/ )
{
    if (A.empty())
        throw yap::exceptions::Exception("mass axes empty", "load_data");

    // set branch addresses
    std::vector<double> m2;
    m2.reserve(A.size());

    std::for_each(A.begin(), A.end(), std::bind(set_address, std::placeholders::_1, std::ref(m2), std::ref(t_mcmc)));
//    std::for_each(std::begin(A), std::end(A),
//                  [&](const yap::MassAxes::value_type& a) { return set_address(a, m2, t_mcmc); });

    if (m2.size() != A.size())
        throw yap::exceptions::Exception("Not all mass axes loaded from TTree", "load_data");

    //
    // load data
    //
    int Phase = -1;
    t_mcmc.SetBranchAddress("Phase", &Phase);
    unsigned Iteration;
    t_mcmc.SetBranchAddress("Iteration", &Iteration);
    unsigned Chain;
    t_mcmc.SetBranchAddress("Chain", &Chain);

    unsigned long long n_entries = t_mcmc.GetEntries();


    if (N < 0)
        // attempt to load all data
        N = n_entries;

    if (lag < 0)
        // calculate lag
        lag = n_entries / N;
    lag = std::max(lag, 1);

    int n_attempted = 0;
    size_t old_size = data.size();

    // Read input data and make sure that they're in the pase space.
    for (unsigned long long n = 0; n < n_entries and n_attempted < N; ++n) {
        t_mcmc.GetEntry(n);

        if (Phase <= 0)
            continue;

        if (Iteration % lag != 0)
            continue;

        // if (fabs(m2[0] - 1.35 * 1.35) > 0.1 or m2[1] > 1.55 or m2[1] < 0.58)
        //     continue;

        ++n_attempted;

        auto P = calculate_four_momenta(initial_mass, M, A, m2);
        if (P.empty()) {
            std::cout << "Point is out of phase space!" << std::endl;
            continue;
        }
        data.push_back(P);
    }

    if (data.empty())
        LOG(INFO) << "No data loaded.";
    else {
        LOG(INFO) << "Loaded " << data.size() - old_size << " data points (" << ((data.size() - old_size) * data[0].bytes() * 1.e-6) << " MB)"
                << " from a tree of size " << n_entries << ", with a lag of " << lag;
        if (old_size != 0)
            LOG(INFO) << "Total data size now " << data.size() << " points (" << (data.bytes() * 1.e-6) << " MB)";
    }

    if (int(data.size() - old_size) < N)
        LOG(WARNING) << "could not load as many data points as requested. Reduce the lag (or set it to -1 to automatically determine the lag).";

    return data.size() - old_size;
}

RootFitData::RootFitData(std::unique_ptr<RootFileHandler> rfh,
                         std::shared_ptr<const FitModel> fit_model) :
    RootFile_(std::move(rfh)),
    FitModel_(fit_model),
    Data_(FitModel_->model()->createDataSet()),
    Partitions_(available_threads(), &Data_)
{
    load_data(Data_, *FitModel_->model(), FitModel_->axes(), FitModel::Dmass(), *RootFile_->mcmcTree(), 1e4, -1);
    Partitions_ = yap::DataPartitionBlock::create(Data_, available_threads());
}
