// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>


#include "DataGenerator.h"
#include "FitModel.h"
#include "parallelism.h"

#include <DataSet.h>
#include <logging.h>

#include <chrono>
#include <ratio>

using namespace std;
using namespace yap;

int main()
{
    plainLogs(el::Level::Info);

    const auto m(std::make_unique<DataGenerator>(std::make_shared<const FitModel>(d3pi())));

    // open log file
    BCLog::OpenLog("try_output/" + m->GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m->SetPrecision(BCEngineMCMC::kMedium);
    m->SetNChains(4);
    m->SetMinimumEfficiency(0.15);
    m->SetMaximumEfficiency(0.35);

    m->SetNIterationsRun(static_cast<int>(1e7 / m->GetNChains()));

    m->WriteMarkovChain("try_output/" + m->GetSafeName() + "_mcmc.root", "RECREATE", true, false);

    // run MCMC, marginalizing posterior
    m->MarginalizeAll(BCIntegrate::kMargMetropolis);

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
