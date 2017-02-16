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
#include "ModelIndependentFitModel.h"

#include "model/d3pi.h"
#include "model/model_independent_d3pi.h"

int main() {
    using namespace std;
    using namespace yap;

//    plainLogs(el::Level::Info);

//    // Generate a binned model with the bin amplitudes obtained from the fit results.
//    const auto m(std::make_unique<DataGenerator>(binned_d3pi_from_file("output/par_fit.txt"), "model_independent"));

    // Generate a binned model with the bin amplitudes guessed from the f0 mass shape.
    const auto m(std::make_unique<DataGenerator>(binned_d3pi("binned_f0"), "binned_f0"));

//    const auto m(std::make_unique<DataGenerator>(std::make_shared<FitModel>(d3pi(), "f0"), "f0"));

    // open log file
    BCLog::OpenLog("try_output/" + m->GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m->SetPrecision(BCEngineMCMC::kMedium);
    m->SetNChains(4);
    m->SetMinimumEfficiency(0.15);
    m->SetMaximumEfficiency(0.35);

    m->SetNIterationsRun(static_cast<int>(2e6 / m->GetNChains()));

    m->WriteMarkovChain("try_output/" + m->GetSafeName() + "_mcmc.root", "RECREATE", true, false);

    // run MCMC, marginalizing posterior
    m->MarginalizeAll(BCIntegrate::kMargMetropolis);

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
