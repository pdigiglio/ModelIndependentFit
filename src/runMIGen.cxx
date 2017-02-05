// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <HelicityFormalism.h>
#include <Model.h>
#include <ZemachFormalism.h>

#include <logging.h>
#include <make_unique.h>

#include "model-independent/d3pi_one_resonance.h"
#include "miFit.h"

#include <chrono>
#include <ratio>

using namespace std;
using namespace yap;

int main()
{
    plainLogs(el::Level::Info);

    const vector<miFit*> test_models = {
        new miFit("D3PI", d3pi_one_resonance(make_model<ZemachFormalism>()))
    };

    for (auto* m : test_models) {

        // open log file
        BCLog::OpenLog("output/" + m->GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

        // set precision
        m->SetPrecision(BCEngineMCMC::kMedium);
        m->SetNChains(4);
        m->SetMinimumEfficiency(0.15);
        m->SetMaximumEfficiency(0.35);

        m->SetNIterationsRun(static_cast<int>(1e6 / m->GetNChains()));

        m->WriteMarkovChain("output/" + m->GetSafeName() + "_mcmc.root", "RECREATE", true, false);

//        // start timing:
//        auto start = chrono::steady_clock::now();

        // run MCMC, marginalizing posterior
        m->MarginalizeAll(BCIntegrate::kMargMetropolis);

//        // end timing
//        auto end = chrono::steady_clock::now();

        // m->PrintAllMarginalized(m->GetSafeName() + "_plots.pdf");

//        // timing:
//        auto diff = end - start;
//        auto ms = chrono::duration<double, micro>(diff).count();
//        auto nevents = (m->GetNIterationsPreRun() + m->GetNIterationsRun()) * m->GetNChains();
//        BCLog::OutSummary(string("Seconds = ") + to_string(ms / 1.e6) + " for " + to_string(nevents) + " iterations, " + to_string(m->likelihoodCalls()) + " calls");
//        BCLog::OutSummary(to_string(ms / nevents) + " microsec / iteration");
//        BCLog::OutSummary(to_string(ms / m->likelihoodCalls()) + " microsec / call");

        // close log file
        BCLog::OutSummary("Exiting");
        BCLog::CloseLog();
    }

    return 0;
}
