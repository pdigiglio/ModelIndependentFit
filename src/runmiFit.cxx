// ***************************************************************
// This file was created using the bat-project script
// for project miFit.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include "Fit.h"
#include "FitIntegrator.h"
#include "FitModel.h"
#include "RootFileHandler.h"
#include "RootFitData.h"

#include <cassert>
#include <iostream>
#include <memory>

int main(int argc, const char *argv[])
{
    if (argc != 4) {
        std::cerr << std::endl;
        std::cerr << "Usage:" << std::endl;
        std::cerr << argv[0] << " <file_path> <file_name> <model_name>" << std::endl;
        std::cerr << std::endl;
        std::abort();
    }

    const auto fit = create_fit(argv[1], argv[2], argv[3]);

    // set nicer style for drawing than the ROOT default
    BCAux::SetStyle();

    // open log file
    BCLog::OpenLog(fit->logFileName(), BCLog::detail, BCLog::detail);

    // set precision
    fit->SetPrecision(BCEngineMCMC::kMedium);

    BCLog::OutSummary("Test model created");
    
    BCLog::OutSummary("Generating integration points on the fly");
    
    //////////////////////////////
    // perform your analysis here

    // Normalize the posterior by integrating it over the full parameter space
    // fit->Normalize();

    // Write Markov Chain to a ROOT file as a TTree
    // fit->WriteMarkovChain(fit->GetSafeName() + "_mcmc.root", "RECREATE");

    fit->SetPrecision(BCEngineMCMC::kMedium);
    fit->SetNIterationsPreRunMax(1e2);
    fit->SetNChains(4);
//    fit->SetInitialPositionAttemptLimit(1e8);
    fit->SetNIterationsRun(static_cast<int>(50e3 / fit->GetNChains()));

    // run MCMC, marginalizing posterior
//    fit->MarginalizeAll(BCIntegrate::kMargMetropolis);


    fit->MCMCUserInitialize();

    const auto guess = guess_parameters(*fit);
    fit->FindMode(guess);//fit->GetBestFitParameters());

    write_fit_result_to_file(*fit);

    // draw all marginalized distributions into a PDF file
//    fit->PrintAllMarginalized(fit->pdfPlotsOutputName());

    // print summary plots
    // fit->PrintParameterPlot(fit->GetSafeName() + "_parameters.pdf");
    // fit->PrintCorrelationPlot(fit->GetSafeName() + "_correlation.pdf");
    // fit->PrintCorrelationMatrix(fit->GetSafeName() + "_correlationMatrix.pdf");
    // fit->PrintKnowledgeUpdatePlots(fit->GetSafeName() + "_update.pdf");

    // print results of the analysis into a text file
    fit->PrintSummary();

    // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
