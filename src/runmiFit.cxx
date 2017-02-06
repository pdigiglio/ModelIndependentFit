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

std::unique_ptr<Fit> create_fit(const char* file_path, const char* file_name, const char* model_name) {
    // Create a model to fit the data.
    // TODO the model name should be a member of the FitModel
    const std::shared_ptr<FitModel> fit_model(make_fit_model());
    assert(fit_model->freeAmplitudes().size() == 22);

    // Create the integrator.
    constexpr unsigned integration_points     = 2e4;
    auto integrator(std::make_unique<FitIntegrator>(std::static_pointer_cast<const FitModel>(fit_model), integration_points));

    auto root_file_handler(std::make_unique<RootFileHandler>(file_path, file_name, model_name));
    auto root_fit_data(std::make_unique<RootFitData>(std::move(root_file_handler), fit_model));

    // Create the BAT model for fitting the data.
    return std::make_unique<Fit>(fit_model, std::move(root_fit_data), std::move(integrator));
}

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

//    // Load data and partition it
//    constexpr double D_mass = 1.86961;
//    load_data(fit->fitData(), *fit->model(), fit->axes(), D_mass, *rfh.mcmcTree(), 1e4, -1);
//    fit->fitPartitions() = yap::DataPartitionBlock::create(fit->fitData(), available_threads());
    
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

//#include "model-independent/d3pi_one_resonance.h"
//#include "miFit.h"
//#include "FitModel.h"
//#include "Integrator.h"
//
////#include "model-independent/d3pi_one_resonance.h"
////#include "miFit.h"
////
////#include <AmplitudeBasis.h>
////#include <Attributes.h>
////#include <BreitWigner.h>
////#include <DataSet.h>
////#include <DecayChannel.h>
////#include <DecayTree.h>
////#include <DecayingParticle.h>
////#include <FinalStateParticle.h>
////#include <Flatte.h>
////#include <FourMomenta.h>
////#include <FourVector.h>
////#include <FreeAmplitude.h>
////#include <HelicityFormalism.h>
////#include <ImportanceSampler.h>
////#include <MassAxes.h>
////#include <MassRange.h>
////#include <MathUtilities.h>
////#include <Model.h>
////#include <ModelIntegral.h>
////#include <PDL.h>
////#include <PHSP.h>
////#include <Parameter.h>
////#include <ParticleCombination.h>
////#include <ParticleTable.h>
////#include <QuantumNumbers.h>
////#include <RelativisticBreitWigner.h>
////#include <SpinAmplitude.h>
////#include <SpinAmplitudeCache.h>
////#include <VariableStatus.h>
////#include <ZemachFormalism.h>
////#include <logging.h>
//
//#include <logging.h>
//#include <make_unique.h>
//
//#include <TCanvas.h>
//#include <TFile.h>
//#include <TTree.h>
//
//#include <cassert>
//#include <future>
//#include <iomanip>
//#include <type_traits>
//
//int main(int argc, const char *argv[])
//{
//    if (argc != 3) {
//        std::cerr << std::endl;
//        std::cerr << "Usage:" << std::endl;
//        std::cerr << argv[0] << " <file_path> <model_name>" << std::endl;
//        std::cerr << std::endl;
//        std::abort();
//    }
//
//    // set nicer style for drawing than the ROOT default
//    BCAux::SetStyle();
//
//    // -- Workflow example ---------------------------------------
//    // Create a model to fit the data.
//    const std::shared_ptr<FitModel> fit_model(make_fit_model());
//
//    // Create the integrator.
//    constexpr unsigned integration_points     = 2e4;
//    auto integrator(std::make_unique<Integrator>(fit_model, integration_points));
//
//    // TODO Create the container for the fit data.
//    // Load data
//    auto root_file_handler(std::make_unique<RootFileHandler>(argv[1], argv[2]));
//    auto root_fit_data(std::make_unique<RootFileData>(root_file_handler, fit_model));
//
//    // Create the BAT model for fitting the data.
//    auto miFit(std::make_unique<miFit>(fit_model, root_fit_data, integrator));
//    // -----------------------------------------------------------
//
//
//    // create new miFit object
//    miFit m(d3pi_binned_fit(rfh.modelName()));
//
//    // open log file
//    BCLog::OpenLog(m.logFileName(), BCLog::detail, BCLog::detail);
//
//    // set precision
//    m.SetPrecision(BCEngineMCMC::kMedium);
//
//    BCLog::OutSummary("Test model created");
//
//    // Load data and partition it
//    constexpr double D_mass = 1.86961;
//    load_data(m.fitData(), *m.model(), m.axes(), D_mass, *rfh.mcmcTree(), 1e4, -1);
//    m.fitPartitions() = yap::DataPartitionBlock::create(m.fitData(), available_threads());
//    
////    const auto m2r = yap::squared(mass_range(D_mass, m.axes(), m.model()->finalStateParticles()));
////    std::mt19937 g(0);
//
////    auto integration_point_generator =
//////    [&](miFit& m, const std::vector<yap::MassRange>& m2r) {
////    [&]() {
////        return yap::phsp(*m.model(), D_mass, m.axes(), m2r, generator, std::numeric_limits<unsigned>::max());
////    };
//
////    m.integrationPointGenerator() = std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), D_mass, m.axes(), m2r, g, std::numeric_limits<unsigned>::max());
////    m.integrationPointGenerator()  = integration_point_generator;
//    BCLog::OutSummary("Generating integration points on the fly");
//    
//    //////////////////////////////
//    // perform your analysis here
//
//    // Normalize the posterior by integrating it over the full parameter space
//    // m.Normalize();
//
//    // Write Markov Chain to a ROOT file as a TTree
//    // m.WriteMarkovChain(m.GetSafeName() + "_mcmc.root", "RECREATE");
//
//    m.SetPrecision(BCEngineMCMC::kMedium);
//    m.SetNIterationsPreRunMax(1e2);
//    m.SetNChains(4);
////    m.SetInitialPositionAttemptLimit(1e8);
//    m.SetNIterationsRun(static_cast<int>(50e3 / m.GetNChains()));
//
////    // Set flat priors for all the parameters
////    m.SetPriorConstantAll();
//
//    // run MCMC, marginalizing posterior
////    m.MarginalizeAll(BCIntegrate::kMargMetropolis);
//
//
//    m.MCMCUserInitialize();
//
//    const auto guess = guess_parameters(m);
//    m.FindMode(guess);//m.GetBestFitParameters());
//
//    write_fit_result_to_file(m);
//
//    // draw all marginalized distributions into a PDF file
////    m.PrintAllMarginalized(m.pdfPlotsOutputName());
//
//    // print summary plots
//    // m.PrintParameterPlot(m.GetSafeName() + "_parameters.pdf");
//    // m.PrintCorrelationPlot(m.GetSafeName() + "_correlation.pdf");
//    // m.PrintCorrelationMatrix(m.GetSafeName() + "_correlationMatrix.pdf");
//    // m.PrintKnowledgeUpdatePlots(m.GetSafeName() + "_update.pdf");
//
//    // print results of the analysis into a text file
//    m.PrintSummary();
//
//    // close log file
//    BCLog::OutSummary("Exiting");
//    BCLog::CloseLog();
//
//    return 0;
//}
