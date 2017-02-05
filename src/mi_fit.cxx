/**
 *
 *    @file  mi_fit.cxx
 *   @brief  
 *
 *    @date  01/11/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "bat_fit.h"
//#include "hist.h"
#include "models/d3pi.h"
#include "models/d4pi.h"
#include "models/dkkpi.h"
#include "tools.h"

#include <Exceptions.h>
#include <HelicityFormalism.h>
#include <MassRange.h>
#include <PHSP.h>
#include <RelativisticBreitWigner.h>
#include <ZemachFormalism.h>
#include <logging.h>
#include <make_unique.h>

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <random>

#include <logging.h>

class ModelHandler {
public:
    explicit ModelHandler(std::string&& model_name) :
        ModelName_(std::move(model_name)),
        RootFile_(TFile::Open(("output/" + ModelName_ + "_mcm.root").data(), "READ")),
        MCMCTree_(RootFile_->GetObject((ModelName_ + "_mcmc").data(), "READ")),
        ParameterTree_(RootFile_->GetObject((ModelName_ + "_parameters").data(), "READ"))
    {
        if (RootFile_->IsZombie())
            throw yap::exceptions::Exception("can't open the file", "ModelHandler::ModelHandler");
        if (!MCMCTree_)
            throw yap::exceptions::Exception("can't retireve mcmc tree", "main");
        if (!ParameterTree_)
            throw yap::exceptions::Exception("can't retireve parameter tree", "main");
    }

    std::string& name() const noexcept
    { return ModelName_; }

    std::string& fitName() const noexcept
    { return name() + "_fit"; }

    TTree* const parameterTree() const noexcept
    { return ParameterTree_; }

private:
    const std::string ModelName_;
    TFile* RootFile_;
    TTree* MCMCTree_;
    TTree* ParameterTree_;
}

int main () {
    yap::plainLogs(el::Level::Info);

//    std::string model_name{"D3PI"};
//
//    auto file = TFile::Open(("output/" + model_name + "_mcmc.root").data(), "READ");
//    if (file->IsZombie())
//        throw yap::exceptions::Exception("can't open the file", "main");
//
//    TTree* t_mcmc{nullptr};
//    file->GetObject((model_name + "_mcmc").data(), t_mcmc);
//    if (!t_mcmc)
//        throw yap::exceptions::Exception("can't retireve mcmc tree", "main");
//
//    TTree* t_pars{nullptr};
//    file->GetObject((model_name + "_parameters").data(), t_pars);
//    if (!t_pars)
//        throw yap::exceptions::Exception("can't retireve mcmc tree", "main");

    const std::unique_ptr mh(new ModelHandler("D3PI"));

    // Fitting model
    const auto model{std::make_unique<bat_fit>(mh->fitName(), yap_model<yap::ZemachFormalism>(), find_mass_axes(*mh->parameterTree()))};
    constexpr double DPlus_mass{1.86961};

    LOG(INFO) << "Loading data...";
    load_data(model->fitData(), *model->model(), model->axes(), model->finalStateParticles());
    LOG(INFO) << "Data loaded.";

    auto m2r {yap::squared(mass_range(DPlus_mass, model->axes(), model->finalStateParticles()))};
    std::mt19937 g{0};
    model->integrationPointGenerator() = [=](void) { return yap::phsp<std::mt19937>(std::cref(*model->model()), DPlus_mass, model->axes(), m2r, g, std::numeric_limits<unsigned>::max()); };
    constexpr unsigned n_integrationPoints = 2e4;
    model->setNIntegrationPoints(n_integrationPoints, 1e4, 4);
    LOG(INFO) << "Generating integration points on the fly";

    // open log file
    BCLog::OpenLog("output/" + model->GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    model->SetPrecision(BCEngineMCMC::kMedium);
    model->SetNIterationsPreRunMax(1e6);
    model->SetNChains(4);

    model->SetNIterationsRun(static_cast<int>(50e3 / model->GetNChains()));

    // start timing:
    auto start = std::chrono::steady_clock::now();

    // run MCMC, marginalizing posterior
    model->MarginalizeAll(BCIntegrate::kMargMetropolis);

    // end timing
    auto end = std::chrono::steady_clock::now();

    model->FindMode(model->GetBestFitParameters());

    model->PrintSummary();
    model->PrintAllMarginalized("output/" + model->GetSafeName() + "_plots.pdf", 2, 2);

    // timing:
    auto diff = end - start;
    auto ms = std::chrono::duration<double, std::micro>(diff).count();
    auto nevents = (model->GetNIterationsPreRun() + model->GetNIterationsRun()) * model->GetNChains();
    BCLog::OutSummary(std::string("Seconds = ") + std::to_string(ms / 1.e6) + " for " + std::to_string(nevents) + " iterations, " + std::to_string(model->likelihoodCalls()) + " calls");
    BCLog::OutSummary(std::to_string(ms / nevents) + " microsec / iteration");
    BCLog::OutSummary(std::to_string(ms / model->likelihoodCalls()) + " microsec / call");

    // Close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
