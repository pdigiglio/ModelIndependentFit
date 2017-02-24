/**
 *
 *    @file  RealTimeParameterPlot.cxx
 *   @brief  
 *
 *    @date  02/21/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "RealTimeParameterPlot.h"

#include "Fit.h"

#include <FreeAmplitude.h>
#include <MathUtilities.h>

#include <TApplication.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TSystem.h>

#include <algorithm>
#include <cassert>
#include <iterator>
#include <memory>
#include <vector>

TROOT gmroot("gmroot","GEM monitor ROOT interface");
TApplication theApp("App", 0, 0);

std::unique_ptr<TGraph> amplitude_graph(const std::vector<double>& bin_low_edges,
                                        const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& fas) noexcept
{
    // Check that there are as many parameters as bins.
    assert(fas.size() == bin_low_edges.size());

    // Vector of free-amplitude amplitudes.
    std::vector<double> amplitudes;
    amplitudes.reserve(fas.size());
    std::transform(std::begin(fas), std::end(fas), std::back_inserter(amplitudes),
                   [](const auto& fa) { return std::abs(fa->value()); });

    return std::make_unique<TGraph>(amplitudes.size(), bin_low_edges.data(), amplitudes.data());
}

std::unique_ptr<TGraph> phase_graph(const std::vector<double>& bin_low_edges,
                                    const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& fas) noexcept
{
    // Check that there are as many parameters as bins.
    assert(fas.size() == bin_low_edges.size());

    // Vector of free-amplitude phases.
    std::vector<double> phases;
    phases.reserve(fas.size());

    const auto fp = yap_to_fit_parameters(fas);
    const auto phase_shift = yap::deg(std::arg(fas[0]->value())) - fp[1];
    std::transform(std::begin(fas), std::end(fas), std::back_inserter(phases),
                   [&](const auto& fa) { return yap::deg(std::arg(fa->value())) - phase_shift; });

    return std::make_unique<TGraph>(phases.size(), bin_low_edges.data(), phases.data());
}

// Helper function to create a TGraph whose y-values are all 0.
std::unique_ptr<TGraph> zero_graph(const std::vector<double>& bin_low_edges) noexcept
{
    const std::vector<double> zero_vector(bin_low_edges.size(), 0);
    return std::make_unique<TGraph>(zero_vector.size(), bin_low_edges.data(), zero_vector.data());
}

RealTimeParameterPlot::RealTimeParameterPlot(const std::vector<double>& bin_low_edges,
                                             const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& init_fas) noexcept :
    AmplitudeCanvas_(std::make_unique<TCanvas>("faa", "Free-amplitude amplitude")),
    FitAmplitude_(zero_graph(bin_low_edges)),
    GuessAmplitude_(amplitude_graph(bin_low_edges, init_fas)),
    PhaseCanvas_(std::make_unique<TCanvas>("fap", "Free-amplitude phase")),
    FitPhase_(zero_graph(bin_low_edges)),
    GuessPhase_(phase_graph(bin_low_edges, init_fas))
{
    // check that #bins == #parameters
    assert(bin_low_edges.size() == init_fas.size());

    // -- Associate plots to canvases ---
    assert(AmplitudeCanvas_);
    AmplitudeCanvas_->cd();
//    GuessAmplitude_->SetTitle("Free-amplitude amplitude"); // XXX disappears after update
    GuessAmplitude_->GetXaxis()->SetTitle("Mass [GeV/c^2]");
    GuessAmplitude_->GetYaxis()->SetTitle("Amplitude");
    GuessAmplitude_->SetLineColor(4);
    GuessAmplitude_->Draw("A L*");
    FitAmplitude_->SetLineColor(2);
    FitAmplitude_->SetMarkerColor(2);
    FitAmplitude_->Draw("SAME L*");

    assert(PhaseCanvas_);
    PhaseCanvas_->cd();
//    GuessPhase_->SetTitle("Free-amplitude phase"); // XXX disappears after update
    GuessPhase_->GetXaxis()->SetTitle("Mass [GeV/c^2]");
    GuessPhase_->GetYaxis()->SetTitle("Phase [deg]");
    GuessPhase_->SetLineColor(4);
    GuessPhase_->Draw("A L*");
    FitPhase_->SetLineColor(2);
    FitPhase_->SetMarkerColor(2);
    FitPhase_->Draw("SAME L*");
    // ----------------------------------

    AmplitudeCanvas_->Modified();
    AmplitudeCanvas_->Update();

    PhaseCanvas_->Modified();
    PhaseCanvas_->Update();

    gSystem->ProcessEvents();
}

RealTimeParameterPlot::~RealTimeParameterPlot() = default;

void RealTimeParameterPlot::setCurrentParameters(const std::vector<std::complex<double>>& p) const noexcept {
    const auto x_axis     = FitAmplitude_->GetX();
    const auto n_elements = FitAmplitude_->GetN();
    const std::vector<double> mass_partition(x_axis, x_axis + n_elements);
    assert(mass_partition.size() == p.size());

    for (size_t i = 0; i < p.size(); ++ i) {
        FitAmplitude_->SetPoint(i, mass_partition[i], std::abs(p[i]));
        FitPhase_->SetPoint(i, mass_partition[i], yap::deg(std::arg(p[i])));
    }

    AmplitudeCanvas_->Modified();
    AmplitudeCanvas_->Update();

    PhaseCanvas_->Modified();
    PhaseCanvas_->Update();

    gSystem->ProcessEvents();
}
