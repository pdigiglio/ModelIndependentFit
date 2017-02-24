/**
 *
 *    @file  root_to_txt.cxx
 *   @brief  
 *
 *    @date  02/24/2017
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "FitModel.h"

#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>

#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>

int main (const int argc, const char *argv[]) {
    if (argc < 2 or argc > 3) {
        std::cerr << "Usage: " << std::endl;
        std::cerr << "./" << argv[0] << "<model_name> [<file_path>]" << std::endl;
        std::abort();
    }

    const std::string model_name(argv[1]);

    // Open the root file.
    std::unique_ptr<TFile> input_root_file(nullptr);
    {
        std::string file_path = ".";
        if (argc == 3)
            file_path = std::string(argv[2]);

        const auto root_file_name = file_path + "/" + model_name + "_mcmc.root";
        input_root_file.reset(TFile::Open(root_file_name.data(), "READ"));
        assert(!input_root_file->IsZombie());
    }

    std::unique_ptr<TTree> mcmc(nullptr);
    {
        TTree* tmp = new TTree();
        input_root_file->GetObject((model_name + "_mcmc").data(), tmp);
        assert(tmp);
        mcmc.reset(tmp);
    }
    
    constexpr size_t x_bins = 100;
    constexpr size_t y_bins = x_bins;

    constexpr auto x2_low  = std::pow((FitModel::PiPlusMass() + FitModel::PiMinusMass()), 2);
    constexpr auto x2_high = std::pow((FitModel::Dmass() - FitModel::PiPlusMass()), 2);
    constexpr auto y2_low  = x2_low;
    constexpr auto y2_high = x2_high;

    std::unique_ptr<TH2D> h(new TH2D("h", "h", x_bins, x2_low, x2_high, y_bins, y2_low, y2_high));
    mcmc->Draw("sqrt(m2_01):sqrt(m2_12)>>h", "", "colz");

    for (size_t i = 0; i < x_bins; ++ i) {
        const auto x = h->GetXaxis()->GetBinCenter(i);

        for(size_t j = 0; j < y_bins; ++ j) {
            const auto bc = h->GetBinContent(i, j);

            if (bc) {
                const auto y = h->GetYaxis()->GetBinCenter(j);
                std::cout << std::setw(10) << x << " "
                          << std::setw(10) << y << " "
                          << bc << std::endl;
            }
        }
    }

	return 0;
}
