/**
 *
 *    @file  parameter_handling.cxx
 *   @brief  
 *
 *    @date  03/01/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "parameter_handling.h"

#include "ModelIndependentFitModel.h"
#include "FreedWave.h"

#include <FreeAmplitude.h>

#include <algorithm>
#include <cassert>
#include <memory>
#include <vector>

std::vector<std::vector<std::complex<double>>> fit_to_yap_parameters(const std::vector<double>& p,
                                                                     const std::shared_ptr<const ModelIndependentFitModel>& m) noexcept {
    // Check that the number of parameters is correct.
    assert(p.size() == std::accumulate(std::begin(m->freedWaves()), std::end(m->freedWaves()), size_t(0), 
                       [](const size_t s, const auto& fw) { return s + 2 * fw->massBins().size(); }));
    
    // The vector of YAP-parameters vector.
    std::vector<std::vector<std::complex<double>>> yp(m->freedWaves().size(),
                                                      std::vector<std::complex<double>>({}));

    size_t offset = 0;
    for (size_t i = 0; i < m->freedWaves().size(); ++ i) {
        // Cache the freed wave.
        const auto& fw = m->freedWaves().at(i);
        // Reserve the memory for the YAP parameters.
        yp.at(i).reserve(fw->massBins().size());

        // Convert the parameters.
        double cumulative_phase = 0.;
        for (size_t j = offset; j < offset + fw->massBins().size(); ++ j) {
            cumulative_phase += p.at(2 * j + 1);
            yp[i].emplace_back(std::polar<double>(p[2 * j], yap::rad(cumulative_phase)));
        }

        // Update the offset.
        assert(fw->massBins().size() == yp[i].size());
        offset += yp[i].size();
    }

    assert(m->freeAmplitudes().size() == offset);
    return yp;
}

std::vector<double> yap_to_fit_parameters(const std::vector<std::complex<double>>& p) noexcept {
    std::vector<double> fp;
    fp.reserve(2 * p.size());

    // TODO don't hard-code this!!
    double cumulative_phase = yap::deg(std::arg(*std::begin(p)));
    for (const auto& yp : p) {
        fp.push_back(std::abs(yp));
        fp.push_back(yap::deg(std::arg(yp)) - cumulative_phase);

        cumulative_phase += fp.back();
    }

    return fp;
}

std::vector<double> yap_to_fit_parameters(const std::vector<std::vector<std::complex<double>>>& p) noexcept {
    std::vector<double> fp;
    fp.reserve(std::accumulate(std::begin(p), std::end(p), 0,
                               [](const size_t s, const auto& pars) { return s + 2 * pars.size(); }));

    std::for_each(std::begin(p), std::end(p), [&](const auto& yps) {
        auto tmp_fp = yap_to_fit_parameters(yps);
        std::copy(std::make_move_iterator(std::begin(tmp_fp)), std::make_move_iterator(std::end(tmp_fp)),
                  std::back_inserter(fp));
        });

    return fp;
}

std::vector<double> yap_to_fit_parameters(const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& fas,
                                          const std::shared_ptr<const ModelIndependentFitModel>& m) noexcept {
    // Check that the number of free amplitudes is correct.
    assert(fas.size() == std::accumulate(std::begin(m->freedWaves()), std::end(m->freedWaves()), size_t(0), 
                         [](const size_t s, const auto& fw) { return s + fw->massBins().size(); }));
    std::vector<std::vector<std::complex<double>>> p(m->freedWaves().size(),
                                                     std::vector<std::complex<double>>({}));

    size_t offset = 0;
    for (size_t i = 0; i < m->freedWaves().size(); ++ i) {
        // Cache the freed wave.
        const auto& fw = m->freedWaves()[i];
        p.at(i).reserve(fw->massBins().size());

        std::transform(std::begin(fas) + offset, std::begin(fas) + offset + fw->massBins().size(),
                       std::back_inserter(p[i]), [] (const auto& fa) { return fa->value(); });

        assert(p[i].size() == fw->massBins().size());

        offset += p[i].size();
    }

    return yap_to_fit_parameters(p);
}

bool check_fas_against_parameters(const std::vector<std::shared_ptr<const yap::FreeAmplitude>>& fas,
                                  const std::vector<std::vector<std::complex<double>>>& p) noexcept {
    // Check that the sizes match.
    assert(fas.size() == std::accumulate(std::begin(p), std::end(p), size_t(0),
                                         [] (const size_t s, const auto& pars) { return s + pars.size(); }));

    size_t offset = 0;
    for (const auto& pars : p) {
        const auto e = std::equal(std::begin(pars), std::end(pars), std::begin(fas) + offset,
                                  [](const auto& par, const auto& fa) { return par == fa->value(); });
        if (!e)
            return e;
        offset += pars.size();
    }
    
    assert(offset == fas.size());
    return true;
}
