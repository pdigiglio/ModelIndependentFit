/**
 *
 *    @file  d3pi_one_resonance.h
 *   @brief  
 *
 *    @date  09/30/16
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  __BAT__D3PI__RHO__ONLY__H
#define  __BAT__D3PI__RHO__ONLY__H

#include <fwd/Model.h>
#include "../miFit.h"

#include <memory>
#include <vector>

// use common radial size for all resonances
constexpr double radial_size = 3.; // [GeV^-1]

///// Creates a decay \f$D^+ \to \pi^+\pi^-\pi^+\f$ with only one resonance.
//std::unique_ptr<yap::Model> d3pi_one_resonance(std::unique_ptr<yap::Model> M);

/// @brief Creates a model with a rho0 and a freed _S_ wave.
std::unique_ptr<yap::Model> d3pi_binned();

miFit d3pi_binned_fit(const std::string& name);

std::vector<double> make_partition(double min, const double max, const unsigned int bins);
const std::vector<double> partition_mass_axis(double low_range, double up_range);

/// @brief Returns the available concurrency on the machine.
const unsigned available_threads() noexcept;

#endif
