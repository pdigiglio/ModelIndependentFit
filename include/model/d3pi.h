/**
 *
 *    @file  d3pi.h
 *   @brief
 * @details  Utility functions to set up a \f$D^+ \to \pi^+\pi^-\pi^+\f$ model.
 *
 *    @date  02/14/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  D3PI_H
#define  D3PI_H

#include <fwd/Model.h>

#include <memory>

/// @brief Returns a YAP model for a \f$D^+ \to \pi^+\pi^-\pi^+\f$ decay.
/// @todo Add a vector of resonance names.
std::unique_ptr<yap::Model> d3pi();

#endif
