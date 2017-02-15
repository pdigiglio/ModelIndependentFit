/**
 *
 *    @file  ModelIndependentFitModel.cxx
 *   @brief  
 *
 *    @date  02/14/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "ModelIndependentFitModel.h"

#include "Fit.h"
#include "FitResultFileIterator.h"
#include "FittedFreeAmplitude.h"
#include "MassBin.h"

#include "model/model_independent_d3pi.h"

#include <Attributes.h>
#include <BreitWigner.h>
#include <DecayTree.h> // for to_string
#include <DecayingParticle.h>
#include <Exceptions.h>
#include <FinalStateParticle.h>
#include <FreeAmplitude.h>
#include <MassAxes.h>
#include <MassRange.h>
#include <Model.h>
#include <PDL.h>
#include <Parameter.h>
#include <ZemachFormalism.h>

#include <memory>
#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <vector>


