/**
 *
 *    @file  RootFileHandler.cxx
 *   @brief  
 *
 *    @date  02/04/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "RootFileHandler.h"

#include <Exceptions.h>

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <memory>

RootFileHandler::RootFileHandler(const char* file_name,
                                 const char* model_name) :
    ModelName_(model_name),
    RootFileName_(file_name),
    RootFile_(TFile::Open((RootFileName_).data(), "READ")),
    MCMCTree_(nullptr),
    ParameterTree_(nullptr)
{
    std::cerr << " > Open input ROOT file: " << RootFileName_ << " ... ";
    if (RootFile_->IsZombie())
        throw yap::exceptions::Exception("Can't open: " + RootFileName_, "ModelHandler::ModelHandler");
    std::cerr << "done!" << std::endl;


    std::cerr << " > Load " + mcmcTreeName() + " tree ... ";
    TTree* tmp = nullptr;
    RootFile_->GetObject(mcmcTreeName().data(), tmp);
    if (!tmp)
        throw yap::exceptions::Exception("Can't retireve " + mcmcTreeName() + " tree",
                "RootFileHandler::RootFileHandler");

    // Initialize the MCMC tree.
    MCMCTree_.reset(tmp);
    std::cerr << "done!" << std::endl;

    std::cerr << " > Load " + parameterTreeName() + " tree ... ";
    tmp = nullptr;
    RootFile_->GetObject(parameterTreeName().data(), tmp);
    if (!tmp)
        throw yap::exceptions::Exception("Can't retireve " + parameterTreeName() + " tree",
                "RootFileHandler::RootFileHandler");

    // Initialize the parameter tree.
    ParameterTree_.reset(tmp);
    std::cerr << "done!" << std::endl;
}
