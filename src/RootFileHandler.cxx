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

#include "TTreeHandler.h"

#include <Exceptions.h>

#include <TFile.h>

#include <iostream>
#include <string>
#include <memory>

RootFileHandler::RootFileHandler(const char* file_path,
                                 const char* file_name,
                                 const char* model_name) :
    Path_(static_cast<std::string>(file_path) + "/./"),
    FileName_(Path_ + static_cast<std::string>(file_name)),
    ModelName_(model_name),
    File_(TFile::Open(FileName_.data(), "READ")),
    MCMCTree_(std::make_unique<TTreeHandler>())
{
    std::cerr << " > Open input ROOT file: " << fileName() << " ... ";
    if (file()->IsZombie())
        throw yap::exceptions::Exception("Can't open: " + fileName(), "RootFileHandler::RootFileHandler");
    std::cerr << "done!" << std::endl;

    // Associate the mcmc TTree to the file tree.
    MCMCTree_->associate(file(), mcmc_tree_name(*this));
}

RootFileHandler::~RootFileHandler() = default;
