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
                                 const char* file_name) :
    Path_(static_cast<std::string>(file_path) + "/"),
    FileName_(Path_ + static_cast<std::string>(file_name)),
    File_(TFile::Open(FileName_.data(), "READ")),
    MCMCTree_(std::make_unique<TTreeHandler>()),
    ParameterTree_(std::make_unique<TTreeHandler>())
{
    std::cerr << " > Open input ROOT file: " << fileName() << " ... ";
    if (file()->IsZombie())
        throw yap::exceptions::Exception("Can't open: " + fileName(), "RootFileHandler::RootFileHandler");
    std::cerr << "done!" << std::endl;

    // Associate the TTrees to the file trees.
    MCMCTree_->associate(file(), mcmcTreeName());
    ParameterTree_->associate(file(), parameterTreeName());
}

RootFileHandler::~RootFileHandler() = default;
