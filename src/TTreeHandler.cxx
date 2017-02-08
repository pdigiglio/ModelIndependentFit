/**
 *
 *    @file  TTreeHandler.cxx
 *   @brief  
 *
 *    @date  02/05/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#include "TTreeHandler.h"

#include <Exceptions.h>

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <memory>

TTreeHandler::TTreeHandler(const std::unique_ptr<TFile>& root_file,
                           const std::string& tree_name) :
    Tree_(std::make_unique<TTree>())
{
    associate(root_file, tree_name);
}

TTreeHandler::~TTreeHandler() = default;

void TTreeHandler::associate(const std::unique_ptr<TFile>& root_file,
                             const std::string& tree_name) {
    if (root_file->IsZombie())
        throw yap::exceptions::Exception("Can't open input file", "TTreeHandler::associate");

    std::cerr << " > Load " + tree_name + " tree ... ";

    // Create an empty branch.
    TTree* tmp = new TTree();

    // Associate the branch with the data in the file.
    root_file->GetObject(tree_name.c_str(), tmp);

    // Check that the pointer is valid.
    if (!tmp)
        throw yap::exceptions::Exception("Can't retireve " + tree_name + " tree",
                                         "TTreeHandler::associate");

    // Store the tree in the internal storage.
    Tree_.reset(tmp);

    std::cerr << "done!" << std::endl;
}

const Int_t TTreeHandler::getEntry(Long64_t entry) const {
    const auto error_value = Tree_->GetEntry(entry);
    if (error_value == 0)
        throw yap::exceptions::Exception("Entry doesn't exist", "TTreeHandler::getEntry");
    if (error_value == -1)
        throw yap::exceptions::Exception("I/O error", "TTreeHandler::getEntry");

    return error_value;
}
