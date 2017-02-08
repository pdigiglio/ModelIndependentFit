/**
 *
 *    @file  TTreeHandler.h
 *   @brief  
 *
 *    @date  02/05/17
 *  @author  Paolo Di Giglio (github.com/pdigiglio),
 *           <p.digiglio91@gmail.com>
 *
 */

#ifndef  TTREE_HANDLER_H
#define  TTREE_HANDLER_H

#include <Exceptions.h>

#include <TFile.h>
#include <TTree.h>

#include <cassert>
#include <memory>
#include <string>

/// Wrapper for a `TTree` with error-checking facilities.
class TTreeHandler {
public:

    /// _Default_ constructor.
    explicit TTreeHandler() :
        Tree_(std::make_unique<TTree>()) 
    {}

    /// _Default_ destructor.
    ~TTreeHandler();

    /// @brief Constructor.
    /// @details Directly associates the internal tree to the file tree.
    explicit TTreeHandler(const std::unique_ptr<TFile>& root_file, const std::string& tree_name);

    /// @brief Associate the tree with the data in the file.
    /// @param root_file The input ROOT file where to look for the tree.
    /// @param tree_name The name of the TTree to look for.
    void associate(const std::unique_ptr<TFile>& root_file, const std::string& tree_name);

    /// @brief Associates a ROOT branch (in the held TTree) with a variable.
    /// @param branch_name The name of the branch to read.
    /// @tparam T The data type contained in the branch.
    /// @return A pointer to the data in the branch.
    template <typename T>
    std::unique_ptr<T> setBranchAddress(const std::string branch_name) const {
        assert(Tree_ != nullptr);

        // Allocate a temporary variable to associate to the branch.
        auto tmp = new T();
        assert(tmp != nullptr);

        // Try to associate the variable to the branch
        const auto error_value = Tree_->SetBranchAddress(branch_name.data(), tmp);
        if (error_value != TTree::ESetBranchAddressStatus::kMatch)
            throw yap::exceptions::Exception("Can't set the branch address for " + branch_name,
                                             "TTreeHandler::setBranchAddress");

        assert(tmp != nullptr);
        return std::unique_ptr<T>(tmp);
    }

    /// @brief Associates a ROOT branch (in the held TTree) with a variable.
    /// @param branch_name The name of the branch to read.
    /// @param v           The variable to associate the branch to.
    /// @tparam T The data type contained in the branch.
    template <typename T>
    void setBranchAddress(const std::string branch_name, T* v) const {
        assert(v != nullptr);
        assert(Tree_ != nullptr);
        assert(branch_name.data() != nullptr);

        const auto error_value = Tree_->SetBranchAddress(branch_name.data(), v);
        if (error_value != TTree::ESetBranchAddressStatus::kMatch)
            throw yap::exceptions::Exception("Can't set the branch address for " + branch_name,
                                             "TTreeHandler::setBranchAddress");
    }

    /// Returns the entries in the tree.
    const Long64_t entries() const noexcept
    { return Tree_->GetEntries(); }

    /// @brief Reads the branches of _entry_.
    /// @param entry The entry to read.
    /// @return The number of bytes read.
    /// @throw yap::exceptions::Exception If the entry doesn't exist or an I/O error occurs.
    const Int_t getEntry(Long64_t entry) const; 

private:
    /// The ROOT TTree pointer.
    std::unique_ptr<TTree> Tree_;

    /// Return the internal tree (useful to call ROOT member functions).
    const std::unique_ptr<TTree>& tree() const noexcept
    { return Tree_; }
};

#endif
