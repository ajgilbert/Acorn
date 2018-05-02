#include "Acorn/Analysis/interface/BranchHandlerBase.h"

namespace ac {
BranchHandlerBase::BranchHandlerBase()
    : branch_ptr_(nullptr), current_(-1), no_overwrite_(false) {}

BranchHandlerBase::~BranchHandlerBase() {}
}