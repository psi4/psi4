#ifndef BRANCHPOLICIES_H
#define BRANCHPOLICIES_H

#include "tuple.h"
#include "tensor.h"

namespace yeti {

class DoNothingMempool {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class ResetMempool {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class RealignMemoryPoolBranchRetrieve {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class NewBranchRetrieve {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class ActionBranchRetrieve {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class SortedBranchRetrieve
{
    protected:
        void retrieve(
            TensorBlock* block,
            TensorBlock* unique_block
        );
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class ConfigureElementComputerBranchRetrieve {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class ResortAndConfigureElementComputerBranchRetrieve :
    public ConfigureElementComputerBranchRetrieve,
    public SortedBranchRetrieve
{
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class DoNothingBranchRetrieve {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class ValidBranchController {
    public:
        void validate(
            TensorBlock* block
        );
};

class AbortReadBranchValidation {
    public:
        void validate(
            TensorBlock* block
        );
};

class AbortWriteBranchValidation {
    public:
        void validate(
            TensorBlock* block
        );
};

class AbortAccumulateBranchValidation {
    public:
        void validate(
            TensorBlock* block
        );
};

class AbortVerbatimBranchValidation {
    public:
        void validate(
            TensorBlock* block
        );
};


class DoNothingBranchRelease {
    public:
        void release(
            TensorBranch* branch
        );
};

class DoNothingBranchFlush {
    public:
        void flush(
            TensorBranch* branch
        );
};

class ClearBranchFlush {
    public:
        void flush(
            TensorBranch* branch
        );
};

class CommitBranchFlush {
    public:
        void flush(
            TensorBranch* branch
        );
};

class SortedAccumulateBranchFlush {
    public:
        void flush(
            TensorBranch* branch
        );
};

class ReuseDataControllers {

    private:
        void retrieve(TensorDataController* controller);

    public:
        void retrieve(
            TensorBranch* branch
        );
};

class ResetDataControllers {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class ReallocateDataControllers {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class DoNothingDataControllerRetrieve {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class MemsetDataControllers {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class SortDataControllers
{
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class DoNothingDataControllerInit {
    public:
        void retrieve(
            TensorBranch* branch
        );
};

class AbortOnObsolete {
    public:
        void obsolete(
            TensorBlock* block
        );
};

class ClearMetaDataOnObsolete {
    public:
        void obsolete(
            TensorBlock* block
        );
};

class DoNothingSync {
    public:
        void sync(
            TensorBlock* block
        );
};

class FlushOnSync {
    public:
        void sync(
            TensorBlock* block
        );
};

class NoUpdate {
    public:
        void update(
            TensorBlock* block
        );
};

class UpdateMaxLog {
    public:
        void update(
            TensorBlock* block
        );
};

}

#endif // BRANCHPOLICIES_H
