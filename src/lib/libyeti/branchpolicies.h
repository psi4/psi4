#ifndef BRANCHPOLICIES_H
#define BRANCHPOLICIES_H

#include "tuple.h"
#include "tensor.h"

namespace yeti {

class DoNothingMempool {
    public:
        void retrieve(
            TensorBlock* block
        );
};

class ResetMempool {
    public:
        void retrieve(
            TensorBlock* block
        );
};

class FlushOldBranchRenew {
    public:
        void renew(
            TensorBlock* block
        );
};

class ZeroBranchRenew {
    public:
        void renew(
            TensorBlock* block
        );
};


class RemoteAccumulateRenew :
    public ZeroBranchRenew
{
    public:
        void renew(
            TensorBlock* block
        );
};

class DoNothingBranchRenew {
    public:
        void renew(
            TensorBlock* block
        );
};

class RealignMemoryPoolBranchRetrieve {
    public:
        void retrieve(
            TensorBlock* block
        );
};

class NewBranchRetrieve {
    public:
        void retrieve(
            TensorBlock* block
        );
};

class ActionBranchRetrieve {
    public:
        void retrieve(
            TensorBlock* block
        );
};

class SortedBranchRetrieve
{
    public:
        void retrieve(
            TensorBlock* block
        );

        static void retrieve(
            TensorBlock* block,
            TensorBlock* unique_block
        );

        static void retrieve(
            TensorBlock* block,
            TensorBlock* unique_block,
            Permutation* perm
        );
};

class ConfigureElementComputerBranchRetrieve {
    public:
        void retrieve(
            TensorBlock* block
        );
};

class ConfigureElementComputerAndSortBranchRetrieve :
    ConfigureElementComputerBranchRetrieve,
    SortedBranchRetrieve
{
    public:
        void retrieve(
            TensorBlock* block
        );
};

class DoNothingBranchRetrieve {
    public:
        void retrieve(
            TensorBlock* block
        );
};

class RemoteBlockBranchRetrieve {
    public:
        void retrieve(
            TensorBlock* block
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
            TensorBlock* block
        );
};

class ThreadAccumulateBranchRelease {
    public:
        void release(
            TensorBlock* block
        );
};

class SetFinalizedBranchRelease {
    public:
        void release(
            TensorBlock* block
        );
};

class CacheBranchRelease {
    public:
        void release(
            TensorBlock* block
        );
};

class RemoteAccumulateBranchRelease {
    public:
        void release(
            TensorBlock* block
        );
};

class DoNothingPreflush {
    public:
        void preflush(
            TensorBlock* block
        );
};


class DoNothingBranchFlush :
    public DoNothingPreflush
{
    public:
        void flush(
            TensorBlock* block
        );
};

class ClearBranchFlush :
    public DoNothingPreflush
{
    public:
        void flush(
            TensorBlock* block
        );
};

class RemoteAccumulateFlush :
    public ClearBranchFlush
{
    public:
        void flush(
            TensorBlock* block
        );

        void preflush(
            TensorBlock* block
        );
};

class CommitBranchFlush :
    public DoNothingPreflush
{
    public:
        void flush(
            TensorBlock* block
        );
};

class SortedAccumulateBranchFlush :
    public DoNothingPreflush
{
    public:
        void flush(
            TensorBlock* block
        );
};

class ReuseDataControllers {

    private:
        void retrieve(TensorDataController* controller);

    public:
        void retrieve(
            TensorBlock* block
        );
};

class ResetDataControllers {
    public:
        void retrieve(
            TensorBlock* block
        );
};

class ReallocateDataControllers {
    public:
        static void retrieve(
            TensorBlock* block
        );
};

class DoNothingDataControllerRetrieve {
    public:
        void retrieve(
            TensorBlock* block
        );
};

class MemsetDataControllers {
    public:
        void retrieve(
            TensorBlock* block
        );
};

class SortDataControllers
{
    public:
        static void retrieve(
            TensorBlock* block,
            TensorBlock* unique_block,
            Permutation* perm
        );

        static void retrieve(
            TensorBlock* block
        );
};

class DoNothingDataControllerInit {
    public:
        void retrieve(
            TensorBlock* block
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

class AbortOnSync {
    public:
        void sync(
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

class RemoteAccumulateSync {
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

class DoNothingOutOfCorePrefetch {
    public:
        void prefetch(TensorBlock* block);

};

class ParentBlockReadPrefetch {
    public:
        void prefetch(TensorBlock* block);
};

class RemoteBlockPrefetch {
    public:
        void prefetch(TensorBlock* block);

};

class DoNothingInCorePrefetch {
    public:
        void prefetch(
            TensorBlock* current_block,
            TensorBlock* prev_block
        );
};

class ResortInCorePrefetch
{

    public:
        void prefetch(
            TensorBlock* current_block,
            TensorBlock* prev_block
        );
};


}

#endif // BRANCHPOLICIES_H
