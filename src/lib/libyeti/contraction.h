#ifndef yeti_contraction_h
#define yeti_contraction_h

#include "class.h"
#include "mapimpl.h"
#include "taskqueue.h"
#include "yetiobject.h"
#include "taskqueue.hpp"
#include "contraction.hpp"
#include "permutation.hpp"
#include "tensor.hpp"
#include "node.hpp"
#include "tensorblock.hpp"
#include "data.hpp"
#include "matrix.hpp"
#include "mallocimpl.h"
#include "thread.h"

#ifdef redefine_size_t
#define size_t custom_size_t
#endif

namespace yeti {

/**
    @class ContractionTask
    Task instance for accumulating product of two matrices into a product matrix
*/
class ContractionTask :
    public Task
{

    private:
        TensorBlock* lblock_;

        TensorBlock* rblock_;

        TensorBlock* pblock_;

        Contraction* cxn_;


    public:
        /**
            @param lmatrix
            @param rmatrix
            @param product_matrix
            @param cxn
        */
        ContractionTask(
            TensorBlock* lblock,
            TensorBlock* rblock,
            TensorBlock* pblock,
            Contraction* cxn
        );

        virtual ~ContractionTask();

        void print(std::ostream& os = std::cout) const;

        /**
            Polymorphic function called by TaskQueue
            @param threadnum The thread number
        */
        void run(uli threadnum);

        void out_of_core_prefetch();

        void in_core_prefetch(Task* prev_task);

        void finalize_task_subset();



};

class ContractionConfiguration
{
    private:
        MatrixConfigurationPtr lconfig_;

        MatrixConfigurationPtr rconfig_;

        MatrixConfigurationPtr pconfig_;

        MatrixIndexPtr lindex_;

        MatrixIndexPtr rindex_;

        MatrixIndexPtr pindex_;

        bool transpose_right_;

        bool transpose_left_;

        Tensor* product_tensor_;

        TensorBlock* tmp_block_;

    public:
        ContractionConfiguration(
            const MatrixConfigurationPtr& lconfig,
            const MatrixConfigurationPtr& rconfig,
            const MatrixConfigurationPtr& pconfig
        );

        ~ContractionConfiguration();

        void clone_product_tensor_for_threads(Tensor* tensor);
        
        Tensor* get_product_tensor();

        void configure_left_block(const uli* sizes, usi depth);

        void configure_right_block(const uli* sizes, usi depth);

        void configure_product_block(const uli* sizes, usi depth);

        void configure_left_block(MetaDataNode* node);

        void configure_right_block(MetaDataNode* node);

        void configure_product_block(MetaDataNode* node);

        void configure_left_block(Tensor* tensor);

        void configure_right_block(Tensor* tensor);

        void configure_product_block(Tensor* tensor);

        void set_tmp_accumulate_block(TensorBlock* pblock);

        TensorBlock* get_product_block(TensorBlock* pblock);

        TensorBlock* get_left_block(Tensor*, uli r, uli c) const;

        TileNode* get_left_node(MetaDataNode* node, uli r, uli c, uli& idx) const;

        TensorBlock* get_right_block(Tensor* tensor, uli r, uli c) const;

        TileNode* get_right_node(MetaDataNode* node, uli r, uli c, uli& idx) const;

        TensorBlock* get_product_block(Tensor* tensor, uli r, uli c) const;

        TileNode* get_product_node(MetaDataNode* node, uli r, uli c) const;

        uli ncxn_rows_left() const;

        uli ncxn_rows_right() const;

        uli ncxn_cols_left() const;

        uli ncxn_cols_right() const;

        uli nrows_left(const uli* sizes) const;

        uli ncols_left(const uli* sizes) const;

        uli nrows_right(const uli* sizes) const;

        uli ncols_right(const uli* sizes) const;

        uli get_left_index(uli r, uli c) const;

        uli get_right_index(uli r, uli c) const;

        uli get_product_index(uli r, uli c) const;

        void reset_contraction_depth(usi depth);

        MatrixConfiguration* get_left_config() const;

        MatrixConfiguration* get_right_config() const;

        MatrixConfiguration* get_product_config() const;

};

/**
    @class Contraction
    Class used for configuring a contraction.
*/
class Contraction :
    public TaskParent
{

    public:
        typedef enum {
            left_tensor = 0,
            right_tensor = 1,
            product_tensor = 2
        } tensor_position_t;

    private:    
        MatrixIndexPtr lindex_;
        
        MatrixIndexPtr rindex_;
        
        MatrixIndexPtr pindex_;
        
        PermutationGroupPtr lcxn_grp_;
        
        PermutationGroupPtr rcxn_grp_;
        
        bool transpose_right_;
        
        bool transpose_left_;
    
        /**
            The matrix configuration built by contraction for the left matrix
        */
        MatrixConfigurationPtr lconfig_;

        /**
            The matrix configuration built by contraction for the right matrix
        */
        MatrixConfigurationPtr rconfig_;

        /**
            The matrix configuration built by contraction for the product matrix
        */
        MatrixConfigurationPtr pconfig_;

        /**
            The final permutational symmetry of the target tensor
            accounting for all symmetries:
            Those imposed by the user, those automatically obtained
            by symmetry of left and right tensors, and special symmetries.
        */
        PermutationGroupPtr final_grp_;

        ContractionConfiguration** cxn_configs_;

        /**
            The set of permutations which permute only the row
            indices in the left matrix and only the col indices
            in the right matrix.  The final product tensor
            will automatically have this permutational symmetry.
            This also includes any "special" symmetries.  These
            are therefore the symmetries obtained by "default"
            without any post-contraction symmetrizations.
        */
        PermutationGroupPtr default_grp_;

        /**
            The set of permutations common to the contraction
            indices for the left matrix and right matrix
        */
        PermutationGroupPtr contraction_grp_;

        Tensor* ltensor_;

        Tensor* rtensor_;

        Tensor* ptensor_;

        Tensor* alpha_tensor_;

        Tensor* beta_tensor_;

        Tensor* gamma_tensor_;

        double scale_;

        tensor_position_t distribution_type_;

        tensor_position_t alpha_type_;

        tensor_position_t beta_type_;

        tensor_position_t gamma_type_;

        bool use_thread_replicated_product_;

        bool flush_product_on_finalize_;

        /**
            Depending on how the matrices are to be interpreted, the contraction
            engine manipulates the data appropriately.  For example, you could
            multiply and of the following cases (^T denotes transpose)
            L * R, L^T * R, L * R^T, L^T * R^T.
        */
        ContractionEngine* engine_;


        bool do_task(
            TensorBlock* lblock,
            TensorBlock* rblock,
            TensorBlock* pblock
        );

        void build_l_r_p();

        void build_l_p_r();

        void build_r_l_p();

        void build_r_p_l();

        void build_p_l_r();

        void build_p_r_l();

        void init_dot_product();
        
        void init_direct_product();
        
        void init_cxn();

        void init_cxn_grp();

    public:
        /**
            @param ltensor
            @param rtensor
            @param lindex. Defines what indices are cxn indices.
            @param rindex. Defines what indices are cxn indices.
            @param required_grp. The permutation group the target tensor
                                is required to have.
            @param lperm.  The left tensor can be sorted on the fly via
                            this permutation
            @param rperm.  The right tensor can be sorted on the fly via
                            this permutation
        */
        Contraction(
            double scale,
            Tensor* ltensor,
            Tensor* rtensor,
            Tensor* ptensor,
            const MatrixIndexPtr& lindex,
            const MatrixIndexPtr& rindex,
            Permutation* lperm,
            Permutation* rperm
        );

        virtual ~Contraction();

        void build();

        void finalize();

        /**
            @return The permutation group for the contraction indices
        */
        PermutationGroup* get_contraction_grp() const;

        /**
            @return The permutation group for the target indices that are automatically obtained
                    from the permutational symmetry of the tensors involved.
        */
        PermutationGroup* get_default_grp() const;

        /**
            @return The set of permutational symmetries the product tensor is required to have.
        */
        PermutationGroup* get_required_grp() const;

        double scale_factor() const;

        /**
            @return The final set of permutational symmetries the product tensor will have.
                    This is the union of the default grp, the required grp, and the special set.
        */
        PermutationGroup* get_final_grp() const;

        /**
            @return The set of permutations that must be applied to the product matrix
                    to make the final product tensor have the required symmetries.
        */
        PermutationSet* get_product_set() const;

        MatrixConfiguration* get_left_matrix_configuration() const;

        MatrixConfiguration* get_right_matrix_configuration() const;

        MatrixConfiguration* get_product_matrix_configuration() const;

        Tensor* get_left_tensor() const;

        Tensor* get_right_tensor() const;

        Tensor* get_product_tensor() const;

        bool is_thread_replicated_product() const;

        bool flush_product_on_finalize() const;
        
        /**
            Run all jobs on the contraction queue
        */
        static void run();

        /**
            Clear all jobs on the contraction queue
        */
        static void clear();

        void print(std::ostream &os = std::cout) const;

        /**
            @return The total number of tasks
        */
        static uli ntasks();

        /**
            @return The contraction engine, configured for particular matrix transpose types,
                    that will be used for actually performing the matrix multiplication
        */
        ContractionEngine* get_engine() const;

        ContractionConfiguration* get_configuration(uli threadnum) const;

        tensor_position_t get_alpha_type() const;

};

struct ContractionEngine {

    public:
        virtual void contract(
            DataNode* ldata,
            DataNode* rdata,
            DataNode* pdata,
            uli nrows,
            uli ncols,
            uli nlink,
            double scale,
            TemplateInfo::type_t dtype
        ) = 0;

};


}

#ifdef redefine_size_t
#undef size_t
#endif

#endif
