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
    public Task,
    public Malloc<ContractionTask>
{


    private:
        friend class Contraction;

        uli row_;

        uli col_;

        uli link_;

        TensorBlock* main_block_;
        
        bool follow_permutations_;

        Contraction* cxn_;

        void run(
            uli threadnum,
            TensorBlock* lblock,
            TensorBlock* rblock,
            TensorBlock* pblock
        );

        void run_product_driven(uli threadnum);

        void run_left_driven(uli threadnum);

        void run_right_driven(uli threadnum);

        void prefetch_product_driven(uli threadnum);

        void prefetch_left_driven(uli threadnum);

        void prefetch_right_driven(uli threadnum);

        void run_p_l_r(uli threadnum);

        void run_p_r_l(uli threadnum);

        void run_r_p_l(uli threadnum);

        void run_r_l_p(uli threadnum);

        void run_l_r_p(uli threadnum);

        void run_l_p_r(uli threadnum);

        void prefetch_p_l_r(uli threadnum);

        void prefetch_p_r_l(uli threadnum);

        void prefetch_r_p_l(uli threadnum);

        void prefetch_r_l_p(uli threadnum);

        void prefetch_l_r_p(uli threadnum);

        void prefetch_l_p_r(uli threadnum);

        TensorBlock* get_product_block(uli threadnum, uli row, uli col);

        TensorBlock* get_left_block(uli threadnum, uli row, uli link);

        TensorBlock* get_right_block(uli threadnum, uli link, uli col);

        ContractionTask(
            uli row,
            uli col,
            uli link,
            TensorBlock* main_block,
            Contraction* cxn
        );

    public:
        /**
            @param idx1 The row index of the main matrix
            @param idx2 The col index of the main matrix
            @param order The cxn ordering
        */
        ContractionTask(
            uli row,
            uli col,
            uli link,
            Contraction* cxn
        );

        virtual ~ContractionTask();

        void print(std::ostream& os = std::cout) const;

        /**
            Polymorphic function called by TaskQueue
            @param threadnum The thread number
        */
        void run(uli threadnum);

        void prefetch(uli threadnum);

        uli append_info(uli* data) const;

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

        TensorBlock* next_tmp_block_;

        Contraction* cxn_;

    public:
        ContractionConfiguration(
            const MatrixConfigurationPtr& lconfig,
            const MatrixConfigurationPtr& rconfig,
            const MatrixConfigurationPtr& pconfig,
            Contraction* cxn
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

        TensorBlock* get_tmp_block() const;

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

        bool all_cols_on_node(Tensor* tensor) const;

        bool all_rows_on_node(Tensor* tensor) const;

        bool transpose_right() const;

        bool transpose_left() const;

        void check_distribution(
            Tensor* ltensor,
            Tensor* rtensor,
            Tensor* ptensor,
            bool& ltensor_all_rows_local,
            bool& ltensor_all_cols_local,
            bool& rtensor_all_rows_local,
            bool& rtensor_all_cols_local,
            bool& ptensor_all_rows_local,
            bool& ptensor_all_cols_local
        );

};

/**
    @class Contraction
    Class used for configuring a contraction.
*/
class Contraction :
    public TaskParent,
    public Malloc<Contraction>
{

    public:
        typedef enum {
            left_tensor = 0,
            right_tensor = 1,
            product_tensor = 2
        } tensor_position_t;

    public:
        typedef enum { 
            p_l_r_cxn, 
            p_r_l_cxn,
            r_l_p_cxn,
            r_p_l_cxn,
            l_r_p_cxn,
            l_p_r_cxn
        } cxn_order_t;

    private:    
        MatrixIndexPtr lindex_;
        
        MatrixIndexPtr rindex_;
        
        MatrixIndexPtr pindex_;
        
        PermutationGroupPtr lcxn_grp_;
        
        PermutationGroupPtr rcxn_grp_;
        
        bool transpose_right_;
        
        bool transpose_left_;

        uli nrows_;
        
        uli ncols_;

        uli nlink_;

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

        uli col_;

        uli row_;

        uli link_;

        tensor_position_t distribution_type_;

        tensor_position_t alpha_type_;

        tensor_position_t beta_type_;

        tensor_position_t gamma_type_;

        bool use_replicated_product_;

        bool product_thread_clash_;

        bool product_node_clash_;

        bool all_tensors_replicated_;
        /**
            Depending on how the matrices are to be interpreted, the contraction
            engine manipulates the data appropriately.  For example, you could
            multiply and of the following cases (^T denotes transpose)
            L * R, L^T * R, L * R^T, L^T * R^T.
        */
        ContractionEngine* engine_;

        cxn_order_t order_;

        Task* get_next_l_r_p();
                                            
        Task* get_next_l_p_r();

        Task* get_next_r_l_p();

        Task* get_next_r_p_l();

        Task* get_next_p_l_r();

        Task* get_next_p_r_l();

        TensorBlock* get_product_block(uli r, uli c);

        TensorBlock* get_left_block(uli r, uli c);

        TensorBlock* get_right_block(uli r, uli c);

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

        void add_dynamic_task(ContractionTask* task);

        void finalize();

        Task* get_next_task();

        cxn_order_t order() const;

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

        void get_row_col_left(uli idx, uli& r, uli& c) const;

        void get_row_col_right(uli idx, uli& r, uli& c) const;

        void get_row_col_product(uli idx, uli& r, uli& c) const;

        uli get_left_index(uli r, uli c) const;

        uli get_right_index(uli r, uli c) const;

        uli get_product_index(uli r, uli c) const;

        Tensor* ptensor() const;

        Tensor* rtensor() const;

        Tensor* ltensor() const;

        bool is_thread_replicated_product() const;
        
        bool product_thread_clash() const;

        bool product_node_clash() const;
        
        /**
            Run all jobs on the contraction queue
        */
        static void run();

        uli nrows() const;

        uli ncols() const;

        uli nlink() const;

        /**
            Clear all jobs on the contraction queue
        */
        void print(std::ostream &os = std::cout) const;

        /**
            @return The contraction engine, configured for particular matrix transpose types,
                    that will be used for actually performing the matrix multiplication
        */
        ContractionEngine* get_engine() const;

        ContractionConfiguration* get_configuration(uli threadnum) const;

        tensor_position_t get_alpha_type() const;

        bool do_task(
            TensorBlock* lblock,
            TensorBlock* rblock,
            TensorBlock* pblock
        );

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
