#ifndef yeti_contraction_h
#define yeti_contraction_h

#include "class.h"
#include "taskqueue.h"
#include "yetiobject.h"

#include "taskqueue.hpp"
#include "matrix.hpp"
#include "contraction.hpp"
#include "permutation.hpp"
#include "tile.hpp"
#include "tensor.hpp"
#include "data.hpp"
#include "mallocimpl.h"

namespace yeti {

/**
    @class ContractionTask
    Task instance for accumulating product of two matrices into a product matrix
*/
class ContractionTask :
    public Task
{

    private:
        Matrix* lmatrix_;

        Matrix* rmatrix_;

        Matrix* product_matrix_;

        Contraction* cxn_;


    public:
        /**
            @param lmatrix
            @param rmatrix
            @param product_matrix
            @param cxn
        */
        ContractionTask(
            Matrix* lmatrix,
            Matrix* rmatrix,
            Matrix* product_matrix,
            Contraction* cxn
        );

        virtual ~ContractionTask();

        void print(std::ostream& os = std::cout) const;

        /**
            Polymorphic function called by TaskQueue
            @param threadnum The thread number
        */
        void run(uli threadnum);

};

/**
    @class Contraction
    Class used for configuring a contraction.
*/
class Contraction :
    public TaskParent
{

    public:
        typedef enum { left_tensor, right_tensor, product_tensor } tensor_position_t;

    private:
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
        MatrixConfigurationPtr product_config_;

        /**
            The "dominant" tensor that will be the owner of all tasks generated
            is at this position.  Communication/recomputation of this tensor
            will be minimized.
        */
        tensor_position_t alphatype_;

        /**
            The "dominant" tensor that will be the owner of all tasks generated.
            Communication/recomputation of this tensor
            will be minimized.
        */
        TensorPtr alpha_tensor_;

        /**
            The final permutational symmetry of the target tensor
            accounting for all symmetries:
            Those imposed by the user, those automatically obtained
            by symmetry of left and right tensors, and special symmetries.
        */
        PermutationGroupPtr final_grp_;

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
            The permutation group the product tensor is required
            to have.
        */
        PermutationGroupPtr required_grp_;

        /**
            The set of permutations that must be applied to the
            product tensor in order to make it have the symmetry
            in #required_grp_
        */
        PermutationSetPtr product_set_;

        /**
            A symmetry in the target indices that does not exist
            in either the left or right tensor
        */
        PermutationSetPtr special_set_;

        /**
            The set of permutations common to the contraction
            indices for the left matrix and right matrix
        */
        PermutationGroupPtr contraction_grp_;

        TensorPtr ltensor_;

        TensorPtr rtensor_;

        TensorPtr product_tensor_;

        MatrixPtr lmatrix_;

        MatrixPtr rmatrix_;

        MatrixPtr pmatrix_;

        /**
            Because all data types in a contraction are dynamic, certain values
            may need to be reallocated on the fly by casting up (float->double) or
            down (double->float).  This temporary space provides a static workspace
            for putting cast values without modifiying the original data. This
            tempory space is for the left tensor
        */
        char* ltmp_;

        /**
            See #ltmp_. Temporary space for right matrix.
        */
        char* rtmp_;

        /**
            See #rtmp_. Temporary space for product matrix.
        */
        char* ptmp_;

        double scale_;

        typedef double max_datasize_type;

        /**
            Depending on how the matrices are to be interpreted, the contraction
            engine manipulates the data appropriately.  For example, you could
            multiply and of the following cases (^T denotes transpose)
            L * R, L^T * R, L * R^T, L^T * R^T.
        */
        ContractionEngine* engine_;

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
            const TensorPtr& ltensor,
            const TensorPtr& rtensor,
            const MatrixIndexPtr& lindex,
            const MatrixIndexPtr& rindex,
            const PermutationGroupPtr& required_grp,
            const PermutationPtr& lperm,
            const PermutationPtr& rperm,
            const TensorPtr& ptensor = 0 //tensor to accumulate into
        );

        virtual ~Contraction();

        /**
            Create a new contraction task and add it to the contraction queue
            @param lmatrix
            @param rmatrix
            @param product_matrix
        */
        void add_task(
            Matrix* lmatrix,
            Matrix* rmatrix,
            Matrix* product_matrix
        );

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

        /**
            @return A set of permutations that target indices will have that is not obviously
                    obtained.  For example, the symmetric nature of the exchange matrix
                    is not a simple, direct result of the permutational symmetry of the
                    integrals and density.
        */
        PermutationSet* get_special_set() const;

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

        Matrix* get_left_matrix() const;

        Matrix* get_right_matrix() const;

        Matrix* get_product_matrix() const;

        Tensor* get_left_tensor() const;

        Tensor* get_right_tensor() const;

        Tensor* get_product_tensor() const;

        /**
            @return The dominant tensor in the contraction, i.e. the one for which
                    communication/recomutation should be minimized
        */
        Tensor* get_alpha_tensor() const;

        void* left_tmp_space() const;

        void* right_tmp_space() const;

        void* product_tmp_space() const;

        /**
            Run all jobs on the contraction queue
        */
        static void run();

        /**
            Clear all jobs on the contraction queue
        */
        static void clear();

        /**
            Optimize job ordering for the most efficient possible contraction
        */
        static void configure();

        /**
        */
        static void print_cxn(std::ostream &os = std::cout);

        void print(std::ostream &os = std::cout) const;

        /**
            @return The number of task owners.  This is the number of distinct alpha tiles
            that are involved in the contraction.
        */
        static uli nowners();

        /**
            @return The total number of tasks
        */
        static uli ntasks();

        /**
            @return The contraction engine, configured for particular matrix transpose types,
                    that will be used for actually performing the matrix multiplication
        */
        ContractionEngine* get_engine() const;

};


}

#endif
