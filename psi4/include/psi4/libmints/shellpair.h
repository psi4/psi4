#ifndef _psi_src_lib_libmints_shellpair_h
#define _psi_src_lib_libmints_shellpair_h

// This should be included from libint2 itself eventually, but a workaround is to include it here
#include <system_error>
#include "psi4/libmints/basisset.h"
#include "libint2/shell.h"
#include "libint2/engine.h"

namespace psi {

    typedef std::vector<std::pair<int, int>> ShellPairBlock;
    using ShellPairData = std::vector<std::shared_ptr<libint2::ShellPair>>;

    class ShellPair {

    public:
    //! Form shell pair information
    static std::tuple<std::vector<ShellPairBlock>, ShellPairData> build_shell_pair_list(std::shared_ptr<BasisSet> bs1,
                                                                                 std::shared_ptr<BasisSet> bs2,
                                                                                 double threshold) {
        const auto nsh1 = bs1->nshell();
        const auto nsh2 = bs2->nshell();
        const auto bs1_equiv_bs2 = (bs1 == bs2);
        std::vector<ShellPairBlock> blocks;
        auto nthreads = Process::environment.get_n_threads();

        // construct the 2-electron repulsion integrals engine
        using libint2::Engine;
        std::vector<Engine> engines;
        engines.reserve(nthreads);
        engines.emplace_back(libint2::Operator::overlap, std::max(bs1->max_nprimitive(), bs2->max_nprimitive()),
                             std::max(bs1->max_am(), bs2->max_am()), 0);
        for (size_t i = 1; i != nthreads; ++i) {
            engines.push_back(engines[0]);
        }

        std::mutex mx;

#pragma omp parallel
        {
            int thread_id = 0;
#ifdef _OPENMP
            thread_id = omp_get_thread_num();
#endif
            auto& engine = engines[thread_id];
            const auto& buf = engine.results();

            // loop over permutationally-unique set of shells
            for (auto s1 = 0l, s12 = 0l; s1 != nsh1; ++s1) {
                auto n1 = bs1->shell(s1).nfunction();

                auto s2_max = bs1_equiv_bs2 ? s1 : nsh2 - 1;
                for (auto s2 = 0; s2 <= s2_max; ++s2, ++s12) {
                    if (s12 % nthreads != thread_id) continue;

                    auto on_same_center = (bs1->shell(s1).center() == bs2->shell(s2).center());
                    bool significant = on_same_center;
                    if (!on_same_center) {
                        auto n2 = bs2->shell(s2).nfunction();
                        engines[thread_id].compute(bs1->l2_shell(s1), bs2->l2_shell(s2));
                        double normsq = std::inner_product(buf[0], buf[0] + n1 * n2, buf[0], 0.0);
                        auto norm = std::sqrt(normsq);
                        significant = (norm >= threshold);
                    }

                    if (significant) {
                        auto block = bs1_equiv_bs2 && bs2->shell(s2).am() > bs1->shell(s1).am()
                                         ? ShellPairBlock{{s2, s1}}
                                         : ShellPairBlock{{s1, s2}};
                        mx.lock();
                        blocks.push_back(block);
                        mx.unlock();
                    }
                }
            }
        }  // end of compute

        // resort shell list in increasing order of angular momentum
        std::sort(blocks.begin(), blocks.end(), [&](auto& l, auto& r) { 
            const auto& lsh1 = bs1->shell(l[0].first);
            const auto& lsh2 = bs2->shell(l[0].second);
            const auto& rsh1 = bs1->shell(r[0].first);
            const auto& rsh2 = bs2->shell(r[0].second);
            const auto lam = lsh1.am() + lsh2.am();
            const auto ram = rsh1.am() + rsh2.am();
            return lam < ram;
        });

        const auto max_engine_precision = std::numeric_limits<double>::epsilon() / 1e10;

        // compute shellpair data assuming that we are computing to default_epsilon
        size_t npairs = blocks.size();
        ShellPairData spdata(npairs);
#pragma omp parallel
        {
            int thread_id = 0;
#ifdef _OPENMP
            thread_id = omp_get_thread_num();
#endif
            for (int pair = 0; pair < npairs; ++pair) {
                if (pair % nthreads == thread_id) {
                    auto s1 = blocks[pair][0].first;
                    auto s2 = blocks[pair][0].second;
                    spdata[pair] = std::make_shared<libint2::ShellPair>(bs1->l2_shell(s1), bs2->l2_shell(s2),
                                                                        std::log(max_engine_precision));
                }
            }
        }
        return std::make_tuple(blocks, spdata);
    }
    };

} // namespace psi4
#endif /* header guard */
