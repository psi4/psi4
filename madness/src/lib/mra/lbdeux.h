/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680

  $Id: lbdeux.h 2173 2011-02-23 21:40:46Z justus.c79@gmail.com $
*/
#ifndef MADNESS_MRA_IBDEUX_H__INCLUDED
#define MADNESS_MRA_IBDEUX_H__INCLUDED

#include <mra/mra.h>
#include <madness_config.h>
#include <map>
#include <queue>

/// \file mra/lbdeux.h
/// \brief Implements (2nd generation) static load/data balancing for functions
/// \ingroup function

namespace madness {

    template <std::size_t NDIM>
    class LBDeuxPmap : public WorldDCPmapInterface< Key<NDIM> > {
        typedef Key<NDIM> keyT;
        typedef std::pair<keyT,ProcessID> pairT;
        typedef std::map<keyT,ProcessID> mapT;
        mapT map;
        typedef typename mapT::const_iterator iteratorT;

    public:
        LBDeuxPmap(const std::vector<pairT>& v) {
            for (unsigned int i=0; i<v.size(); ++i) {
                map.insert(v[i]);
            }
        }

        // If level 0 is not entered as a node this will
        // be an infinite loop.
        ProcessID owner(const keyT& key) const {
            while (key.level() >= 0) {
                iteratorT it = map.find(key);
                if (it == map.end()) {
                    return owner(key.parent());
                }
                else {
                    return it->second;
                }
            }
            madness::print("Mon Dieux!", key);
            throw "LBDeuxPmap: lookup failed";
        }

        void print() const {
            madness::print("LBDeuxPmap");
        }
    };



    template <std::size_t NDIM>
    class LBNodeDeux {
        static const int nchild = (1<<NDIM);
        typedef Key<NDIM> keyT;
        typedef LBNodeDeux<NDIM> nodeT;
        typedef WorldContainer<keyT,nodeT> treeT;
        volatile double child_cost[nchild];
        volatile double my_cost;
        volatile double total_cost;
        volatile bool gotkids;
        volatile int nsummed;

        /// Computes index of child key in this node using last bit of translations
        int index(const keyT& key) {
            int ind = 0;
            for (std::size_t d=0; d<NDIM; ++d) ind += ((key.translation()[d])&0x1) << d;
            return ind;
        }

    public:
        LBNodeDeux()
                : my_cost(0.0), total_cost(0.0), gotkids(false), nsummed(0) {
            for (int i=0; i<nchild; ++i) child_cost[i] = 0.0;
        }

        bool has_children() const {
            return gotkids;
        }

        double get_total_cost() const {
            return total_cost;
        }

        /// Accumulates cost into this node
        Void add(double cost, bool got_kids) {
            total_cost = (my_cost += cost);
            gotkids = gotkids || got_kids;
            return None;
        }

        /// Accumulates cost up the tree from children
        Void sum(const treeT& tree, const keyT& child, double value) {
            child_cost[index(child)] = value;
            ++nsummed;
            if (nsummed == nchild) {
                for (int i=0; i<nchild; ++i) total_cost += child_cost[i];
                if (child.level() > 1) {
                    keyT key = child.parent();
                    keyT parent = key.parent();
                    const_cast<treeT&>(tree).task(parent, &nodeT::sum, tree, key, double(total_cost));
                }
            }
            return None;
        }


        /// Logically deletes this node by setting cost to -1

        /// Cannot actually erase this node from the container since the send() handler
        /// is holding an accessor to it.
        Void deleter(const treeT& tree, const keyT& key) {
            total_cost = my_cost = -1.0;
            if (has_children()) {
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT child = kit.key();
                    const_cast<treeT&>(tree).task(child, &nodeT::deleter, tree, child);
                }
            }
            return None;
        }

        /// Descends tree deleting all except internal nodes and sub-tree parents
        Void partition(const treeT& tree, const keyT& key, double avg) {
            if (has_children()) {
                // Sort children in descending cost order
                keyT keys[nchild];
                double vals[nchild];
                for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                    const keyT child = kit.key();
                    int ind = index(child);
                    keys[ind] = child;
                    vals[ind] = child_cost[ind];
                }
                for (int i=0; i<nchild; ++i) {
                    for (int j=i+1; j<nchild; ++j) {
                        if (vals[i] < vals[j]) {
                            std::swap(vals[i],vals[j]);
                            std::swap(keys[i],keys[j]);
                        }
                    }
                }

                // Split off subtrees in decreasing cost order
                for (int i=0; i<nchild; ++i) {
                    if (total_cost <= avg) {
                        const_cast<treeT&>(tree).task(keys[i], &nodeT::deleter, tree, keys[i]);
                    }
                    else {
                        total_cost -= vals[i];
                        const_cast<treeT&>(tree).task(keys[i], &nodeT::partition, tree, keys[i], avg);
                    }
                }
            }
            return None;
        }

        template <typename Archive>
        void serialize(Archive& ar) {
            ar & archive::wrap_opaque(this,1);
        }
    };


    template <std::size_t NDIM>
    class LoadBalanceDeux {
        typedef Key<NDIM> keyT;
        typedef LBNodeDeux<NDIM> nodeT;
        typedef WorldContainer<keyT,nodeT> treeT;
        typedef typename treeT::iterator iteratorT;
        typedef typename treeT::const_iterator const_iteratorT;
        World& world;
        treeT tree;


        template <typename T, typename costT>
        struct add_op {
            LoadBalanceDeux* lb;
            const costT& costfn;
            add_op(LoadBalanceDeux* lb, const costT& costfn) : lb(lb), costfn(costfn) {}
            void operator()(const keyT& key, const FunctionNode<T,NDIM>& node) const {
                lb->tree.send(key, &nodeT::add, costfn(key,node), node.has_children());
            }
        };

        /// Sums costs up the tree returning to everyone the total cost
        double sum() {
            world.gop.fence();
            const_iteratorT end = tree.end();
            for (const_iteratorT it=tree.begin(); it!=end; ++it) {
                const keyT& key = it->first;
                const nodeT& node = it->second;
                if (!node.has_children() && key.level() > 0) {
                    tree.task(key.parent(), &nodeT::sum, tree, key, node.get_total_cost());
                }
            }
            world.gop.fence();
            double total;
            keyT key0(0);
            if (world.rank() == tree.owner(key0)) {
                total = tree.find(key0).get()->second.get_total_cost();
            }
            world.gop.broadcast(total, tree.owner(key0));
            world.gop.fence();

            return total;
        }

        /// Used to sort results into descending order
        static bool compare(const std::pair<keyT,double>& a, const std::pair<keyT,double>& b) {
            return a.second < b.second;
        }


    public:
        LoadBalanceDeux(World& world)
                : world(world)
                , tree(world, FunctionDefaults<NDIM>::get_pmap()) {
            world.gop.fence();
        };

        /// Accumulates cost from a function
        template <typename T, typename costT>
        void add_tree(const Function<T,NDIM>& f, const costT& costfn, bool fence=false) {
            const_cast<Function<T,NDIM>&>(f).unaryop_node(add_op<T,costT>(this,costfn), fence);
        }

        /// Printing for the curious
        void print_tree(const keyT& key = keyT(0)) {
            Future<iteratorT> futit = tree.find(key);
            iteratorT it = futit.get();
            if (it != tree.end()) {
                for (int i=0; i<key.level(); ++i) std::cout << "  ";
                print(key, it->second.get_total_cost());

                if (it->second.has_children()) {
                    for (KeyChildIterator<NDIM> kit(key); kit; ++kit) {
                        print_tree(kit.key());
                    }
                }
            }
        }

        struct CostPerProc {
            double cost;
            int proc;
            CostPerProc() : cost(0.0), proc(0) {}
            CostPerProc(double cost, int proc) : cost(cost), proc(proc) {}
            bool operator<(const CostPerProc& other) const {
                return cost > other.cost;  // Want ascending order
            }
        };

        /// Actually does the partitioning of the tree
        std::shared_ptr< WorldDCPmapInterface<keyT> > load_balance(double fac = 1.0, bool printstuff=false) {
            world.gop.fence();
            // Compute full tree of costs
            double avg = sum()/(world.size()*fac);
            //if (world.rank() == 0) print_tree();
            world.gop.fence();

            // Create partitioning
            keyT key0(0);
            if (world.rank() == tree.owner(key0)) {
                tree.send(key0, &nodeT::partition, tree, key0, avg*1.1);
            }
            world.gop.fence();

            // Collect entire vector onto node0
            std::vector< std::pair<keyT,double> > results;
            const_iteratorT end = tree.end();
            for (const_iteratorT it=tree.begin(); it!=end; ++it) {
                if (it->second.get_total_cost() >= 0) {
                    results.push_back(std::make_pair(it->first,it->second.get_total_cost()));
                }
            }
            results = world.gop.concat0(results, 128*1024*1024);
            world.gop.fence();

            std::vector< std::pair<keyT,ProcessID> > map;

            if (world.rank() == 0) {

                std::sort(results.begin(), results.end(), compare);
                if (printstuff) {
                    print("THESE ARE THE INITIAL SUBTREES");
                    for (unsigned int i=0; i<results.size(); ++i) print(i,results[i]);
                }

                // Now use bin packing to cram the results together
                map.reserve(results.size());

                // Shove the first nproc entries directly into the queue
                unsigned int nproc = world.size();
                std::priority_queue<CostPerProc> costs;
                for (unsigned int p=0; p<nproc && !results.empty(); ++p) {
                    const std::pair<keyT,double>& f = results.back();
                    costs.push(CostPerProc(f.second,p));
                    map.push_back(std::make_pair(f.first,p));
                    results.pop_back();
                }

                // Process the remainder using the sorting maintained by the priority queue
                while (!results.empty()) {
                    const std::pair<keyT,double>& f = results.back();
                    CostPerProc top = costs.top();
                    costs.pop();
                    top.cost += f.second;
                    costs.push(top);
                    map.push_back(std::make_pair(f.first,top.proc));
                    results.pop_back();
                }
                if (printstuff) {
                    print("THIS IS THE MAP");
                    print(map);
                    print("THESE ARE THE COSTS PER PROCESSOR");
                    while (!costs.empty()) {
                        print(costs.top().proc,costs.top().cost);
                        costs.pop();
                    }
                }
            }

            world.gop.fence();
            world.gop.broadcast_serializable(map, 0);
            world.gop.fence();

            // Return the Procmap

            return std::shared_ptr< WorldDCPmapInterface<keyT> >(new LBDeuxPmap<NDIM>(map));
        }
    };
}


#endif // MADNESS_MRA_IBDEUX_H__INCLUDED

