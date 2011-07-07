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


  $Id: loadbal.h 2343 2011-06-03 17:11:25Z justus.c79@gmail.com $
*/

/// \file loadbal.h
/// \brief Declares and partially implements MyPmap, LoadBalImpl and associated load balancing classes.


#ifndef MADNESS_MRA_LOADBAL_H__INCLUDED
#define MADNESS_MRA_LOADBAL_H__INCLUDED

#include <world/world.h>
#include <misc/misc.h>
#include <tensor/tensor.h>
#include <misc/ran.h>
#include <mra/key.h>
#include <mra/funcimpl.h>

namespace madness {
    template <typename T, std::size_t NDIM> class FunctionImpl;

    typedef unsigned long Cost;
    typedef double CompCost;

    inline Cost default_cost_fun() {
        return 1;
    }


    /// Finds exponent k such that d^k <= me < d^{k+1}
    inline int nearest_power(int me, int d) {
        int k = 0;
        while (me != 0) {
            if (me%d == 0) {
                ++k;
                me/=d;
            }
            else {
                break;
            }
        }
        return k;
    }

    template <int D> class LBNode;
    template <int D> class TreeCoords;
    template <int D> class MyPmap;
    template <int D> class LBTree;
    class NodeData;
    template <int D> class LoadBalImpl;

    /// Convenient typedef shortcuts

    /// Makes it easier to handle these unwieldy templated types
    template <int D>
    struct DClass {
//     typedef Key<D> KeyD;
//     typedef const Key<D> KeyDConst;
//     //typedef TreeCoords<D> TreeCoords;
//     typedef LBNode<D> NodeD;
//     typedef const LBNode<D> NodeDConst;
//     typedef MyPmap<D> MyPmap;
//     typedef LBTree<D> treeT;
        typedef std::vector< std::vector< madness::TreeCoords<D> > > vvTreeCoords;
    };



    template <int D>
    class PartitionInfo {
    public:
        Cost maxcost;
        Cost cost_left;
        Cost skel_cost;
        unsigned int partition_number;
        unsigned int step_num;
        double facter;
        PartitionInfo(double f=1.1) :
                maxcost(0)
                , cost_left(0)
                , skel_cost(0)
                , partition_number(0)
                , step_num(0)
                , facter(f) { };

        void reset(unsigned int p=1) {
            maxcost = 0;
            cost_left = skel_cost;
            partition_number = p;
            ++step_num;
        }

        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & maxcost & cost_left & skel_cost & partition_number & step_num & facter;
        }
    };

    template <int D>
    std::ostream& operator<<(std::ostream& s, const PartitionInfo<D>& pi) {
        s << "maxcost = " << pi.maxcost << ", cost_left = " << pi.cost_left <<
        ", skel_cost = " << pi.skel_cost << ", partition_number = " <<
        pi.partition_number << ", step_num = " << pi.step_num <<
        ", facter = " << pi.facter;
        return s;
    }



    /// Diagnostic data contained in fascimile tree
    /// Diagnostic data, including the cost of the node and the subtree headed by that node,
    /// along with a bool flag used during depth-first partitioning
    class NodeData {
        friend std::ostream& operator<<(std::ostream& s, const NodeData& nd);
    public:
        int cost;
        int subcost;
        bool is_taken;
        NodeData(int c = 1, int s = 1, bool i = false) : cost(c), subcost(s), is_taken(i) {};
        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & cost & subcost & is_taken;
        }
        void print() {
            std::cout << "cost = " << cost << ", subcost = " << subcost << ", is_taken = " << is_taken << std::endl;
        }
        template <typename functionT>
        void set_data(functionT function) {
            functionT tmp = (*function);
            this->cost = (int) tmp;
        }

        template <typename functionT, typename arg1T>
        void set_data(functionT function, const arg1T& arg1) {
            functionT tmp = (*function)(arg1);
            this->cost = (int) tmp;
        }
    };


    inline std::ostream& operator<<(std::ostream& s, const NodeData& nd) {
        s << "cost " << nd.cost << ", subcost " << nd.subcost << ", is_taken " << nd.is_taken;
        return s;
    }


    /// The node that is used in the fascimile copy of the tree to be load balanced

    /// The node used in the tree that is operated upon and load balanced in LoadBalImpl.
    template <int D>
    class LBNode {
    private:
        NodeData data;
        std::vector<bool> c; /// Existence of each child individually

        void all_children(bool status=false) {
            c.clear();
            c.assign(dim, status);
        }

    public:
        mutable KeyChildIterator<D> rpit;
        static int dim; /// Number of children in standard tree (e.g. 2^D)
        int nrecvd;

        LBNode() {
            rpit = KeyChildIterator<D>();
            data = NodeData();
            all_children();
            nrecvd = 0;
        }

        LBNode(const NodeData& d, bool children=false, int n=0) : data(d), nrecvd(n) {
            rpit = KeyChildIterator<D>();
            all_children(children);
        }

        LBNode(const LBNode& node) : data(node.data), c(node.c), rpit(node.rpit), nrecvd(node.nrecvd) { };


        /// Determines whether node has any children at all
        bool has_children() const {
            for (int i = 0; i < dim; ++i)
                if (c[i]) return true;
            return false;
        }

        bool has_child(unsigned int i) const {
            return c[i];
        }

        bool has_child(int i) const {
            return c[i];
        }

        int get_num_children() const {
            int nkids = 0;
            for (int i=0; i < dim; ++i) {
                if (has_child(i)) ++nkids;
            }
            return nkids;
        }

        void set_child(int i, bool setto = true) {
            c[i] = setto;
        }

        void set_all_children(bool setto = true) {
            all_children(setto);
        }

        void set_data(const NodeData& d) {
            data = d;
        }

        template <typename functionT>
        void set_cost(functionT function) {
            data.set_data<functionT>(function);
        }

        template <typename functionT, typename arg1T>
        void set_cost(functionT function, const arg1T& arg1) {
            data.set_data<functionT, arg1T>(function, arg1);
        }

        NodeData get_data() const {
            return data;
        }

        std::vector<bool> get_c() const {
            return c;
        }

        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & data & c & rpit;
        }
    };


    template <int D>
    std::ostream& operator<<(std::ostream& s, const LBNode<D>& node) {
        s << "data = " << node.get_data() << ", c = " << node.get_c();
        if (node.rpit) {
            s  << ", key_iterator = " << node.rpit.key();
        }
        else {
            s << ", key_iterator = <EMPTY>";
        }
        return s;
    }

//     template <int D>
//     std::ostream& operator<<(std::ostream& s, const LBNode<D>& node) {
//         s << "data = " << node.get_data() << ", c = " << node.get_c();
// 	if (node.rpit) {
// 	  s << ", key_iterator = " << node.rpit.key();
// 	} else {
// 	  s << ", key_iterator = <EMPTY>";
// 	}
//         return s;
//     };


    template <int D>
    int LBNode<D>::dim = power<D>();


    /// Key + owner, struct used to determine mapping of tree nodes
    template <int D>
    class TreeCoords {
    public:
        Key<D> key;
        ProcessID owner;

        TreeCoords(const Key<D>& k, ProcessID o) : key(Key<D>(k)), owner(o) {};
        TreeCoords(const TreeCoords& t) : key(Key<D>(t.key)), owner(t.owner) {};
        TreeCoords() : key(Key<D>()), owner(-1) {};
        void print() const {
            madness::print(key, "   owner =", owner);
        }

        bool operator< (const TreeCoords t) const {
            return (this->key < t.key);
        }

        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & key & owner;
        }
    };

    template <int D>
    inline std::ostream& operator<<(std::ostream& s, const TreeCoords<D>& tc) {
        s << tc.key << "   owner = " << tc.owner;
        return s;
    }



    template<int D>
    class ProcMapImpl {
    public:
        typedef std::map< const Key<D>,ProcessID> Mapinfo;
        typedef typename Mapinfo::iterator iterator;
        typedef const iterator iterator_const;
        typedef std::pair< const Key<D>, ProcessID > pairT;

        ProcMapImpl() {};
        ProcMapImpl(std::vector< TreeCoords<D> > v) {
            int vlen = v.size();
            for (int i = 0; i < vlen; ++i) {
                themap.insert(std::make_pair(v[i].key, v[i].owner));
            }
        }

        ProcMapImpl(const TreeCoords<D>& t) {
            themap.insert(std::make_pair(t.key, t.owner));
        }
        void insert(const TreeCoords<D>& t) {
            themap.insert(std::make_pair(t.key, t.owner));
        }
        void erase(const TreeCoords<D>& t) {
            themap.erase(t.key);
        }

        ProcessID find_owner(const Key<D>& key) const {
            typename std::map< const Key<D>,ProcessID>::const_iterator it = themap.find(key);
            if (it != themap.end()) {
                return it->second;
            }
            else if (key.level() == 0) {
                madness::print("find_owner: owner of ", key, "not found but returning 0");
                return 0;
            }
            else {
                return this->find_owner(key.parent());
            }
        }

        void print() {
            for (iterator it = themap.begin(); it != themap.end(); ++it) {
                madness::print(it->first, "   ", it->second);
            }
        }

    private:
        Mapinfo themap;
    };


    /// Procmap implemented using Tree of TreeCoords

    template <int D>
    class MyPmap : public WorldDCPmapInterface< Key<D> > {
    private:
        unsigned int map_type; // 0 = simple map, 1 = gaussian distributed map, 2 = treecoords list
        const int nproc;
        const int n;
        std::shared_ptr< ProcMapImpl<D> > tree_map; // for map_type 2
        Tensor<ProcessID> simple_key_map; // map of keys at level n, for map_type 1
        typedef Key<D> KeyD;

        /// private method that builds the Tree underlying the procmap
        void build_tree_map(std::vector< TreeCoords<D> > v) {
            tree_map = std::shared_ptr< ProcMapImpl<D> > (new ProcMapImpl<D>(v));
        }

        ProcessID simple_hash(const KeyD& key) const {
	    if (key.level() == 0) return 0;
            KeyD parent = (key.level() > n) ? key.parent(key.level()-n) : key;
            return (parent.hash()%nproc);
        }
        ProcessID not_so_simple_hash(const KeyD& key) const {
            KeyD parent = (key.level() > n) ? key.parent(key.level()-n) : key;
            return simple_key_map((const long *) &(parent.translation()[0]));
        }

        void prepare_not_so_simple_map(World& world) {
	    std::vector<long> vdim(D);
            for (int i=0; i<D; ++i) vdim[i] = 1L<<n;
            simple_key_map = Tensor<ProcessID>(vdim);

            std::list< std::pair<KeyD,double> > costmap;
            Vector<Translation,D> l;
            long cent = (1L<<n) / 2;
            for (TensorIterator<ProcessID> iter=simple_key_map.unary_iterator(0,false,false); iter._p0; ++iter) {
                double dist = 0.01;
                for (int i=0; i<D; ++i) {
                    l[i] = iter.ind[i];
                    dist += (l[i] - cent)*(l[i] - cent);
                }
                double cost = 1.0/dist; // actually dist squared
                cost *= (1.0 + 0.001*RandomValue<double>()); // To shuffle (nearly) equal values
                costmap.push_back(std::pair<KeyD,double>(KeyD(n,l),cost));
            }
            costmap.sort(costmapcmp);
//            if (world.rank() == 0) {
//                for (typename std::list< std::pair<KeyD,double> >::iterator it=costmap.begin(); it!=costmap.end(); ++it) {
//                    madness::print("costmap", it->first, it->second);
//                }
//            }
            ProcessID p = 0;
            for (typename std::list< std::pair<KeyD,double> >::iterator it=costmap.begin(); it!=costmap.end(); ++it) {
                const long *l = (const long *) &(it->first.translation()[0]);
                simple_key_map(l)  = p;
                ++p;
                if (p == world.size()) p = 0;
            }
//            if (world.rank() == 0) {
//                madness::print("SIMPLE MAP", D,"\n", simple_key_map);
//            }
        }

    public:
        MyPmap() : map_type(2) {};

        static bool costmapcmp(const std::pair<KeyD,double>& a, const std::pair<KeyD,double>& b) {
            return a.second > b.second;
        }
        MyPmap(World& world)
                : map_type(1)
                , nproc(world.nproc())
                , n(int((std::log((double)world.size())/std::log(2.0)+3)/D) + 2) { // 16*nproc = 2^(nD)
            // We set n to have about 16 tasks per processor and we try to
            // give each process a mix of large, medium, and small
            // tasks.  Currently estimate cost as inversely
            // proportional to distance from center but we could
            // enable the user to provide a function.

            //if (world.rank() == 0) madness::print("DIM",D,"N IN MAP IS",n);
            prepare_not_so_simple_map(world);
        }

        MyPmap(World& world, unsigned int map_type, int n=100)
                : map_type(map_type)
                , nproc(world.nproc())
                , n(n) {
            if (map_type==1) {
                n =int((std::log((double)world.size())/std::log(2.0)+3)/D) + 2; // 16*nproc = 2^(nD)
                prepare_not_so_simple_map(world);
            }
        }

        MyPmap(World& world, std::vector<TreeCoords<D> > v) : map_type(2), nproc(world.nproc()), n(0) {
            build_tree_map(v);
        }

        MyPmap(const MyPmap<D>& other) : map_type(other.map_type), nproc(other.nproc), n(other.n), tree_map(other.tree_map) {};

        MyPmap<D>& operator=(const MyPmap<D>& other) {
            if (this != &other) {
                map_type = other.map_type;
		simple_key_map = other.simple_key_map; // shallow copy
                nproc = other.nproc;
                n = other.n;
                tree_map = other.tree_map;
            }
            return *this;
        }

        void print() const {
            if (map_type == 2) {
                tree_map->print();
            } else if (map_type == 1) {
	        madness::print("MyPmap: gaussian distributed map with n =", n);
	    } else {
                madness::print("MyPmap: simple map with n =", n);
            }
        }

        /// Find the owner of a given key
        ProcessID owner(const KeyD& key) const {
	    if (map_type == 0) {
                return simple_hash(key);
	    } else if (map_type == 1) {
	        return not_so_simple_hash(key);
            } else {
                return tree_map->find_owner(key);
            }
        }
    };


//     template <int D>
//     class MyPmap : public WorldDCPmapInterface< Key<D> > {
//     private:
//         Tensor<ProcessID> simple_key_map; // map of keys at level n
//         bool simplemap;
//         const int nproc;
//         const ProcessID me;
//         const int n;
//         std::shared_ptr< ProcMapImpl<D> > tree_map;
//         typedef Key<D> KeyD;

//         /// private method that builds the Tree underlying the procmap
//         void build_tree_map(std::vector< TreeCoords<D> > v) {
//             tree_map.reset(new ProcMapImpl<D>(v));
//         };

//         ProcessID simple_hash(const KeyD& key) const {
//             KeyD parent = (key.level() > n) ? key.parent(key.level()-n) : key;
//             return simple_key_map((const long *) &(parent.translation()[0]));
//         };

//     public:
//         MyPmap() : simplemap(false) {};

//         static bool costmapcmp(const std::pair<KeyD,double>& a, const std::pair<KeyD,double>& b) {
//             return a.second > b.second;
//         };

//         MyPmap(World& world)
//                 : simplemap(true)
//                 , nproc(world.nproc())
//                 , me(world.rank())
//                 , n(int((std::log((double)world.size())/std::log(2.0)+3)/D) + 2) { // 16*nproc = 2^(nD)
//             // We set n to have about 16 tasks per processor and we try to
//             // give each process a mix of large, medium, and small
//             // tasks.  Currently estimate cost as inversely
//             // proportional to distance from center but we could
//             // enable the user to provide a function.

//             //if (world.rank() == 0) madness::print("DIM",D,"N IN MAP IS",n);

//             std::vector<long> vdim(D);
//             for (int i=0; i<D; ++i) vdim[i] = 1L<<n;
//             simple_key_map = Tensor<ProcessID>(vdim);

//             std::list< std::pair<KeyD,double> > costmap;
//             Vector<Translation,D> l;
//             long cent = (1L<<n) / 2;
//             for (TensorIterator<ProcessID> iter=simple_key_map.unary_iterator(0,false,false); iter._p0; ++iter) {
//                 double dist = 0.01;
//                 for (int i=0; i<D; ++i) {
//                     l[i] = iter.ind[i];
//                     dist += (l[i] - cent)*(l[i] - cent);
//                 }
//                 double cost = 1.0/dist; // actually dist squared
//                 cost *= (1.0 + 0.001*RandomValue<double>()); // To shuffle (nearly) equal values
//                 costmap.push_back(std::pair<KeyD,double>(KeyD(n,l),cost));
//             }
//             costmap.sort(costmapcmp);
// //             if (world.rank() == 0) {
// //                 for (typename std::list< std::pair<KeyD,double> >::iterator it=costmap.begin(); it!=costmap.end(); ++it) {
// //                     madness::print("costmap", it->first, it->second);
// //                 }
// //             }
//             ProcessID p = 0;
//             for (typename std::list< std::pair<KeyD,double> >::iterator it=costmap.begin(); it!=costmap.end(); ++it) {
//                 const long *l = (const long *) &(it->first.translation()[0]);
//                 simple_key_map(l)  = p;
//                 ++p;
//                 if (p == world.size()) p = 0;
//             }
// //             if (world.rank() == 0) {
// //                 madness::print("SIMPLE MAP", D,"\n", simple_key_map);
// //             }
//         };

//         MyPmap(World& world, vector<TreeCoords<D> > v) : simplemap(false), nproc(world.nproc()), me(world.rank()), n(0) {
//             build_tree_map(v);
//         };

//         MyPmap(const MyPmap<D>& other) : simplemap(other.staticmap), nproc(other.nproc), me(other.me), n(other.n), tree_map(other.tree_map) {};

//         MyPmap<D>& operator=(const MyPmap<D>& other) {
//             if (this != &other) {
//                 simple_key_map = other.simple_key_map; // shallow copy
//                 simplemap = other.simplemap;
//                 nproc = other.nproc;
//                 me = other.me;
//                 n = other.n;
//                 tree_map = other.tree_map;
//             }
//             return *this;
//         };

//         void print() const {
//             if (!simplemap) {
//                 tree_map->print();
//             }
//             else {
//                 madness::print("MyPmap: simple map with n =", n);
//             }
//         };

//         /// Find the owner of a given key
//         ProcessID owner(const KeyD& key) const {
//             if (simplemap)
//                 return simple_hash(key);
//             else {
//                 return tree_map->find_owner(key);
//             }
//         };
//     };

    /// The container in which the fascimile tree with its keys mapping to LBNodes is stored
    template <int D>
    class LBTree : public WorldObject< LBTree<D> > {
    public:
        typedef WorldObject<LBTree<D> > woT;
        typedef WorldContainer<Key<D>,LBNode<D> > dcT;
        typedef typename dcT::iterator iterator;
        typedef typename dcT::const_iterator const_iterator;

        World& world;

        static const Key<D> root;
        static typename DClass<D>::vvTreeCoords list_of_list;
        static std::vector<Cost> cost_list;

        PartitionInfo<D> partition_info;
        //std::vector<typename TreeCoords<D> > temp_list;
        std::vector< TreeCoords<D> > temp_list;
        Cost(*cost_fun)();

    private:
        dcT impl;

    public:
        LBTree(World& world, const std::shared_ptr< WorldDCPmapInterface< Key<D> > >& pmap, Cost(*cost_f)()=&default_cost_fun) : woT(world)
                , world(world)
                //, cost_fun(cost_fun)
                , impl(world,pmap) {
            impl.process_pending();
            this->process_pending();
            this->cost_fun=cost_f;
        }

        virtual ~LBTree() { }

        /// Initialize the LBTree by converting a FunctionImpl to a LBTree
        template <typename T, typename costfunT>
        inline void init_tree(const std::shared_ptr< FunctionImpl<T,D> >& f, const costfunT& costfun) {
            typename FunctionImpl<T,D>::dcT::const_iterator end = f->coeffs.end();
            for (typename FunctionImpl<T,D>::dcT::const_iterator it = f->coeffs.begin(); it != end; ++it) {
                // convert Node to LBNode
                NodeData nd;
                const Key<D>& key = it->first;
                const typename FunctionImpl<T,D>::nodeT&  node = it->second;

                //		nd.cost = Cost(1e6*costfun(key, node));
                nd.cost = Cost(costfun(key, node));
                nd.subcost = nd.cost;
                LBNode<D> lbnode(nd,node.has_children());
                impl.replace(key, lbnode);

//             	  if (!(it->second.has_children())) {
//                    nd.cost = (*cost_fun)();
//                    nd.subcost = nd.cost;
//                    LBNode<D> lbnode(nd,false);
// 	                  // insert into impl
//                    impl.insert(key, lbnode);
//                } else {
// 	                  nd.cost = (*cost_fun)();
// 	                  nd.subcost = nd.cost;
// 	                  LBNode<D> lbnode(nd,true);
// 	                  // insert into impl
// 	                  impl.insert(key, lbnode);
//                }
            }
        }


        // Methods:

        template <typename T, typename costfunT>
        inline void add_tree(const std::shared_ptr< FunctionImpl<T,D> >& f, const costfunT& costfun) {
            typename FunctionImpl<T,D>::dcT::const_iterator end = f->coeffs.end();
            for (typename FunctionImpl<T,D>::dcT::const_iterator it = f->coeffs.begin(); it != end; ++it) {
                // convert Node to LBNode
                NodeData nd;
                Key<D> key = it->first;
                typename LBTree<D>::const_iterator tree_it = impl.find(key);
                Cost new_cost = Cost(costfun(it->first,it->second));
                //		Cost new_cost = Cost(1e6*costfun(it->first,it->second));
                if (tree_it != impl.end()) {
                    LBNode<D> lbnode = tree_it->second;
                    if (it->second.has_children()) {
                        lbnode.set_all_children(true);
                    }
                    NodeData nd=lbnode.get_data();
                    nd.cost+=new_cost;
                    nd.subcost+=new_cost;
                    lbnode.set_data(nd);
                    impl.replace(key, lbnode);
                }
                else {
                    nd.cost = new_cost;
                    nd.subcost = nd.cost;
                    LBNode<D> lbnode(nd,it->second.has_children());
                    impl.replace(key, lbnode);
                }
            }
        }

        ProcessID owner(const Key<D>& key) {
            return impl.owner(key);
        }

        void print(const Key<D>& key) {
            typename LBTree<D>::iterator it = impl.find(key);
            if (it == impl.end()) return;
            for (Level i = 0; i < key.level(); ++i) std::cout << "  ";
            madness::print(key, it->second);
            for (KeyChildIterator<D> kit(key); kit; ++kit) {
                print(kit.key());
            }
        };

        void find_partitions(PartitionInfo<D>& pi);
        bool verify_partition(std::vector<TreeCoords<D> >& part_list);
        Void launch_make_partition(PartitionInfo<D> pi, bool first_time);
        Void meld_all(bool first_time=false);

        Cost fix_cost();


        void init_fix_cost();
        void fix_cost_spawn();
        Void fix_cost_sum(const Key<D>& key, Cost c);


        void rollup();

        void reset(bool taken);

        void meld(LBTree<D>::iterator it);

    private:
        // This function was created because there is a bug in the Intel 12.0.0
        // compiler that prevents it from correctly generating the type for
        // member function pointers in template classes
        Void make_partition_internal(const Key<D>& key, Cost partition_size,
                            Cost used_up, PartitionInfo<D> pi, bool downward);
    public:

        Void make_partition(const Key<D>& key, Cost partition_size,
                            Cost used_up, PartitionInfo<D> pi, bool downward = false) {
            return make_partition_internal(key, partition_size, used_up, pi, downward);
        }

        Void totally_reset(PartitionInfo<D> pi);
        Void add_to_partition(TreeCoords<D> p);

        Key<D> first_child(const Key<D>& key, const LBNode<D>& node);
        Key<D> next_sibling(const Key<D>& key);

        bool reset_partition(Cost& partition_size, Cost& used_up, PartitionInfo<D>& pi);


        MyPmap<D>& get_mypmap() {
            return *static_cast< MyPmap<D>* >(impl.get_pmap().get());
        }

        template <typename Archive>
        void serialize(const Archive& ar) {
            ar & root & list_of_list & cost_list & impl;
        }


    };

    template <int D>
    const Key<D> LBTree<D>::root(0);

    template <int D>
    typename DClass<D>::vvTreeCoords LBTree<D>::list_of_list;

    template <int D>
    typename std::vector<Cost> LBTree<D>::cost_list;

    /// Implementation of load balancing

    /// Implements the load balancing algorithm upon the tree underlying a function.
    template <int D>
    class LoadBalImpl {
    private:
        typedef MyPmap<D> Pmap;
        int k;
        double comm_bandw;
        double comm_latency;
        double flop_time;
        std::shared_ptr<LBTree<D> > skeltree;
        World& world;

        template<typename T, typename costfunT>
        void construct_skel(const std::shared_ptr<FunctionImpl<T,D> >& f, const costfunT& costfun) {
            skeltree.reset(new LBTree<D>(f->world, f->coeffs.get_pmap()));
//            madness::print("about to initialize tree");
            skeltree->template init_tree<T>(f,costfun);
//            madness::print("just initialized tree");
            pi.skel_cost = skeltree->fix_cost();
            //	    madness::print("construct_skel: pi.skel_cost =", pi.skel_cost);
            pi.cost_left = pi.skel_cost;
        }

    public:
        PartitionInfo<D> pi;
        //Constructors
        template <typename T, typename costfunT>
        LoadBalImpl(Function<T,D>& f, const costfunT& costfun,
                    double a=1e-8, double b=1e-5, double c=5e-10, double facter=1.1) : k(f.get_impl()->k)
                , comm_bandw(a)
                , comm_latency(b)
                , flop_time(c)
                , world(f.get_impl()->world)
                , pi(PartitionInfo<D>(facter)) {
            construct_skel(f.get_impl(), costfun);
            pi.partition_number = f.get_impl()->world.mpi.nproc()-1;
        }

        ~LoadBalImpl() {};

        //Methods

        /// Returns a shared pointer to a new process map, which can then be used to redistribute the function
        std::shared_ptr< WorldDCPmapInterface< Key<D> > > load_balance() {
            return std::shared_ptr< WorldDCPmapInterface< Key<D> > >(new MyPmap<D>(world, find_best_partition()));
        }

        std::vector< TreeCoords<D> > find_best_partition();
        std::vector< std::vector< TreeCoords<D> > > find_all_partitions();

        CompCost compute_comp_cost(Cost c, int n);

        template <typename T, typename costfunT>
        void add_tree(Function<T,D> f, const costfunT& costfun) {
            skeltree->add_tree(f.get_impl(), costfun);
            pi.skel_cost = skeltree->fix_cost();
            pi.cost_left = pi.skel_cost;

        }

    };


    Cost compute_partition_size(Cost cost, unsigned int parts);


}

#endif // MADNESS_MRA_LOADBAL_H__INCLUDED
