"""Functions to enumerate all perfect and maximum matchings in bipartite graph.

Implemented following the algorithms in the paper "Algorithms for Enumerating
All Perfect, Maximum and Maximal Matchings in Bipartite Graphs" by Takeaki Uno,
using numpy and networkx modules of python.

NOTICE: optimization needed.

Author: guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
Update time: 2017-06-29 14:41:56.

Updated Dec 2017 LAB for pep8, py3, more tests, starter_match, and simpler interface

"""
from __future__ import print_function

try:
    import networkx as nx
except ImportError:
    raise ImportError("""Python module networkx not found. Solve by installing it: `conda install networkx` or `pip install networkx`""")
from networkx import bipartite
import numpy as np


def _plotGraph(graph):
    """Plot graph using nodes as position number."""

    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111)

    pos = [(ii[1], ii[0]) for ii in graph.nodes()]
    pos_dict = dict(zip(graph.nodes(), pos))
    nx.draw(graph, pos=pos_dict, ax=ax, with_labels=True)
    plt.show(block=False)
    return


def _formDirected(g, match):
    """Form directed graph D from G and matching M.

    Parameters
    ----------
    g : 
        Undirected bipartite graph. Nodes are separated by their
        'bipartite' attribute.
    match : 
        List of edges forming a matching of `g`. 
    
    Returns
    -------
    networkx.DiGraph
	    Directed graph, with edges in `match` pointing from set-0
	    (bipartite attribute==0) to set-1 (bipartite attrbiute==1), and
	    the other edges in `g` but not in `match` pointing from set-1 to
	    set-0.

    """
    d = nx.DiGraph()

    for ee in g.edges():
        if ee in match or (ee[1], ee[0]) in match:
            if g.node[ee[0]]['bipartite'] == 0:
                d.add_edge(ee[0], ee[1])
            else:
                d.add_edge(ee[1], ee[0])
        else:
            if g.node[ee[0]]['bipartite'] == 0:
                d.add_edge(ee[1], ee[0])
            else:
                d.add_edge(ee[0], ee[1])

    return d


def _enumMaximumMatching(g, starter_match=None):
    """Find all maximum matchings in an undirected bipartite graph `g`.

    Parameters
    ----------
    g : networkx.Graph
        Undirected bipartite graph. Nodes are separated by their
        'bipartite' attribute.
    starter_match : dict, optional
        Single perfect match to inaugurate Uno's algorithm.

    Returns
    -------
    list
        Each is a list of edges forming a maximum matching of `g`. 

    Author
    ------
    guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
    Update time: 2017-05-21 20:04:51.

    """
    all_matches = []

    #----------------Find one matching M----------------
    if starter_match is None:
        match = bipartite.hopcroft_karp_matching(g)
    else:
        match = starter_match

    #---------------Re-orient match arcs---------------
    match2 = []
    for kk, vv in match.items():
        if g.node[kk]['bipartite'] == 0:
            match2.append((kk, vv))
    match = match2
    all_matches.append(match)

    #-----------------Enter recursion-----------------
    all_matches = _enumMaximumMatchingIter(g, match, all_matches, None)

    return all_matches


def _enumMaximumMatchingIter(g, match, all_matches, add_e=None):
    """Recurively search maximum matchings.

    Parameters
    ----------
    g : 
        Undirected bipartite graph. Nodes are separated by their
        'bipartite' attribute.
    match : 
        List of edges forming one maximum matching of `g`.
    all_matches : 
	    List, each is a list of edges forming a maximum matching of `g`.
	    Newly found matchings will be appended into this list.
    add_e : tuple, optional
        Edge used to form subproblems. If not `None`, will be added to each
        newly found matchings.

    Returns
    -------
    list
        Updated list of all maximum matchings.

    Author
    ------
    guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
    Update time: 2017-05-21 20:09:06.

    """
    #---------------Form directed graph D---------------
    d = _formDirected(g, match)

    #-----------------Find cycles in D-----------------
    cycles = list(nx.simple_cycles(d))

    if len(cycles) == 0:

        #---------If no cycle, find a feasible path---------
        all_uncovered = set(g.node).difference(set([ii[0] for ii in match]))
        all_uncovered = all_uncovered.difference(set([ii[1] for ii in match]))
        all_uncovered = list(all_uncovered)

        #--------------If no path, terminiate--------------
        if len(all_uncovered) == 0:
            return all_matches

        #----------Find a length 2 feasible path----------
        idx = 0
        uncovered = all_uncovered[idx]
        while True:

            if uncovered not in nx.isolates(g):
                paths = nx.single_source_shortest_path(d, uncovered, cutoff=2)
                len2paths = [vv for kk, vv in paths.items() if len(vv) == 3]

                if len(len2paths) > 0:
                    reversed = False
                    break

                #----------------Try reversed path----------------
                paths_rev = nx.single_source_shortest_path(d.reverse(), uncovered, cutoff=2)
                len2paths = [vv for kk, vv in paths_rev.items() if len(vv) == 3]

                if len(len2paths) > 0:
                    reversed = True
                    break

            idx += 1
            if idx > len(all_uncovered) - 1:
                return all_matches

            uncovered = all_uncovered[idx]

        #-------------Create a new matching M'-------------
        len2path = len2paths[0]
        if reversed:
            len2path = len2path[::-1]
        len2path = list(zip(len2path[:-1], len2path[1:]))

        new_match = []
        for ee in d.edges():
            if ee in len2path:
                if g.node[ee[1]]['bipartite'] == 0:
                    new_match.append((ee[1], ee[0]))
            else:
                if g.node[ee[0]]['bipartite'] == 0:
                    new_match.append(ee)

        if add_e is not None:
            for ii in add_e:
                new_match.append(ii)

        all_matches.append(new_match)

        #---------------------Select e---------------------
        e = set(len2path).difference(set(match))
        e = list(e)[0]

        #-----------------Form subproblems-----------------
        g_plus = g.copy()
        g_minus = g.copy()
        g_plus.remove_node(e[0])
        g_plus.remove_node(e[1])

        g_minus.remove_edge(e[0], e[1])

        add_e_new = [
            e,
        ]
        if add_e is not None:
            add_e_new.extend(add_e)

        all_matches = _enumMaximumMatchingIter(g_minus, match, all_matches, add_e)
        all_matches = _enumMaximumMatchingIter(g_plus, new_match, all_matches, add_e_new)

    else:
        #----------------Find a cycle in D----------------
        cycle = cycles[0]
        cycle.append(cycle[0])
        cycle = list(zip(cycle[:-1], cycle[1:]))

        #-------------Create a new matching M'-------------
        new_match = []
        for ee in d.edges():
            if ee in cycle:
                if g.node[ee[1]]['bipartite'] == 0:
                    new_match.append((ee[1], ee[0]))
            else:
                if g.node[ee[0]]['bipartite'] == 0:
                    new_match.append(ee)

        if add_e is not None:
            for ii in add_e:
                new_match.append(ii)

        all_matches.append(new_match)

        #-----------------Choose an edge E-----------------
        e = set(match).intersection(set(cycle))
        e = list(e)[0]

        #-----------------Form subproblems-----------------
        g_plus = g.copy()
        g_minus = g.copy()
        g_plus.remove_node(e[0])
        g_plus.remove_node(e[1])
        g_minus.remove_edge(e[0], e[1])

        add_e_new = [
            e,
        ]
        if add_e is not None:
            add_e_new.extend(add_e)

        all_matches = _enumMaximumMatchingIter(g_minus, new_match, all_matches, add_e)
        all_matches = _enumMaximumMatchingIter(g_plus, match, all_matches, add_e_new)

    return all_matches


def _enumMaximumMatching2(g):
    """Find all maximum matchings in an undirected bipartite graph `g`.
    Similar to _enumMaximumMatching but implemented using adjacency matrix
    of graph for slight speed boost.

    Parameters
    ----------
    g: 
        Undirected bipartite graph. Nodes are separated by their
        'bipartite' attribute.

    Returns
    -------
    list
        Each is a list of edges forming a maximum matching of `g`. 

    Author
    ------
    guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
    Update time: 2017-05-21 20:04:51.

    """
    from scipy import sparse

    s1 = set(n for n, d in g.nodes(data=True) if d['bipartite'] == 0)
    s2 = set(g) - s1
    n1 = len(s1)
    nodes = list(s1) + list(s2)

    adj = nx.adjacency_matrix(g, nodes).tolil()
    all_matches = []

    #----------------Find one matching----------------
    match = bipartite.hopcroft_karp_matching(g)

    matchadj = np.zeros(adj.shape).astype('int')
    for kk, vv in match.items():
        matchadj[nodes.index(kk), nodes.index(vv)] = 1
    matchadj = sparse.lil_matrix(matchadj)

    all_matches.append(matchadj)

    #-----------------Enter recursion-----------------
    all_matches = _enumMaximumMatchingIter2(adj, matchadj, all_matches, n1, None, True)

    #---------------Re-orient match arcs---------------
    all_matches2 = []
    for ii in all_matches:
        match_list = sparse.find(ii[:n1] == 1)
        m1 = [nodes[jj] for jj in match_list[0]]
        m2 = [nodes[jj] for jj in match_list[1]]
        match_list = zip(m1, m2)

        all_matches2.append(match_list)

    print('got all')
    return all_matches2


def _enumMaximumMatchingIter2(adj, matchadj, all_matches, n1, add_e=None, check_cycle=True):
    """Recurively search maximum matchings.
    Similar to _enumMaximumMatching but implemented using adjacency matrix
    of graph for a slight speed boost.

    Parameters
    ----------
#    g : 
#        Undirected bipartite graph. Nodes are separated by their
#        'bipartite' attribute.
#    match : 
#        List of edges forming one maximum matching of `g`.
#    all_matches : 
#	    List, each is a list of edges forming a maximum matching of `g`.
#	    Newly found matchings will be appended into this list.
    add_e : tuple, optional
        Edge used to form subproblems. If not `None`, will be added to each
        newly found matchings.

    Returns
    -------
    list
        Updated list of all maximum matchings.

    Author
    ------
    guangzhi XU (xugzhi1987@gmail.com; guangzhi.xu@outlook.com)
    Update time: 2017-05-21 20:09:06.

    """
    from scipy import sparse

    #-------------------Find cycles-------------------
    if check_cycle:
        d = matchadj.multiply(adj)
        d[n1:, :] = adj[n1:, :] - matchadj[n1:, :].multiply(adj[n1:, :])

        dg = nx.from_numpy_matrix(d.toarray(), create_using=nx.DiGraph())
        cycles = list(nx.simple_cycles(dg))
        if len(cycles) == 0:
            check_cycle = False
        else:
            check_cycle = True

    if check_cycle:
        cycle = cycles[0]
        cycle.append(cycle[0])
        cycle = zip(cycle[:-1], cycle[1:])

        #--------------Create a new matching--------------
        new_match = matchadj.copy()
        for ee in cycle:
            if matchadj[ee[0], ee[1]] == 1:
                new_match[ee[0], ee[1]] = 0
                new_match[ee[1], ee[0]] = 0
                e = ee
            else:
                new_match[ee[0], ee[1]] = 1
                new_match[ee[1], ee[0]] = 1

        if add_e is not None:
            for ii in add_e:
                new_match[ii[0], ii[1]] = 1

        all_matches.append(new_match)

        #-----------------Form subproblems-----------------
        g_plus = adj.copy()
        g_minus = adj.copy()
        g_plus[e[0], :] = 0
        g_plus[:, e[1]] = 0
        g_plus[:, e[0]] = 0
        g_plus[e[1], :] = 0
        g_minus[e[0], e[1]] = 0
        g_minus[e[1], e[0]] = 0

        add_e_new = [
            e,
        ]
        if add_e is not None:
            add_e_new.extend(add_e)

        all_matches = _enumMaximumMatchingIter2(g_minus, new_match, all_matches, n1, add_e, check_cycle)
        all_matches = _enumMaximumMatchingIter2(g_plus, matchadj, all_matches, n1, add_e_new, check_cycle)

    else:
        #---------------Find uncovered nodes---------------
        uncovered = np.where(np.sum(matchadj, axis=1) == 0)[0]

        if len(uncovered) == 0:
            return all_matches

        #---------------Find feasible paths---------------
        paths = []
        for ii in uncovered:
            aa = adj[ii, :].dot(matchadj)
            if aa.sum() == 0:
                continue
            paths.append((ii, int(sparse.find(aa == 1)[1][0])))
            if len(paths) > 0:
                break

        if len(paths) == 0:
            return all_matches

        #----------------------Find e----------------------
        feas1, feas2 = paths[0]
        e = (feas1, int(sparse.find(matchadj[:, feas2] == 1)[0]))

        #----------------Create a new match----------------
        new_match = matchadj.copy()
        new_match[feas2, :] = 0
        new_match[:, feas2] = 0
        new_match[feas1, e[1]] = 1
        new_match[e[1], feas1] = 1

        if add_e is not None:
            for ii in add_e:
                new_match[ii[0], ii[1]] = 1

        all_matches.append(new_match)

        #-----------------Form subproblems-----------------
        g_plus = adj.copy()
        g_minus = adj.copy()
        g_plus[e[0], :] = 0
        g_plus[:, e[1]] = 0
        g_plus[:, e[0]] = 0
        g_plus[e[1], :] = 0
        g_minus[e[0], e[1]] = 0
        g_minus[e[1], e[0]] = 0

        add_e_new = [
            e,
        ]
        if add_e is not None:
            add_e_new.extend(add_e)

        all_matches = _enumMaximumMatchingIter2(g_minus, matchadj, all_matches, n1, add_e, check_cycle)
        all_matches = _enumMaximumMatchingIter2(g_plus, new_match, all_matches, n1, add_e_new, check_cycle)

    if len(all_matches) % 1000 == 0:
        print('len', len(all_matches))

    print('another')
    return all_matches


def _findCycle(adj, n1):
    from scipy import sparse

    path = []
    visited = set()

    def visit(v):
        if v in visited:
            return False
        visited.add(v)
        path.append(v)
        neighbours = sparse.find(adj[v, :] == 1)[1]
        for nn in neighbours:
            if nn in path or visit(nn):
                return True
        path.remove(v)
        return False

    nodes = range(n1)
    result = any(visit(v) for v in nodes)
    return result, path


def uno(edges, match=None, verbose=1):
    """Perform Uno's algorithm to enumerate all equivalent matchings among
    the bipartite graph defined by `edges`. Optionally given a single
    known `match`.

    http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.107.8179&rep=rep1&type=pdf

    "Algorithms for Enumerating All Perfect, Maximum and Maximal Matchings
    in Bipartite Graphs" by Takeaki UNO

    """
    if match is None:
        p_match = None
    else:
        p_match = {(1, m[0]): (0, m[1]) for m in match}
        p_m_inv = {(0, m[1]): (1, m[0]) for m in match}
        p_match.update(p_m_inv)

    p_edges = [[(1, e[0]), (0, e[1])] for e in edges]
    if verbose >= 2:
        print('Edges:')
        for e in p_edges:
            print('\t', e)

    g = nx.Graph()
    for e in p_edges:
        g.add_node(e[0], bipartite=0)
        g.add_node(e[1], bipartite=1)
        if verbose >= 2:
            print('  Node:', e[0], e[1])
    g.add_edges_from(p_edges)

    all_matches = _enumMaximumMatching(g, starter_match=p_match)
    p_all_matches = [sorted([(pt[0][1], pt[1][1]) for pt in am]) for am in all_matches]

    return p_all_matches


def example4(alg=1):

    edges = [(0, 0),
            (0, 1),
            (1, 5),
            (1, 6),
            (2, 0),
            (2, 1),
            (3, 5),
            (3, 6),
            (4, 2),
            (4, 3),
            (5, 2),
            (5, 3),
            (6, 4),
            (6, 7),
            (7, 4),
            (7, 7)]
    match = [(0, 0), (2, 1), (4, 2), (5, 3), (6, 4), (3, 5), (1, 6), (7, 7)]

    ref = [[(0, 0), (2, 1), (4, 2), (5, 3), (6, 4), (3, 5), (1, 6), (7, 7)],  # ----
           [(0, 1), (2, 0), (4, 2), (5, 3), (6, 4), (3, 5), (1, 6), (7, 7)],  # *---

           [(0, 0), (2, 1), (4, 3), (5, 2), (6, 4), (3, 5), (1, 6), (7, 7)],  # -*--
           [(0, 1), (2, 0), (4, 3), (5, 2), (6, 4), (3, 5), (1, 6), (7, 7)],  # **--

           [(0, 0), (2, 1), (4, 2), (5, 3), (6, 7), (3, 5), (1, 6), (7, 4)],  # --*-
           [(0, 1), (2, 0), (4, 2), (5, 3), (6, 7), (3, 5), (1, 6), (7, 4)],  # *-*-

           [(0, 0), (2, 1), (4, 2), (5, 3), (6, 4), (3, 6), (1, 5), (7, 7)],  # ---*
           [(0, 1), (2, 0), (4, 2), (5, 3), (6, 4), (3, 6), (1, 5), (7, 7)],  # *--*

           [(0, 0), (2, 1), (4, 3), (5, 2), (6, 7), (3, 5), (1, 6), (7, 4)],  # -**-
           [(0, 1), (2, 0), (4, 3), (5, 2), (6, 7), (3, 5), (1, 6), (7, 4)],  # ***-

           [(0, 0), (2, 1), (4, 3), (5, 2), (6, 4), (3, 6), (1, 5), (7, 7)],  # -*-*
           [(0, 1), (2, 0), (4, 3), (5, 2), (6, 4), (3, 6), (1, 5), (7, 7)],  # **-*

           [(0, 0), (2, 1), (4, 3), (5, 2), (6, 7), (3, 6), (1, 5), (7, 4)],  # -***
           [(0, 1), (2, 0), (4, 3), (5, 2), (6, 7), (3, 6), (1, 5), (7, 4)],  # ****

           [(0, 0), (2, 1), (4, 2), (5, 3), (6, 7), (3, 6), (1, 5), (7, 4)],  # --**
           [(0, 1), (2, 0), (4, 2), (5, 3), (6, 7), (3, 6), (1, 5), (7, 4)],  # *-**
        ]
    ref = [sorted(r) for r in ref]

#cost:
# [[ 0.000  0.000  83.505  83.505  53.406  3.378  3.378  53.406]
# [ 3.398  3.398  53.169  53.169  29.828  0.000  0.000  29.828]
# [ 0.000  0.000  83.293  83.293  53.237  3.336  3.336  53.237]
# [ 3.359  3.359  53.323  53.323  29.944  0.000  0.000  29.944]
# [ 83.559  83.559  0.000  0.000  3.372  53.380  53.380  3.372]
# [ 83.297  83.297  0.000  0.000  3.320  53.171  53.171  3.320]
# [ 53.240  53.240  3.379  3.379  0.000  29.830  29.830  0.000]
# [ 53.468  53.468  3.322  3.322  0.000  30.001  30.001  0.000]]
#ptsCR [(0, 0), (2, 1), (4, 2), (5, 3), (6, 4), (3, 5), (1, 6), (7, 7)]

    ans = uno(edges, verbose=2)
    _check('Example 4a (internal match)', ans, ref, verbose=2)

    ans = uno(edges, verbose=2, match=match)
    _check('Example 4b (provided match)', ans, ref, verbose=2)


def example3(alg=1):

    match = [(1, 2), (3, 4), (5, 6), (7, 8)]
    edges = [(1, 2),
             (1, 4),
             (1, 6),
             (3, 4),
             (3, 6),
             (3, 8), 
             (5, 6),
             (5, 8),
             (5, 2),
             (7, 8),
             (7, 2),
             (7, 4)]

    ref = [ [(1, 2), (3, 6), (5, 8), (7, 4)],
            [(1, 2), (3, 4), (5, 6), (7, 8)],
            [(1, 2), (3, 8), (5, 6), (7, 4)],
            [(1, 4), (3, 6), (5, 2), (7, 8)],
            [(1, 4), (3, 6), (5, 8), (7, 2)],
            [(1, 4), (3, 8), (5, 6), (7, 2)],
            [(1, 6), (3, 4), (5, 2), (7, 8)],
            [(1, 6), (3, 4), (5, 8), (7, 2)],
            [(1, 6), (3, 8), (5, 2), (7, 4)]]

    ans = uno(edges, verbose=2)
    _check('Example 3a (internal match)', ans, ref)

    ans = uno(edges, verbose=2, match=match)
    _check('Example 3b (provided match)', ans, ref, verbose=2)


def _check(msg, ans, ref, verbose=1):

    tans = [tuple(qw) for qw in ans]
    tref = [tuple(qw) for qw in ref]
    extra_answers = set(tans).difference(set(tref))
    missd_answers = set(tref).difference(set(tans))
    if verbose >= 2:
        for a in tans:
            print('Computed:', a)
        for a in tref:
            print('Supplied:', a)

    try:
        assert (extra_answers == set())
        assert (missd_answers == set())
    except AssertionError as err:
        print(msg, 'failed:')
        if extra_answers != set():
            for a in extra_answers:
                print('Incomplete Ref:', a)
        if missd_answers != set():
            for a in missd_answers:
                print('Incomplete Soln:', a)
        raise err
    else:
        print(msg, 'passed')


def example2(alg=1):
    """https://mathematica.stackexchange.com/questions/77410/find-all-perfect-matchings-of-a-graph/82893#82893"""

    g = nx.Graph()
    edges = [[(1, 1), (0, 2)],
             [(1, 1), (0, 4)],
             [(1, 1), (0, 6)],
             [(1, 3), (0, 4)],
             [(1, 3), (0, 6)],
             [(1, 3), (0, 8)], 
             [(1, 5), (0, 6)],
             [(1, 5), (0, 8)],
             [(1, 5), (0, 2)],
             [(1, 7), (0, 8)],
             [(1, 7), (0, 2)],
             [(1, 7), (0, 4)]]

#1 <-> 2, 3 <-> 6, 4 <-> 7, 5 <-> 8
#1 <-> 2, 3 <-> 4, 5 <-> 6, 7 <-> 8
#1 <-> 2, 3 <-> 8, 4 <-> 7, 5 <-> 6
#1 <-> 4, 2 <-> 5, 3 <-> 6, 7 <-> 8
#1 <-> 4, 2 <-> 7, 3 <-> 6, 5 <-> 8
#1 <-> 4, 2 <-> 7, 3 <-> 8, 5 <-> 6
#1 <-> 6, 2 <-> 5, 3 <-> 4, 7 <-> 8
#1 <-> 6, 2 <-> 7, 3 <-> 4, 5 <-> 8
#1 <-> 6, 2 <-> 5, 3 <-> 8, 4 <-> 7

#Match2: [(1, 2), (3, 6), (5, 8), (7, 4)]
#Match2: [(1, 2), (3, 4), (5, 6), (7, 8)]
#Match2: [(1, 2), (3, 8), (5, 6), (7, 4)]
#Match2: [(1, 4), (3, 6), (5, 2), (7, 8)]
#Match2: [(1, 4), (3, 6), (5, 8), (7, 2)]
#Match2: [(1, 4), (3, 8), (5, 6), (7, 2)]
#Match2: [(1, 6), (3, 4), (5, 2), (7, 8)]
#Match2: [(1, 6), (3, 4), (5, 8), (7, 2)]
#Match2: [(1, 6), (3, 8), (5, 2), (7, 4)]

    for ii in edges:
        g.add_node(ii[0], bipartite=0)
        g.add_node(ii[1], bipartite=1)

    g.add_edges_from(edges)
    #plotGraph(g)

    if alg == 1:
        all_matches = _enumMaximumMatching(g)
    elif alg == 2:
        all_matches = _enumMaximumMatching2(g)

    ref = [ [(1, 2), (3, 6), (5, 8), (7, 4)],
            [(1, 2), (3, 4), (5, 6), (7, 8)],
            [(1, 2), (3, 8), (5, 6), (7, 4)],
            [(1, 4), (3, 6), (5, 2), (7, 8)],
            [(1, 4), (3, 6), (5, 8), (7, 2)],
            [(1, 4), (3, 8), (5, 6), (7, 2)],
            [(1, 6), (3, 4), (5, 2), (7, 8)],
            [(1, 6), (3, 4), (5, 8), (7, 2)],
            [(1, 6), (3, 8), (5, 2), (7, 4)]]

    for mm in all_matches:
        ans = sorted([(ii[0][1], ii[1][1]) for ii in mm])
        if ans in ref:
            ref.remove(ans)
        print('Match2:', ans)
        g_match = nx.Graph()
        for ii in mm:
            g_match.add_edge(ii[0], ii[1])
        #plotGraph(g_match)

    assert (ref == [])
    print('Example 2 passed')


def example1(alg=1):
    g=nx.Graph()
    edges=[
            [(1,0), (0,0)],
            [(1,0), (0,1)],
            [(1,0), (0,2)],
            [(1,1), (0,0)],
            [(1,2), (0,2)],
            #[(1,2), (0,5)],
            [(1,3), (0,2)],
            #[(1,3), (0,3)],
            [(1,4), (0,3)],
            [(1,4), (0,5)],
            [(1,5), (0,2)],
            [(1,5), (0,4)],
            #[(1,5), (0,6)],
            [(1,6), (0,1)],
            [(1,6), (0,4)],
            [(1,6), (0,6)]
            ]

    for ii in edges:
        g.add_node(ii[0], bipartite=0)
        g.add_node(ii[1], bipartite=1)
        print('  Node:', ii[0], ii[1])

    g.add_edges_from(edges)
    #plotGraph(g)

    if alg == 1:
        all_matches = _enumMaximumMatching(g)
    elif alg == 2:
        all_matches = _enumMaximumMatching2(g)

    for mm in sorted(all_matches):
        ans = [(ii[0][1], ii[1][1]) for ii in mm]
        #print('Match:', mm)
        print('Match2:', sorted(ans))
        g_match = nx.Graph()
        for ii in mm:
            g_match.add_edge(ii[0], ii[1])
        #plotGraph(g_match)


#-------------Main---------------------------------
if __name__ == '__main__':
    import time

    t0 = time.time()
    example1(alg=1)
    t1 = time.time()
    example2(alg=1)
    t2 = time.time()
    #example1(alg=2)
    t3 = time.time()
    example2(alg=2)
    t4 = time.time()
    example3(alg=1)
    t5 = time.time()
    example4(alg=1)
    t6 = time.time()
    print('ex1 alg1:', t1 - t0)
    print('ex2 alg1:', t2 - t1)
    #print('ex1 alg2:', t3 - t2)
    print('ex2 alg2:', t4 - t3)
    print('ex3 alg1:', t5 - t4)
    print('ex4 alg1:', t6 - t5)
