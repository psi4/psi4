#
#   Direct construction of the recurrence DAG (no symbolic expansion, no CSE)
#
#   The Obara-Saika VRR and the Head-Gordon-Pople HRR are evaluated as a directed
#   acyclic graph of intermediates.  Every intermediate can usually be reached by
#   several recurrence steps (which cartesian direction to reduce, and for the VRR
#   which center); PathSearch picks the steps so that the total number of floating
#   point operations of the whole kernel is as small as it can make it.  Graph then
#   holds the chosen intermediates as hash-consed linear combinations and schedules
#   them for emission as straight-line C++.
#

import random
import re

###########################################################################
#                                                                         #
#   Recurrence path search                                                #
#                                                                         #
#   A "rule" maps a key to the list of ways it can be computed:           #
#       rule(key) -> None                    (key is an input / leaf)     #
#       rule(key) -> [(cost, children), ...] (one entry per option)       #
#   The search minimises the summed cost of all keys that are reachable   #
#   from the targets through the chosen options.                          #
#                                                                         #
###########################################################################

class PathSearch:
    def __init__(self, rule):
        self.rule = rule
        self.opts = {}      # key -> list of (cost, children) or None
        self.choice = {}    # key -> index of the chosen option
        self.ref = {}       # key -> number of users among the active keys
        self.policy = None  # key -> option indices to consider when a key is first needed

    def options(self, key):
        o = self.opts.get(key, False)
        if o is False:
            o = self.rule(key)
            self.opts[key] = o
        return o

    def _default(self, key, opts):
        # first choice for a key: the candidate option that adds the least to the total cost
        cand = range(len(opts))
        if self.policy is not None:
            cand = self.policy(key)
        if len(cand) == 1:
            return cand[0]
        best, bestc = cand[0], None
        for n in cand:
            cost, children = opts[n]
            for ch in children:
                cost += self.acquire(ch)
            for ch in children:
                self.release(ch)
            if bestc is None or cost < bestc:
                best, bestc = n, cost
        return best

    def acquire(self, key):
        """Add one user of key; returns the cost of everything that became active."""
        total = 0
        stack = [key]
        while stack:
            k = stack.pop()
            opts = self.options(k)
            if opts is None:
                continue
            r = self.ref.get(k, 0)
            self.ref[k] = r + 1
            if r > 0:
                continue
            if k not in self.choice:
                self.choice[k] = self._default(k, opts)
            cost, children = opts[self.choice[k]]
            total += cost
            stack.extend(children)
        return total

    def release(self, key):
        """Remove one user of key; returns the cost of everything that became inactive."""
        total = 0
        stack = [key]
        while stack:
            k = stack.pop()
            opts = self.options(k)
            if opts is None:
                continue
            r = self.ref[k] - 1
            self.ref[k] = r
            if r > 0:
                continue
            cost, children = opts[self.choice[k]]
            total += cost
            stack.extend(children)
        return total

    def _switch(self, key, new):
        """Change the option of an active key; returns the change of the total cost."""
        opts = self.opts[key]
        old = self.choice[key]
        delta = opts[new][0] - opts[old][0]
        for ch in opts[new][1]:
            delta += self.acquire(ch)
        for ch in opts[old][1]:
            delta -= self.release(ch)
        self.choice[key] = new
        return delta

    def _group_moves(self, group, tags):
        """Try to give all active keys of a group the same kind of step; returns the gain."""
        members = {}
        for k, r in self.ref.items():
            if r > 0 and len(self.opts[k]) > 1:
                members.setdefault(group(k), []).append(k)
        total = 0
        for grp in sorted(members):
            keys = members[grp]
            if len(keys) < 2:
                continue
            for tag in sorted({t for k in keys for t in tags[k]}):
                undo, delta = [], 0
                for k in keys:
                    if self.ref.get(k, 0) == 0 or tag not in tags[k]:
                        continue
                    n = tags[k].index(tag)
                    if n != self.choice[k]:
                        undo.append((k, self.choice[k]))
                        delta += self._switch(k, n)
                if delta < 0:
                    total += delta
                else:
                    for (k, n) in reversed(undo):
                        if self.ref.get(k, 0) > 0:
                            self._switch(k, n)
                        else:
                            self.choice[k] = n
        return total

    def _local_search(self, order, group, tags, max_sweeps=50):
        """Improve the choices until no single change (or change of a whole group) helps."""
        total = 0
        for sweep in range(max_sweeps):
            improved = False
            if group is not None:
                gain = self._group_moves(group, tags)
                total += gain
                improved = gain < 0
            for key in sorted((k for k, r in self.ref.items() if r > 0), key=order):
                if self.ref.get(key, 0) == 0:
                    continue
                opts = self.opts[key]
                if len(opts) < 2:
                    continue
                cur = self.choice[key]
                best, bestd = cur, 0
                for n in range(len(opts)):
                    if n == cur:
                        continue
                    d = self._switch(key, n)
                    self._switch(key, cur)
                    if d < bestd:
                        best, bestd = n, d
                if best != cur:
                    total += self._switch(key, best)
                    improved = True
            if not improved:
                break
        return total

    def optimize(self, targets, order, group=None, tags=None, policies=(None,), kicks=0, seed=0, verbose=False):
        """Multi-start greedy construction and local search, then iterated local search.

        order(key) must sort parents before children (it is used to visit the keys top-down).
        group(key) and tags[key] (one tag per option) enable moves that switch all keys of
        a group to options with the same tag at once.  Each policy restricts the options
        that are considered when a key is first needed, which gives structurally different
        starting points.  The search is deterministic."""
        targets = sorted(set(targets), key=order)
        best = None
        for npol, policy in enumerate(policies):
            self.choice, self.ref, self.policy = {}, {}, policy
            total = 0
            for t in targets:
                total += self.acquire(t)
            start = total
            total += self._local_search(order, None, tags)
            if group is not None:
                # the two neighbourhoods end in different local minima; keep the better
                plain = (total, dict(self.choice), dict(self.ref))
                total += self._local_search(order, group, tags)
                if plain[0] < total:
                    total, self.choice, self.ref = plain
            if verbose:
                print(f"  path search, start {npol}: cost {start:g} -> {total:g}")
            if best is None or total < best[0]:
                best = (total, dict(self.choice), dict(self.ref))
        total, self.choice, self.ref = best
        self.policy = None
        rng = random.Random(seed)
        for kick in range(kicks):
            saved = (total, dict(self.choice), dict(self.ref))
            movable = sorted((k for k, r in self.ref.items() if r > 0 and len(self.opts[k]) > 1), key=order)
            if not movable:
                break
            for key in rng.sample(movable, min(len(movable), 1 + len(movable)//50)):
                if self.ref.get(key, 0) > 0:
                    n = rng.randrange(len(self.opts[key]) - 1)
                    total += self._switch(key, n if n < self.choice[key] else n + 1)
            total += self._local_search(order, group if (kick % 2) else None, tags)
            if total >= saved[0]:
                total, self.choice, self.ref = saved
            elif verbose:
                print(f"  path search: cost {total:g} after kick {kick}")
        return total

    def chosen(self, key):
        return self.choice[key]

###########################################################################
#                                                                         #
#   Expression graph                                                      #
#                                                                         #
#   A node is either a leaf (a C++ expression that is read from memory)   #
#   or a linear combination  sum_i k_i * s_i * v_i * node_i  where k is   #
#   an integer, s an optional loop-invariant scalar and v an optional     #
#   per-gridpoint array.                                                  #
#                                                                         #
###########################################################################

# Relative cost of an addition: multiplications and fused multiply-adds compete for the
# same (two) execution ports of current x86 cores, additions can also run elsewhere.
ADD_COST = 0.25

def term_counts(terms):
    """(multiplications, additions) of a linear combination when multiply-adds are fused."""
    unit = sum(1 for (k, s, v, _) in terms if s is None and v is None and abs(k) == 1)
    mult = len(terms) - unit
    return mult, max(0, unit - 1)

def term_ops(terms):
    mult, add = term_counts(terms)
    return mult + ADD_COST*add

class Graph:
    def __init__(self):
        self.defs = []           # id -> ('leaf', cstr) or ('lin', terms)
        self.memo = {}
        self.scalar_exprs = {}   # name -> C++ expression of a loop-invariant scalar
        self.used = set()        # names of the scalars that were emitted

    def define_scalar(self, name, expr):
        self.scalar_exprs.setdefault(name, expr)

    def _intern(self, d):
        i = self.memo.get(d)
        if i is None:
            i = len(self.defs)
            self.defs.append(d)
            self.memo[d] = i
        return i

    def leaf(self, cstr):
        return self._intern(('leaf', cstr))

    def lin(self, terms):
        merged = {}
        for (k, s, v, n) in terms:
            if k != 0:
                merged[(s, v, n)] = merged.get((s, v, n), 0) + k
        terms = tuple(sorted(((k, s, v, n) for (s, v, n), k in merged.items() if k != 0),
                             key=lambda t: (t[3], t[1] or '', t[2] or '', t[0])))
        if len(terms) == 0:
            raise RuntimeError("dag.Graph.lin(): empty linear combination.")
        if len(terms) == 1 and terms[0][:3] == (1, None, None):
            return terms[0][3]
        return self._intern(('lin', terms))

    def factor(self, roots):
        """Use  v*x + s*(v*y) = v*(x + s*y)  where the product v*y has no other user.

        v is a per-gridpoint factor and s a scalar.  This saves a multiplication whenever a
        VRR step and one of its differences carry the same factor PG_i, which happens for
        same center integrals.  Returns the nodes that replace roots."""
        nuses = {}
        for r in roots:
            nuses[r] = nuses.get(r, 0) + 1
        order = self._subtree_all(roots)
        for n in order:
            for o in self.operands(n):
                nuses[o] = nuses.get(o, 0) + 1
        new = {}
        for n in order:
            if self.defs[n][0] == 'leaf':
                new[n] = n
                continue
            terms = list(self.defs[n][1])
            for a, (ka, sa, va, na) in enumerate(terms):
                # a: the term s*(v*y), its operand is the product v*y
                if va is not None or self.defs[na][0] != 'lin' or nuses[na] != 1:
                    continue
                product = self.defs[na][1]
                if len(product) != 1 or product[0][1] is not None or product[0][2] is None:
                    continue
                (kp, _, v, y) = product[0]
                # b: the term v*x
                b = next((i for i, t in enumerate(terms) if i != a and t[1] is None and t[2] == v), None)
                if b is None:
                    continue
                (kb, _, _, x) = terms[b]
                inner = self.lin([(kb, None, None, new[x]), (ka*kp, sa, None, new[y])])
                terms = [t for i, t in enumerate(terms) if i not in (a, b)]
                terms = [(k, s, w, new[o]) for (k, s, w, o) in terms] + [(1, None, v, inner)]
                new[n] = self.lin(terms)
                break
            else:
                new[n] = self.lin([(k, s, w, new[o]) for (k, s, w, o) in terms])
        return [new[r] for r in roots]

    def operands(self, i):
        d = self.defs[i]
        return [t[3] for t in d[1]] if d[0] == 'lin' else []

    def count(self, roots):
        """(multiplications, additions) per gridpoint needed for the given roots."""
        mult = add = 0
        for i in self._subtree_all(roots):
            if self.defs[i][0] == 'lin':
                m, a = term_counts(self.defs[i][1])
                mult, add = mult + m, add + a
        return mult, add

    #
    #   Emission
    #

    def _scalar(self, k, s):
        if k != 1:
            # integer multiples of a scalar are hoisted out of the gridpoint loop
            self.define_scalar(f"{s}_{k}", f"{k}*{s}")
            s = f"{s}_{k}"
        self.used.add(s)
        return s

    def scalar_definitions(self):
        """Definitions of the used scalars (and of what they depend on), in dependency order."""
        lines, done = [], set()
        def visit(name):
            if name in done or name not in self.scalar_exprs:
                return
            done.add(name)
            for dep in re.findall(r"[A-Za-z_]\w*", self.scalar_exprs[name]):
                visit(dep)
            lines.append(f"const T {name} = {self.scalar_exprs[name]};")
        for name in sorted(self.used):
            visit(name)
        return lines

    def expr(self, i, names):
        """C++ expression for node i; names maps already emitted node ids to variable names."""
        d = self.defs[i]
        if d[0] == 'leaf':
            return d[1]
        pos, neg = [], []
        for (k, s, v, n) in d[1]:
            f = []
            if s is not None:
                f.append(self._scalar(abs(k), s))
            elif abs(k) != 1:
                f.append(str(abs(k)))
            if v is not None:
                f.append(v)
            f.append(names[n])
            (pos if k > 0 else neg).append("*".join(f))
        if pos:
            return " + ".join(pos) + "".join(" - " + t for t in neg)
        return "-" + neg[0] + "".join(" - " + t for t in neg[1:])

    def _subtree(self, root, done):
        """Nodes below root (inclusive) that are not in done, in depth-first postorder."""
        out = []
        if root in done:
            return out
        seen = {root}
        stack = [(root, iter(self.operands(root)))]
        while stack:
            node, it = stack[-1]
            for ch in it:
                if ch not in seen and ch not in done:
                    seen.add(ch)
                    stack.append((ch, iter(self.operands(ch))))
                    break
            else:
                stack.pop()
                out.append(node)
        return out

    def schedule(self, roots, reorder=True):
        """Order the nodes needed for roots; returns (order, nuses).

        An output is finished before the next one is started (depth-first), so that an
        intermediate is emitted close to its first user.  With reorder, the next output is
        the one that changes the number of live values the least: the kernels have many
        more live intermediates than there are vector registers, and every value that is
        alive across many statements costs a spill to the stack."""
        nuses = {}
        users = {}
        for r in roots:
            nuses[r] = nuses.get(r, 0) + 1
        for n in self._subtree_all(roots):
            for o in self.operands(n):
                nuses[o] = nuses.get(o, 0) + 1
                users.setdefault(o, set()).add(n)

        order, done = [], set()
        pending = list(dict.fromkeys(roots))
        while pending:
            if not reorder:
                best, bsub = pending[0], self._subtree(pending[0], done)
            else:
                best, bkey, bsub = None, None, None
                for r in pending:
                    sub = self._subtree(r, done)
                    inside = set(sub)
                    # new values that stay alive after this output, minus values that die with it
                    stay = sum(1 for n in sub if any(u not in inside for u in users.get(n, ())))
                    die = sum(1 for o in {o for n in sub for o in self.operands(n)}
                              if o in done and all(u in inside or u in done for u in users[o]))
                    key = (stay - die, len(sub))
                    if bkey is None or key < bkey:
                        best, bkey, bsub = r, key, sub
            pending.remove(best)
            order.extend(bsub)
            done.update(bsub)
        return order, nuses

    def _subtree_all(self, roots):
        done, out = set(), []
        for r in roots:
            sub = self._subtree(r, done)
            out.extend(sub)
            done.update(sub)
        return out

###########################################################################
#                                                                         #
#   Obara-Saika VRR                                                       #
#                                                                         #
#   key = (ax, ay, az, bx, by, bz, N, d) is the d-th finite difference    #
#   of the intermediate [a|b]^(N) with respect to the auxiliary index,    #
#                                                                         #
#       D^0 [a|b]^(N) = [a|b]^(N)                                         #
#       D^d [a|b]^(N) = D^(d-1) [a|b]^(N) - D^(d-1) [a|b]^(N+1)           #
#                                                                         #
#   The VRR only ever needs lower intermediates of the same auxiliary     #
#   index through such differences,                                       #
#                                                                         #
#       [a+1_i|b]^(N) = PA_i [a|b]^(N) - PG_i [a|b]^(N+1)                 #
#                     + (1/2p) (a_i D[a-1_i|b]^(N) + b_i D[a|b-1_i]^(N))  #
#                                                                         #
#   and because its coefficients do not depend on N, the differences      #
#   satisfy the same recurrence with d -> d+1 in the last two terms.      #
#   A difference can therefore be formed either by subtracting two        #
#   intermediates or by recursion; the path search picks the cheaper.     #
#                                                                         #
###########################################################################

XYZ = "xyz"

class VRRBuilder:
    def __init__(self, same_center, graph, differences=True, grouped=True, kicks=0):
        self.sc = same_center
        self.g = graph
        self.differences = differences
        self.grouped = grouped
        self.kicks = kicks
        self.search = PathSearch(self.rule)
        self.steps = {}     # key -> list of term lists (one per option), parallel to the rule
        self.tags = {}      # key -> list of step labels (one per option)
        self.node = {}

        graph.define_scalar("tpi", "0.5 / v_p[prim]")
        for c in XYZ:
            graph.define_scalar(f"PA{c}", f"v_PA{c}[prim]")
            graph.define_scalar(f"PB{c}", f"v_PA{c}[prim] + AB{c}")

    def target(self, a, b):
        if self.sc:
            # same center: (x-A)^a (x-A)^b = (x-A)^(a+b)
            a, b = (a[0]+b[0], a[1]+b[1], a[2]+b[2]), (0, 0, 0)
        return tuple(a) + tuple(b) + (0, 0)

    def rule(self, key):
        idx, N, d = list(key[0:6]), key[6], key[7]
        L = sum(idx)
        if L == 0 and d == 0:
            return None
        opts, steps, tags = [], [], []
        if d > 0:
            terms = [(1, None, None, tuple(idx) + (N, d-1)), (-1, None, None, tuple(idx) + (N+1, d-1))]
            opts.append((term_ops(terms), tuple(t[3] for t in terms)))
            steps.append(terms)
            tags.append(-1)
        if L > 0 and (d == 0 or self.differences):
            for c in (0, 1):
                for i in range(3):
                    p = 3*c + i
                    if idx[p] == 0:
                        continue
                    low = list(idx)
                    low[p] -= 1
                    terms = []
                    if not self.sc:
                        terms.append((1, f"P{'AB'[c]}{XYZ[i]}", None, tuple(low) + (N, d)))
                    terms.append((-1, None, f"PG{XYZ[i]}[g]", tuple(low) + (N+1, d)))
                    for cc in (0, 1):
                        q = 3*cc + i
                        if low[q] > 0:
                            low2 = list(low)
                            low2[q] -= 1
                            terms.append((low[q], "tpi", None, tuple(low2) + (N, d+1)))
                    opts.append((term_ops(terms), tuple(t[3] for t in terms)))
                    steps.append(terms)
                    tags.append(p)
        self.steps[key] = steps
        self.tags[key] = tags
        return opts

    @staticmethod
    def order(key):
        # parents before children: decreasing total L, then decreasing difference order
        return (-sum(key[0:6]), -key[7], key)

    @staticmethod
    def group(key):
        # the same intermediate for all auxiliary indices
        return key[0:6] + (key[7],)

    def policy(self, center, largest):
        """Starting point: reduce the given center first, along the direction with the
        largest (or smallest) angular momentum."""
        def candidates(key):
            tags = self.tags[key]
            cand = [n for n, t in enumerate(tags) if t < 0]
            rec = [n for n, t in enumerate(tags) if t >= 0]
            if rec:
                first = [n for n in rec if tags[n]//3 == center] or rec
                cand.append(min(first, key=lambda n: (-key[tags[n]] if largest else key[tags[n]], tags[n])))
            return cand
        return candidates

    def build(self, targets, verbose=False):
        """Returns the graph node ids of the target keys."""
        policies = [None] + [self.policy(c, big) for c in (1, 0) for big in (False, True)]
        self.search.optimize(targets, self.order, self.group if self.grouped else None, self.tags,
                             policies=policies, kicks=self.kicks, verbose=verbose)
        return [self._materialize(t) for t in targets]

    def _materialize(self, key):
        n = self.node.get(key)
        if n is not None:
            return n
        if self.search.options(key) is None:
            n = self.g.leaf(f"scaled_boys_{key[6]}[g]")
        else:
            terms = self.steps[key][self.search.chosen(key)]
            n = self.g.lin([(k, s, v, self._materialize(ch)) for (k, s, v, ch) in terms])
        self.node[key] = n
        return n

###########################################################################
#                                                                         #
#   Head-Gordon-Pople HRR                                                 #
#                                                                         #
#   key = (ax, ay, az, bx, by, bz, LA, LB); the base integrals are        #
#   (e,0,LA,LB) for LA >= LB ("left-to-right") and (0,e,LA,LB) otherwise  #
#                                                                         #
###########################################################################

class HRRBuilder:
    def __init__(self, same_center, graph, hrr_blist):
        self.sc = same_center
        self.g = graph
        self.blist = hrr_blist
        self.search = PathSearch(self.rule)
        self.steps = {}
        self.node = {}
        self.base = {}      # graph node of a base integral -> its index

    def base_index(self, node):
        return self.base[node]

    def rule(self, key):
        idx, LA, LB = list(key[0:6]), key[6], key[7]
        ltr = (LA >= LB)
        src, dst = (3, 0) if ltr else (0, 3)   # move quanta from src center to dst center
        if not any(idx[src:src+3]):
            return None
        opts, steps = [], []
        for i in range(3):
            if idx[src+i] == 0:
                continue
            low = list(idx)
            low[src+i] -= 1
            up = list(low)
            up[dst+i] += 1
            terms = [(1, None, None, tuple(up) + (LA, LB))]
            if not self.sc:
                terms.append((1 if ltr else -1, f"AB{XYZ[i]}", None, tuple(low) + (LA, LB)))
            opts.append((term_ops(terms), tuple(t[3] for t in terms)))
            steps.append(terms)
        self.steps[key] = steps
        return opts

    @staticmethod
    def order(key):
        LA, LB = key[6], key[7]
        moved = sum(key[3:6]) if LA >= LB else sum(key[0:3])
        return (-moved, key)

    def build(self, targets, verbose=False):
        self.search.optimize(targets, self.order, verbose=verbose)
        return [self._materialize(t) for t in targets]

    def _materialize(self, key):
        n = self.node.get(key)
        if n is not None:
            return n
        if self.search.options(key) is None:
            n = self.g.leaf(f"base_integral_{self.blist[key]}[g]")
            self.base[n] = self.blist[key]
        else:
            terms = self.steps[key][self.search.chosen(key)]
            n = self.g.lin([(k, s, v, self._materialize(ch)) for (k, s, v, ch) in terms])
        self.node[key] = n
        return n
