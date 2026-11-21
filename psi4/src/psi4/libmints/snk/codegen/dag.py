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

import math
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
        self.schedules = {}      # (roots, reorder) -> result of schedule()
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
        """Algebraic simplifications of the graph; returns the nodes that replace roots.

        1. v*x + s*(v*y) = v*(x + s*y) where the product v*y has no other user (v is a
           per-gridpoint factor, s a scalar): a VRR step and one of its differences carry
           the same factor PG_i, which happens for same center integrals.
        2. A product s*x with a scalar is merged into the steps that use it, unless it is
           the plain addend of at least two two-term steps (where it costs nothing extra).
        3. Integer multipliers are folded into the scalar products below them.
        4. s*x + s*y = s*(x + y) for terms of a step with the same coefficient.
        5. Intermediates with only negative terms are negated, and their users with them.

        The node numbers of before the call are no longer valid."""
        nuses = {}
        for r in roots:
            nuses[r] = nuses.get(r, 0) + 1
        order = sorted(self._subtree_all(roots))
        for n in order:
            for o in self.operands(n):
                nuses[o] = nuses.get(o, 0) + 1
        old_defs = self.defs
        self.defs, self.memo = [], {}
        new = {}
        for n in order:
            if old_defs[n][0] == 'leaf':
                new[n] = self.leaf(old_defs[n][1])
                continue
            terms = list(old_defs[n][1])
            for a, (ka, sa, va, na) in enumerate(terms):
                # a: the term s*(v*y), its operand is the product v*y
                if va is not None or old_defs[na][0] != 'lin' or nuses[na] != 1:
                    continue
                product = old_defs[na][1]
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
        roots = [new[r] for r in roots]
        roots = self._fold_and_inline_scalars(roots)
        roots = self._factor_same_coeff(roots)
        roots = self._normalize_signs(roots)
        self.schedules = {}
        return roots

    def _fold_and_inline_scalars(self, roots):
        for _ in range(2):
            order = sorted(self._subtree_all(roots))
            users = {n: [] for n in order}
            for r in roots:
                users.setdefault(r, []).append((None, 1, None, None))
            for n in order:
                if self.defs[n][0] == 'lin':
                    for (k, s, v, ch) in self.defs[n][1]:
                        users[ch].append((n, k, s, v))
            old_defs = self.defs
            self.defs, self.memo = [], {}
            new = {}
            for n in order:
                if old_defs[n][0] == 'leaf':
                    new[n] = self.leaf(old_defs[n][1])
                    continue
                terms = []
                for (k, s, v, ch_orig) in old_defs[n][1]:
                    ch = new[ch_orig]
                    if s is None and v is None and self.defs[ch][0] == 'lin' and len(self.defs[ch][1]) == 1:
                        (k0, s0, v0, x0) = self.defs[ch][1][0]
                        if s0 is not None and v0 is None:
                            two_term_unit_users = [
                                u for (u, ku, su, vu) in users.get(ch_orig, [])
                                if u is not None and abs(ku) == 1 and su is None and vu is None
                                and len(old_defs[u][1]) == 2
                                and any(st is not None or vt is not None for (_, st, vt, _) in old_defs[u][1])
                            ]
                            other_keep = [u for (u, ku, su, vu) in users.get(ch_orig, [])
                                          if u is None or not (su is None and vu is None)]
                            if len(two_term_unit_users) <= 1 and len(other_keep) == 0:
                                terms.append((k * k0, s0, None, x0))
                                continue
                    terms.append((k, s, v, ch))
                new[n] = self.lin(terms)
            roots = [new[r] for r in roots]

            order = sorted(self._subtree_all(roots))
            users = {n: [] for n in order}
            for n in order:
                if self.defs[n][0] == 'lin':
                    for (k, s, v, ch) in self.defs[n][1]:
                        users[ch].append((n, k, s, v))
            root_set = set(roots)
            can_scale = {}
            for n in order:
                if self.defs[n][0] == 'leaf':
                    can_scale[n] = False
                    continue
                can_scale[n] = all(s is not None or (can_scale.get(ch, False) and len(users[ch]) == 1)
                                   for (k, s, v, ch) in self.defs[n][1])
            scale = {}
            for n in reversed(order):
                if n in root_set or not can_scale.get(n, False) or not users[n]:
                    continue
                if all(s is None and v is None for (_, _, s, v) in users[n]):
                    c = 0
                    for (u, k, _, _) in users[n]:
                        c = math.gcd(c, abs(k))
                    if c > 1:
                        scale[n] = c
                elif len(users[n]) == 1 and users[n][0][2] is None:
                    c = abs(users[n][0][1])
                    if c > 1:
                        scale[n] = c
                if n in scale:
                    c = scale[n]
                    new_terms = [(k * c, s, v, ch) for (k, s, v, ch) in self.defs[n][1]]
                    self.defs[n] = ('lin', tuple(new_terms))
                    for (k, s, v, ch) in new_terms:
                        users[ch] = [(u, k if u == n and s_u == s and v_u == v else k_u, s_u, v_u)
                                     for (u, k_u, s_u, v_u) in users[ch]]
            old_defs = self.defs
            self.defs, self.memo = [], {}
            new = {}
            for n in order:
                if old_defs[n][0] == 'leaf':
                    new[n] = self.leaf(old_defs[n][1])
                    continue
                terms = [((k // (scale.get(ch, 1) if s is None else 1)), s, v, new[ch])
                         for (k, s, v, ch) in old_defs[n][1]]
                new[n] = self.lin(terms)
            roots = [new[r] for r in roots]
        return roots

    def _factor_same_coeff(self, roots):
        order = sorted(self._subtree_all(roots))
        old_defs = self.defs
        self.defs, self.memo = [], {}
        new = {}
        for n in order:
            if old_defs[n][0] == 'leaf':
                new[n] = self.leaf(old_defs[n][1])
                continue
            terms = [(k, s, v, new[ch]) for (k, s, v, ch) in old_defs[n][1]]
            groups = {}
            for i, (k, s, v, ch) in enumerate(terms):
                if s is not None or v is not None:
                    groups.setdefault((abs(k), s, v), []).append(i)
            used_idx = set()
            new_terms = []
            for (ak, s, v), idxs in groups.items():
                if len(idxs) >= 2:
                    inner = self.lin([(1 if terms[i][0] > 0 else -1, None, None, terms[i][3]) for i in idxs])
                    new_terms.append((ak, s, v, inner))
                    used_idx.update(idxs)
            if used_idx:
                for i, t in enumerate(terms):
                    if i not in used_idx:
                        new_terms.append(t)
                new[n] = self.lin(new_terms)
            else:
                new[n] = self.lin(terms)
        return [new[r] for r in roots]

    def _normalize_signs(self, roots):
        order = sorted(self._subtree_all(roots))
        root_set = set(roots)
        flip = set()
        old_defs = self.defs
        self.defs, self.memo = [], {}
        new = {}
        for n in order:
            if old_defs[n][0] == 'leaf':
                new[n] = self.leaf(old_defs[n][1])
                continue
            terms = [(-k if ch in flip else k, s, v, new[ch]) for (k, s, v, ch) in old_defs[n][1]]
            if n not in root_set and all(k < 0 for (k, _, _, _) in terms):
                terms = [(-k, s, v, ch) for (k, s, v, ch) in terms]
                flip.add(n)
            new[n] = self.lin(terms)
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
        alive across many statements costs a spill to the stack.  Several variants of this
        greedy order (which values count as live, how the next output is rated, which
        output comes first) are then tried and the one with the fewest live values is kept."""
        cached = self.schedules.get((tuple(roots), reorder))
        if cached is None:
            cached = self.schedules[(tuple(roots), reorder)] = self._schedule(roots, reorder)
        return cached

    def _schedule(self, roots, reorder):
        all_nodes = self._subtree_all(roots)
        nuses = {}
        users = {n: set() for n in all_nodes}
        for r in roots:
            nuses[r] = nuses.get(r, 0) + 1
        for n in all_nodes:
            for o in self.operands(n):
                nuses[o] = nuses.get(o, 0) + 1
                users[o].add(n)

        orig_order, orig_done = [], set()
        pending = list(dict.fromkeys(roots))
        while pending:
            if not reorder:
                best, bsub = pending[0], self._subtree(pending[0], orig_done)
            else:
                best, bkey, bsub = None, None, None
                for r in pending:
                    sub = self._subtree(r, orig_done)
                    inside = set(sub)
                    stay = sum(1 for n in sub if any(u not in inside for u in users.get(n, ())))
                    die = sum(1 for o in {o for n in sub for o in self.operands(n)}
                              if o in orig_done and all(u in inside or u in orig_done for u in users[o]))
                    key = (stay - die, len(sub))
                    if bkey is None or key < bkey:
                        best, bkey, bsub = r, key, sub
            pending.remove(best)
            orig_order.extend(bsub)
            orig_done.update(bsub)
        if not reorder or any(users[r] for r in roots):
            return orig_order, nuses

        weight = {}
        for n in all_nodes:
            if self.defs[n][0] == 'leaf':
                weight[n] = 0
            else:
                weight[n] = 1 + sum(weight[o] for o in set(self.operands(n)))

        non_leaves = [n for n in all_nodes if self.defs[n][0] != 'leaf']
        position = {n: i for i, n in enumerate(non_leaves)}
        unique_roots = list(dict.fromkeys(roots))
        root_set = set(roots)

        def liveness_score(order):
            last_use = {}
            for i, n in enumerate(order):
                for o in self.operands(n):
                    last_use[o] = max(last_use.get(o, 0), i)
                if n in root_set:
                    last_use[n] = max(last_use.get(n, 0), i)
            live_all, live_nl = set(), set()
            max_all = sum_all = max_nl = sum_nl = 0
            for i, n in enumerate(order):
                live_all.add(n)
                if self.defs[n][0] != 'leaf':
                    live_nl.add(n)
                for o in list(live_all):
                    if last_use.get(o, -1) <= i:
                        live_all.discard(o)
                        live_nl.discard(o)
                if self.defs[n][0] != 'leaf':
                    max_all = max(max_all, len(live_all))
                    sum_all += len(live_all)
                    max_nl = max(max_nl, len(live_nl))
                    sum_nl += len(live_nl)
            return (max_all + max_nl, sum_all + sum_nl)

        def run_schedule(inc_leaf, w_stay, w_prog, start_root=None):
            rem_users = {n: set(u) for n, u in users.items()}
            order, done, live = [], set(), set()
            pending = list(unique_roots)
            # intermediates whose operands have all been emitted
            waiting = {n: len(set(self.operands(n))) for n in non_leaves}
            ready = set()

            def mark_done(n):
                done.add(n)
                order.append(n)
                ready.discard(n)
                for u in users[n]:
                    waiting[u] -= 1
                    if waiting[u] == 0:
                        ready.add(u)

            def subtree_ordered(root):
                if root in done:
                    return []
                seen = {root}
                out = []
                def dfs(node):
                    ops = [ch for ch in self.operands(node) if ch not in done and ch not in seen]
                    if len(ops) > 1:
                        ops.sort(key=lambda ch: (
                            -sum(1 for o in self.operands(ch) if o in live and len(rem_users.get(o, ())) <= 2),
                            -sum(1 for o in self.operands(ch) if o in live),
                            -weight[ch],
                            ch
                        ))
                    for ch in ops:
                        if ch not in done and ch not in seen:
                            seen.add(ch)
                            dfs(ch)
                    out.append(node)
                dfs(root)
                return out

            def emit_node(n):
                if n in done:
                    return
                for o in self.operands(n):
                    if o not in done and self.defs[o][0] == 'leaf':
                        mark_done(o)
                        if inc_leaf and rem_users.get(o):
                            live.add(o)
                mark_done(n)
                if n in pending:
                    pending.remove(n)
                if (inc_leaf or self.defs[n][0] != 'leaf') and rem_users.get(n):
                    live.add(n)
                for o in self.operands(n):
                    if o in rem_users:
                        rem_users[o].discard(n)
                        if not rem_users[o]:
                            live.discard(o)

            def collapse_ready():
                changed = True
                while changed:
                    changed = False
                    best_n, best_d = None, 0
                    for n in sorted(ready, key=position.__getitem__):
                        stay = 1 if rem_users.get(n) else 0
                        die = sum(1 for o in set(self.operands(n)) if o in live and rem_users.get(o) == {n})
                        if stay - die < best_d:
                            best_n, best_d = n, stay - die
                    if best_n is not None:
                        emit_node(best_n)
                        changed = True

            first = start_root
            while pending:
                collapse_ready()
                if not pending:
                    break
                if first is not None and first in pending:
                    bsub = subtree_ordered(first)
                    first = None
                else:
                    bkey, best = None, None
                    for r in pending:
                        # the rating only needs the set of nodes below r, not their order
                        sub = self._subtree(r, done)
                        inside = set(sub)
                        stay = sum(1 for n in sub if (inc_leaf or self.defs[n][0] != 'leaf')
                                   and any(u not in inside and u not in done for u in users[n]))
                        touch_ops = {o for n in sub for o in self.operands(n) if o in live}
                        die = sum(1 for o in touch_ops if all(u in inside or u in done for u in users[o]))
                        progress = math.fsum(sum(1 for u in rem_users[o] if u in inside) / len(rem_users[o])
                                             for o in touch_ops if rem_users.get(o))
                        key = (w_stay * stay - die - w_prog * progress, -len(touch_ops), len(sub))
                        if bkey is None or key < bkey:
                            bkey, best = key, r
                    bsub = subtree_ordered(best)
                for n in bsub:
                    emit_node(n)
                    collapse_ready()
            return order

        best_order, best_score = orig_order, liveness_score(orig_order)
        step = max(1, len(unique_roots) // 6)
        starts = [None] + [unique_roots[i] for i in range(0, len(unique_roots), step)]
        for inc_leaf in (True, False):
            for (w_stay, w_prog) in ((1.0, 0.0), (1.0, 0.8), (1.5, 1.0)):
                for s in starts:
                    cand = run_schedule(inc_leaf, w_stay, w_prog, s)
                    score = liveness_score(cand)
                    if score < best_score:
                        best_score, best_order = score, cand
        return best_order, nuses

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
#   key = (ax, ay, az, bx, by, bz, N, d) is the intermediate [a|b]^(N)     #
#   for d = 0 and for d > 0 its d-th finite difference with respect to    #
#   the auxiliary index, weighted with (1/2p)^d,                          #
#                                                                         #
#       W^0 [a|b]^(N) = [a|b]^(N)                                         #
#       W^d [a|b]^(N) = (1/2p) (W^(d-1) [a|b]^(N) - W^(d-1) [a|b]^(N+1))  #
#                                                                         #
#   The VRR only ever needs lower intermediates of the same auxiliary     #
#   index through such differences,                                       #
#                                                                         #
#       [a+1_i|b]^(N) = PA_i [a|b]^(N) - PG_i [a|b]^(N+1)                 #
#                     + a_i W[a-1_i|b]^(N) + b_i W[a|b-1_i]^(N)           #
#                                                                         #
#   and because its coefficients do not depend on N, the differences      #
#   satisfy the same recurrence with d -> d+1 in the last two terms.      #
#   A difference can therefore be formed either by subtracting two        #
#   intermediates or by recursion; the path search picks the cheaper.     #
#                                                                         #
#   The factor 1/2p can be applied once per difference (which the steps   #
#   that use it then add) or in every step that uses the difference.      #
#   Which is cheaper depends on how many steps share a difference and     #
#   is settled per difference by Graph.factor(); the search has to        #
#   assume one of the two, selected by shared_differences.  In that       #
#   model, keys with d < 0 are the plain |d|-th differences of the Boys   #
#   functions, which are then weighted with a single product.             #
#                                                                         #
###########################################################################

XYZ = "xyz"

class VRRBuilder:
    def __init__(self, same_center, graph, shared_differences=False, differences=True, grouped=True, kicks=0):
        self.sc = same_center
        self.g = graph
        self.shared = shared_differences
        self.differences = differences
        self.grouped = grouped
        self.kicks = kicks
        self.search = PathSearch(self.rule)
        self.steps = {}     # key -> list of term lists (one per option), parallel to the rule
        self.tags = {}      # key -> list of step labels (one per option)
        self.node = {}
        self.unweighted = {}

        graph.define_scalar("tpi", "0.5 / v_p[prim]")
        for c in XYZ:
            graph.define_scalar(f"PA{c}", f"v_PA{c}[prim]")
            graph.define_scalar(f"PB{c}", f"v_PA{c}[prim] + AB{c}")

    def _tpi_pow(self, d):
        if d == 1:
            return "tpi"
        name = f"tpi_p{d}"
        prev = self._tpi_pow(d - 1)
        self.g.define_scalar(name, f"tpi*{prev}")
        return name

    def target(self, a, b):
        if self.sc:
            # same center: (x-A)^a (x-A)^b = (x-A)^(a+b)
            a, b = (a[0]+b[0], a[1]+b[1], a[2]+b[2]), (0, 0, 0)
        return tuple(a) + tuple(b) + (0, 0)

    def rule(self, key):
        idx, N, d = list(key[0:6]), key[6], key[7]
        L = sum(idx)
        if d < 0:
            ad = -d
            c0 = tuple(idx) + (N, 0 if ad == 1 else -(ad - 1))
            c1 = tuple(idx) + (N + 1, 0 if ad == 1 else -(ad - 1))
            terms = [(1, None, None, c0), (-1, None, None, c1)]
            self.steps[key] = [terms]
            self.tags[key] = [-1]
            return [(ADD_COST, (c0, c1))]
        if L == 0 and d == 0:
            return None
        opts, steps, tags = [], [], []
        if d > 0:
            if self.shared and L == 0:
                u_key = tuple(idx) + (N, -d)
                terms = [(1, self._tpi_pow(d), None, u_key)]
                opts.append((1.0, (u_key,)))
            else:
                c0 = tuple(idx) + (N, d-1)
                c1 = tuple(idx) + (N+1, d-1)
                terms = [(1, "tpi", None, c0), (-1, "tpi", None, c1)]
                opts.append(((1.0 + ADD_COST) if self.shared else ADD_COST, (c0, c1)))
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
                    terms, cost_terms = [], []
                    if not self.sc:
                        terms.append((1, f"P{'AB'[c]}{XYZ[i]}", None, tuple(low) + (N, d)))
                        cost_terms.append((1, f"P{'AB'[c]}{XYZ[i]}", None, tuple(low) + (N, d)))
                    terms.append((-1, None, f"PG{XYZ[i]}[g]", tuple(low) + (N+1, d)))
                    cost_terms.append((-1, None, f"PG{XYZ[i]}[g]", tuple(low) + (N+1, d)))
                    for cc in (0, 1):
                        q = 3*cc + i
                        if low[q] > 0:
                            low2 = list(low)
                            low2[q] -= 1
                            terms.append((low[q], None, None, tuple(low2) + (N, d+1)))
                            cost_terms.append((low[q], None if self.shared else "tpi", None, tuple(low2) + (N, d+1)))
                    opts.append((term_ops(cost_terms), tuple(t[3] for t in terms)))
                    steps.append(terms)
                    tags.append(p)
        self.steps[key] = steps
        self.tags[key] = tags
        return opts

    @staticmethod
    def order(key):
        # parents before children: decreasing total L, then decreasing difference order
        d = key[7]
        return (-sum(key[0:6]), -abs(d), 0 if d >= 0 else 1, key)

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

    def _can_unweight(self, key):
        d = key[7]
        if d <= 0:
            return True
        if self.search.options(key) is None:
            return True
        choice = self.search.chosen(key)
        if self.tags[key][choice] != -1:
            return False
        return all(self._can_unweight(t[3]) for t in self.steps[key][choice])

    def _get_unweighted(self, key):
        if key[7] <= 0:
            return self._materialize(key)
        u = self.unweighted.get(key)
        if u is not None:
            return u
        terms = self.steps[key][self.search.chosen(key)]
        if len(terms) == 1:
            u = self._materialize(terms[0][3])
        else:
            u = self.g.lin([(1, None, None, self._get_unweighted(terms[0][3])),
                            (-1, None, None, self._get_unweighted(terms[1][3]))])
        self.unweighted[key] = u
        return u

    def _materialize(self, key):
        n = self.node.get(key)
        if n is not None:
            return n
        if self.search.options(key) is None:
            n = self.g.leaf(f"scaled_boys_{key[6]}[g]")
        else:
            choice = self.search.chosen(key)
            terms = self.steps[key][choice]
            if key[7] > 0 and self.tags[key][choice] == -1:
                if self._can_unweight(key):
                    u = self._get_unweighted(key)
                    n = self.g.lin([(1, self._tpi_pow(key[7]), None, u)])
                else:
                    diff = self.g.lin([(1, None, None, self._materialize(terms[0][3])),
                                       (-1, None, None, self._materialize(terms[1][3]))])
                    n = self.g.lin([(1, "tpi", None, diff)])
            else:
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
    def __init__(self, same_center, graph, hrr_blist, kicks=0):
        self.sc = same_center
        self.g = graph
        self.blist = hrr_blist
        self.kicks = kicks
        self.search = PathSearch(self.rule)
        self.steps = {}
        self.tags = {}
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
        opts, steps, tags = [], [], []
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
            tags.append(src + i)
        self.steps[key] = steps
        self.tags[key] = tags
        return opts

    @staticmethod
    def order(key):
        LA, LB = key[6], key[7]
        moved = sum(key[3:6]) if LA >= LB else sum(key[0:3])
        return (-moved, key)

    def policy(self, largest):
        def candidates(key):
            tags = self.tags[key]
            return [min(range(len(tags)), key=lambda n: (-key[tags[n]] if largest else key[tags[n]], tags[n]))]
        return candidates

    def build(self, targets, verbose=False):
        policies = [None, self.policy(False), self.policy(True)]
        self.search.optimize(targets, self.order, tags=self.tags, policies=policies,
                             kicks=self.kicks, verbose=verbose)
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
