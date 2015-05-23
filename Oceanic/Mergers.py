#!/usr/bin/env python3

from copy import deepcopy
from collections import defaultdict
from math import sqrt
from operator import itemgetter

from . import Raft as r
from . import GenomeData as gd

class Merger():
    
    def order(self, start, end):
        if start < end:
            return start, end
        else:
            return end, start
    
    def do_join(self, a, b, bridge):
        a.replace(bridge)
        b.empty()

    def get_hit_scaffolds(self, a, b, hits):
        a_hit_scaffolds = [scaffold for scaffold in a.scaffolds if scaffold in hits]
        if not a_hit_scaffolds:
            return [], []

        b_hit_scaffolds = [b_scaffold for b_scaffold in b.scaffolds for a_scaffold in a_hit_scaffolds if b_scaffold in hits[a_scaffold]]
        if not b_hit_scaffolds:
            return [], []
        
        return a_hit_scaffolds, b_hit_scaffolds


class MarkerMerge(Merger):
    
    def __init__(self, this, other, options=None):
        self.joins = []

    def merge(self):
        pass

    def bridge(self, a, b):
        if a.scaffold != b.scaffold:
            return None
        bridge = []
        scaffold = a.genome.blocks[a.scaffold]
    
        start, target, a_forward, b_forward = self.orient(a,b)
        this_cm = scaffold[start].cm
        direction = 1 if start < target else -1
        
        cm_dir = 0

        while True:
            if direction == 1:
                start = scaffold[start].next_block
            else:
                start = scaffold[start].prev_block

            if start == 0:
                break
            
            block_i_chr = scaffold[start].chromosome
            block_i_cm = scaffold[start].cm
            if block_i_chr == '0' or block_i_cm == -1 and block_i_chr == a.chromosome.name:
                bridge.append((a.scaffold, start, direction))
                continue
        
            if block_i_chr != a.chromosome.name: # This check comes after the cM check because chr==0 is OK
                break
        
            merge = False
            if block_i_cm == this_cm:
                merge = True
        
            if block_i_cm == a.chromosome.markers[this_cm].next_cm and cm_dir in (0, 1):
                this_cm = a.chromosome.markers[this_cm].next_cm
                cm_dir = 1
                merge = True
        
            
            if block_i_cm == a.chromosome.markers[this_cm].prev_cm and cm_dir in (0, -1):
                this_cm = a.chromosome.markers[this_cm].prev_cm
                cm_dir = -1
                merge = True
                
            
            if merge:
                if start == target:
                    bridge = self.merge_logs(bridge, a.logs)
                    bridge = self.merge_logs(bridge, b.logs)
                    self.do_join(a, b, bridge)
                    return True
                else:
                    bridge.append((a.scaffold, start, direction))
                
            if bridge and bridge[-1][1] == start:
                continue
    
            break
    
        return None

    def orient(self, a, b):
        edge_blocks = [a.start, a.last, b.start, b.last]
        if a.start < b.start:
            a_far, a_near, b_near, b_far = sorted(edge_blocks)
        else:
            b_far, b_near, a_near, a_far = sorted(edge_blocks)
        
        a_forward = a.start <= a.last
        b_forward = b.start <= b.last
        return a.get_log_start(a_near), b.get_log_start(b_near), a_forward, b_forward
    
    def merge_logs(self, this, other):
        if not this:
            return deepcopy(other)
    
        new = deepcopy(this) + deepcopy(other)
        new.sort(key=lambda x: x[1]) # Sort by start position
    
        this_dir = this[0][2]
        new = [(scaffold, start, this_dir) for scaffold, start, direction in new]
        if this_dir == -1:
            new.reverse()
    
        return new

class RopeMerge(Merger):
    
    def __init__(self, this, other, options=None):
        self.chromosome = this.chromosome
        self.springlines = defaultdict(lambda:defaultdict(bool))

        self.this = other
        self.other = other
        print(self.this)
        if self.this != self.other:
            print(self.other)

        this.turn_hooks('last')
        other.turn_hooks('first')
        
        self.ropes = defaultdict(self.Rope)

        self.tie_ropes(this, other)

        self.untangle_ropes()

    class PBHit:
        class ScaffoldHit:
            def __init__(self, scaffold, start, end, hitlen, seqlen, pccov):
                self.scaffold = scaffold
                self.start = start
                self.end = end
                self.hitlen = hitlen
                self.seqlen = seqlen
                self.pccov = pccov

        def __init__(self, hit):
            rstart, rend, qstart, qend, rhitlen, qhitlen, self.pcid, rseqlen, qseqlen, \
            rpccov, qpccov, rscaffold, qscaffold, self.hittype = hit
        
            self.pb = self.ScaffoldHit(rscaffold, rstart, rend, rhitlen, rseqlen, rpccov)
            self.g  = self.ScaffoldHit(qscaffold, qstart, qend, qhitlen, qseqlen, qpccov)
    
        def __repr__(self):
            return '''{0:16s} ({1:7d} bp) {2:7d} - {3:7d} ({4:10s} {5:7d} - {6:7d} ({7:7d} bp); {8:6.2f}%)'''.format(
                self.g.scaffold, self.g.seqlen, self.g.start, self.g.end, self.pb.scaffold, self.pb.start, self.pb.end, self.pb.hitlen, self.pcid
            )

    class SpringLine():
        def __init__(self, this_knot, other_knot):
            self.this_knot = this_knot
            self.other_knot = other_knot

        def __repr__(self):
            return repr(self.this_knot) + '\t----\t' + repr(self.other_knot)

    class Rope():
        def __init__(self, hit):
            self.name = hit.pb.scaffold
            self.length = hit.pb.seqlen
            self.knots = []
            self.lookup = defaultdict(r.Knot)
        
        def __repr__(self):
            out = '{}\t{}\n'.format(self.name, self.length)
            for knot in self.knots:
                out += '{}\n'.format(repr(self.lookup[knot]))
            return out

        @property
        def matching_hooks(self):
            hooks = defaultdict(int)
            for knot in self.knots:
                hooks[self.lookup[knot].hook] += 1

            return list(hooks.keys())

        def tie(self, hit, pool_knot):
            rope_knot = r.Knot(hit.scaffold, hit.start, hit.end)
            self.knots.append(rope_knot)
            self.lookup[rope_knot] = pool_knot
            self.lookup[pool_knot] = rope_knot

        def untie(self):
            for knot in self.knots:
                self.lookup[knot].untie()

        def get_springlines(self, knot, raft, raft_end):
            springlines = [RopeMerge.SpringLine(knot, other_knot) for other_knot in raft.hook_knots(raft_end, self.name)]
            return springlines

    def tie_ropes(self, this, other):
        pool_scaffolds = this.scaffolds + other.scaffolds
        statement = 'select * from pacbio_overlaps where qscaffold in ({})'.format(','.join('?'*len(pool_scaffolds)))
    
        for args in self.chromosome.overlaps.execute(statement, pool_scaffolds):
            hit = self.PBHit(args)
            self.tie_rope(self.ropes, hit, this, other)

    def tie_rope(self, ropes, hit, this, other):

        if hit.g.scaffold in this.scaffolds:
            hit_pool = this
        elif hit.g.scaffold in other.scaffolds:
            hit_pool = other
        else:
            print("No scaffold match!")
            return

        if hit.pb.scaffold not in ropes:
            ropes[hit.pb.scaffold] = self.Rope(hit)

        pool_knot = hit_pool.tie(hit.g, ropes[hit.pb.scaffold])

        if pool_knot:
            ropes[hit.pb.scaffold].tie(hit.pb, pool_knot)

    def untangle_ropes(self):
        to_remove = []
        for name in self.ropes:
            if len(self.ropes[name].knots) <= 1:
                to_remove.append(name)
            elif len(self.ropes[name].matching_hooks) == 1:
                to_remove.append(name)

        for name in to_remove:
            self.ropes[name].untie()
            del self.ropes[name]

    def get_springline(self, a, b):
        a_springline = a.get_springline('last', b)
        b_springline = b.get_springline('first', a)
        if a_springline == b_springline and a_springline is not None:
            return a_springline
        else:
            return False

    def bridge(self, a, b):
        if a == b or (a in self.springlines and b in self.springlines[a]):
            return False
            
        print("Searching for bridge between A:")
        print(a)
        print("and B:")
        print(b)
        springline = self.get_springline(a, b)
        print(springline)
        self.springlines[a][b] = springline
        self.springlines[b][a] = springline
        return True


    def merge(self):
        for a in self.springlines:
            for b in self.springlines[a]:
                print(a)
                print(b)
                print(self.springlines[a][b])


class OKMerge(Merger):
    def __init__(self, this, other, options=None):
        self.this = this
        self.other = other
        self.genome = this.genome
        
    def bridge(self, a, b):
        if self.this.pooltype != 'ok' or self.other.pooltype != 'ok':
            return False

        i = j = None
        for n, p in enumerate(self.this.chromosome.pools):
            if self.this is p:
                i = n
            elif self.other is p:
                j = n

        for n in range(i+1, j):
            if len(self.this.chromosome.pools[n]) != 0 or self.this.chromosome.pools[n].pooltype in ['other', 'orient']:
                return False
        
        newgap = self.this.genome.add_gap()
        a.append(newgap.scaffold, newgap.start, 1)
        
        a.merge(b)
        b.empty()

        return True
        
    def merge(self):
        pass


if __name__ == '__main__':
    pass