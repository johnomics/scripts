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
        
            if block_i_cm == a.chromosome.markers[this_cm].next_cm:
                this_cm = a.chromosome.markers[this_cm].next_cm
                merge = True
        
            
            if block_i_cm == a.chromosome.markers[this_cm].prev_cm:
                this_cm = a.chromosome.markers[this_cm].prev_cm
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
#        if springline.contains:
#            springline.contains.raft.discard(springline.contains.hook)
#        else:
        self.springlines[a][b] = springline
        self.springlines[b][a] = springline
        return True


    def merge(self):
        for a in self.springlines:
            for b in self.springlines[a]:
                print(a)
                print(b)
                print(self.springlines[a][b])



class OverlapMerge(Merger):
    
    def __init__(self, this, other, options=None):
        
        self.hits = defaultdict(lambda:defaultdict(list))
        statement = '''select * from genome_overlaps
                       where rscaffold in ({})
                       and   qscaffold in ({})
                    '''.format(','.join('?'*len(this.scaffolds)), ','.join('?'*len(other.scaffolds)))

        for *args, rscaffold, qscaffold, hittype in this.chromosome.overlaps.execute(statement, this.scaffolds + other.scaffolds):
            self.hits[rscaffold][qscaffold].append([rscaffold, qscaffold, hittype] + args)

    class Hit:
        class ScaffoldHit:
            def __init__(self, scaffold, start, end, hitlen, seqlen, pccov):
                self.scaffold = scaffold
                self.start = start
                self.end = end
                self.hitlen = hitlen
                self.seqlen = seqlen
                self.pccov = pccov
                
        def __init__(self, hit):
            rscaffold, qscaffold, self.hittype, rstart, rend, qstart, qend, \
            rhitlen, qhitlen, self.pcid, rseqlen, qseqlen, rpccov, qpccov = hit
            
            self.r = self.ScaffoldHit(rscaffold, rstart, rend, rhitlen, rseqlen, rpccov)
            self.q = self.ScaffoldHit(qscaffold, qstart, qend, qhitlen, qseqlen, qpccov)

    def merge(self):
        pass

    class Container():
        def __init__(self, self_range):
            self.scaffold = self_range.scaffold
            self.start = self_range.start
            self.end = self_range.end
            if self.end < self.start:
                self.start, self.end = self.end, self.start
            self.len = self.end - self.start + 1
            self.positions = []
            self.coverage = 0
            self.pccov = 0
            self.contained = None

        def add(self, hit):
            start, end = hit.start, hit.end
            if start > end:
                start, end = end, start
            if end < self.start or start > self.end:
                return
            if start < self.start:
                start = self.start
            if end > self.end:
                end = self.end
            self.positions.append((start, 'start'))
            self.positions.append((end, 'end'))

        def __repr__(self):
            return '{:7d} {:5.2f} %. Contained? {} {}'.format(self.coverage, self.pccov, self.contained, self.overlaps)

        @property
        def overlaps(self):
            overlaps = []
            within = 0
            start = 0
            for position in sorted(self.positions, key = lambda x:x[0]):
                if position[1] == 'start':
                    if not within:
                        start = position[0]
                    within += 1
                elif position[1] == 'end':
                    within -= 1
                    if not within:
                        overlaps.append((start, position[0]))
            return overlaps
            
        def get_coverage(self):
            begin = False
            end = False
            for overlap in self.overlaps:
                if overlap[0] <= self.start + self.len*0.15:
                    begin = True
                if overlap[1] >= self.end - self.len*0.15:
                    end = True
                self.coverage += overlap[1] - overlap[0] + 1

            self.pccov = self.coverage / self.len * 100
            self.contained = begin and end and self.pccov > 75
            return self.coverage, self.pccov, self.contained

    def pass_scaffold_pair(self, scaffold1, scaffold2):
        return scaffold1 == scaffold2 or scaffold1 not in self.hits or scaffold2 not in self.hits[scaffold1]

    
    def within(self, hit, raft):
        h_start, h_end = self.order(hit.start, hit.end)
        r_start, r_end = self.order(raft.start, raft.end)
        if h_end < r_start or r_end < h_start:
            return False
        else:
            return True
        

    def bridge(self, a, b):
        if a == b:
            return False

        a_hooks = a.get_hooks('last')
        b_hooks = b.get_hooks('first')
        print(a)
        print(b)
        for a_hook in a_hooks:
            for b_hook in b_hooks:
                if self.pass_scaffold_pair(a_hook.scaffold, b_hook.scaffold):
                    continue
                
                a_container = self.Container(a_hook)
                b_container = self.Container(b_hook)

                for hit in (self.Hit(h) for h in self.hits[a_hook.scaffold][b_hook.scaffold]):
                    if self.within(hit.r, a_hook) and self.within(hit.q, b_hook):
                        a_container.add(hit.r)
                        b_container.add(hit.q)
                if not a_container.positions and not b_container.positions:
                    continue

                a_container.get_coverage()
                b_container.get_coverage()
                
                print(a_container)
                print(b_container)
                if a_container.contained:
                    a.discard(b_hook)
                    return False
                if b_container.contained:
                    b.discard(a_hook)
                    return False

        return False


class PacBioMerge(Merger):
    def __init__(self, this, other, tablename):

        self.this = this
        self.other = other
        self.genome = this.genome
        self.chromosome = this.chromosome

        self.tablename = tablename if tablename is not None else 'pacbio_overlaps'
        self.bridges = defaultdict(lambda: defaultdict(list))

        self.genome_scaffolds = this.scaffolds + other.scaffolds

        self.hits, self.pb_lengths = self.get_hits()

        self.bridges = self.get_bridges()

    def get_hits(self):
        hits = defaultdict(lambda: defaultdict(list))
        pb_lengths = {}
        statement = 'select * from {} where qscaffold in ({})'.format(self.tablename, ','.join('?'*len(self.genome_scaffolds)))
        for args in self.chromosome.overlaps.execute(statement, self.genome_scaffolds):
            hit = self.Hit(args)
            pb_lengths[hit.pb.scaffold] = hit.pb.seqlen
            hits[hit.pb.scaffold][hit.g.scaffold].append(hit)
        
        to_remove = [pacbio for pacbio in hits if len(hits[pacbio]) <= 1]
        for pacbio in to_remove:
            del hits[pacbio]

        return hits, pb_lengths

    class Hit:
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
            return '''{0:16s} ({1:7d} bp) {2:7d} - {3:7d} ({4:7d} - {5:7d} ({6:7d} bp); {7:6.2f}%)'''.format(
                self.g.scaffold, self.g.seqlen, self.g.start, self.g.end, self.pb.start, self.pb.end, self.pb.hitlen, self.pcid
            )

    def get_bridges(self):
        bridges = defaultdict(lambda: defaultdict(list))
        for pacbio in self.hits:
            for gs1 in self.hits[pacbio]:
                for gs2 in self.hits[pacbio]:
                    for hit1 in self.hits[pacbio][gs1]:
                        for hit2 in self.hits[pacbio][gs2]:
                                bridges[gs1][gs2].append(self.Bridge(pacbio, self.pb_lengths[pacbio], hit1, hit2))
        return bridges


    class Bridge:
        
        def __init__(self, pacbio, pb_len, hit1, hit2):
            self.pacbio = pacbio
            self.pb_len = pb_len
            self.hit1 = hit1
            self.hit2 = hit2
            if hit1.pb.start > hit2.pb.start:
                self.hit1, self.hit2 = self.hit2, self.hit1
            self.gap = self.hit2.pb.start - self.hit1.pb.end - 1
            self._score = None
    
        @property
        def score(self):
            if self._score is None:
                self._score = sqrt((self.hit1.pb.hitlen * self.hit1.pcid/100) ** 2 + (self.hit2.pb.hitlen * self.hit2.pcid/100)  ** 2)
            
            return self._score

        def __repr__(self):
            out = ("{0} ({1:6d} bp) {2} - {3}\n\t{4}\n\t{5}".format(
                       self.pacbio, self.pb_len, self.hit1.pb.start, self.hit2.pb.end,
                       self.hit1, self.hit2
                   ))

            if self.hit1.pb.end < self.hit2.pb.start:
                out += "\tFill gap with {} {}-{}\t{} bp".format(self.pacbio, self.hit1.pb.end+1, self.hit2.pb.start-1, self.gap)
            elif self.hit1.pb.start <= self.hit2.pb.start and self.hit1.pb.end >= self.hit2.pb.end:
                out += "\tSecond hit contained in first hit"
            elif self.hit2.pb.start <= self.hit1.pb.start and self.hit2.pb.end >= self.hit1.pb.end:
                out += "\tFirst hit contained in second hit"
            else:
                out += "\tOverlap of {} bases".format(self.hit1.pb.end - self.hit2.pb.start + 1)
            return out
            
        def order(self, start, end):
            if start < end:
                return start, end
            else:
                return end, start
        
        def get_hit(self, end_range):

            r_start, r_end = self.order(end_range.start, end_range.end)
            hit = None
            
            if self.hit1.g.scaffold == end_range.scaffold:
                h_start, h_end = self.order(self.hit1.g.start, self.hit1.g.end)
                if r_start < h_end and h_start < r_end:
                    hit = self.hit1

            elif self.hit2.g.scaffold == end_range.scaffold:
                h_start, h_end = self.order(self.hit2.g.start, self.hit2.g.end)
                if r_start < h_end and h_start < r_end:
                    hit = self.hit2
            
            return hit
        
        def get_hit_overhangs(self, range1, range2):
            hit1_overhang = hit2_overhang = None
            if (range1.start < range1.end and range1.end <= self.hit1.g.end or
                range1.start > range1.end and range1.end >= self.hit1.g.end):
                    hit1_overhang = abs(self.hit1.g.end - range1.end)
            else:
                    hit1_overhang = abs(range1.end - self.hit1.g.end) * -1

            if (range2.start < range2.end and range2.end <= self.hit2.g.start or
                range2.start > range2.end and range2.end >= self.hit2.g.start):
                    hit2_overhang = abs(self.hit2.g.start - range2.end)
            else:
                    hit2_overhang = abs(range2.end - self.hit2.g.start) * -1
            
            if self.hit1.pb.end < self.hit2.pb.start:
                pb_range = r.Range(self.hit1.pb.scaffold, self.hit1.pb.end+1, self.hit2.pb.start-1)
            else:
                pb_range = r.Range(self.hit1.pb.scaffold, self.hit1.pb.end, self.hit2.pb.start)
            pb_overhang = self.hit2.pb.start - self.hit1.pb.end - 1
    
            return hit1_overhang, hit2_overhang, pb_overhang, pb_range

    def skip_range_pair(self, range1, range2):
        r1start, r1end = self.order(range1.start, range1.end)
        r2start, r2end = self.order(range2.start, range2.end)
        if range1.scaffold == range2.scaffold and r1start == r1end and r2start == r2end:
            return True
        if range1.scaffold not in self.bridges or range2.scaffold not in self.bridges[range1.scaffold]:
            return True
        return False

    def get_range_dict(self, ranges):
        range_dict = defaultdict(list)
        for r in ranges:
            new = deepcopy(r)
            if new.end < new.start:
                new.reverse()
            range_dict[new].append(r)
        return range_dict
        
    def add_pacbio_block(self, pb_range, raft):
        new_pb = gd.Block(pb_range.scaffold, pb_range.start, pb_range.end)
        self.genome.blocks[pb_range.scaffold][pb_range.start] = new_pb

        new_pb_direction = 1 if pb_range.start < pb_range.end else -1
        raft.append(new_pb.scaffold, new_pb.start, new_pb_direction)
        return

    def bridge_rafts(self, a, b, rb):
        a.trim(rb.self_range, rb.self_overhang)
        b.trim(rb.other_range, rb.other_overhang)

        self.add_pacbio_block(rb.pb_range, a)

        a.merge(b)
        b.empty()
        return

    
    def bridge(self, a, b):
        if a == b:
            return None

        a_ranges = self.get_range_dict(a.get_ranges('last'))
        b_ranges = self.get_range_dict(b.get_ranges('first'))
        
        for a_range_class in a_ranges:
            for a_range in a_ranges[a_range_class]:
                for b_range_class in b_ranges:
                    for b_range in b_ranges[b_range_class]:
                        if self.skip_range_pair(a_range, b_range):
                            continue
                
                        for bridge in self.bridges[a_range.scaffold][b_range.scaffold]:
                            a_hit = bridge.get_hit(a_range)
                            b_hit = bridge.get_hit(b_range)
                            if not a_hit or not b_hit:
                                continue
                            a.add_bridge(b, a_hit, a_range, b_range, bridge)
        return None

    def merge(self):
        for a in self.this.rafts:
            for a_range in a.bridges:
                for b_range in a.bridges[a_range]:
                    for pb_range in a.bridges[a_range][b_range]:
                        print(a.bridges[a_range][b_range][pb_range])
                

class OKMerge(Merger):
    def __init__(self, this, other, options=None):
        self.this = this
        self.other = other
        self.genome = this.genome
        
    def bridge(self, a, b):
        if self.this.pooltype != 'ok' or self.other.pooltype != 'ok':
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