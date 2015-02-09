#!/usr/bin/env python3

from copy import deepcopy
from collections import defaultdict
from math import sqrt
from operator import itemgetter

from . import Raft as r

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

    def merge(self, a):
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
            
            if start == target:
                bridge = self.merge_logs(bridge, a.logs)
                bridge = self.merge_logs(bridge, b.logs)
                self.do_join(a, b, bridge)
                return True
                
            if start == 0:
                break
            
            block_i_chr = scaffold[start].chromosome
            block_i_cm = scaffold[start].cm
            if block_i_chr == '0' or block_i_cm == -1 and block_i_chr == a.chromosome.name:
                bridge.append((a.scaffold, start, direction))
                continue
        
            if block_i_chr != a.chromosome.name: # This check comes after the cM check because chr==0 is OK
                break
        
            if block_i_cm == this_cm:
                bridge.append((a.scaffold, start, direction))
        
            if block_i_cm == a.chromosome.markers[this_cm].next_cm:
                this_cm = a.chromosome.markers[this_cm].next_cm
                bridge.append((a.scaffold, start, direction))
        
            
            if block_i_cm == a.chromosome.markers[this_cm].prev_cm:
                this_cm = a.chromosome.markers[this_cm].prev_cm
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
        return a_near, b_near, a_forward, b_forward
    
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

    def merge(self, a):
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

        a_ranges = a.get_ranges('last')
        b_ranges = b.get_ranges('first')

        for a_range in a_ranges:
            for b_range in b_ranges:
                if self.pass_scaffold_pair(a_range.scaffold, b_range.scaffold):
                    continue
                
                a_container = self.Container(a_range)
                b_container = self.Container(b_range)

                for hit in (self.Hit(h) for h in self.hits[a_range.scaffold][b_range.scaffold]):
                    if self.within(hit.r, a_range) and self.within(hit.q, b_range):
                        a_container.add(hit.r)
                        b_container.add(hit.q)
                if not a_container.positions and not b_container.positions:
                    continue

                a_container.get_coverage()
                b_container.get_coverage()
                
                if a_container.contained:
                    a.discard(b_range)
                    return False
                if b_container.contained:
                    b.discard(a_range)
                    return False

        return False


class PacBioMerge(Merger):
    def __init__(self, this, other, tablename):

        self.this = this
        self.other = other
        
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
            range_dict[repr(new)].append(r)
        return range_dict
        
        
    def bridge(self, a, b):
        if a == b:
            return None

        a_ranges = self.get_range_dict(a.get_ranges('last'))
        b_ranges = self.get_range_dict(b.get_ranges('first'))
        
        for a_range_i in a_ranges:
            for a_range in a_ranges[a_range_i]:
                for b_range_i in b_ranges:
                    for b_range in b_ranges[b_range_i]:
                        if self.skip_range_pair(a_range, b_range):
                            continue
                
                        for bridge in self.bridges[a_range.scaffold][b_range.scaffold]:
                            a_hit = bridge.get_hit(a_range)
                            b_hit = bridge.get_hit(b_range)
                            if not a_hit or not b_hit:
                                continue
                            a.add_bridge(a_hit, a_range, a_range_i, b_range, b_range_i, bridge)

        for a_range in sorted(a.bridges):
            merged = False
            for b_range in sorted(a.bridges[a_range]):
                ab = a.bridges[a_range][b_range]
                best_hit = ab[sorted(ab, key = lambda x:(-ab[x].coverage, ab[x].overhang))[0]]
                if best_hit.bridge.hit1.g.hitlen < 6000 or best_hit.bridge.hit2.g.hitlen < 6000 and best_hit.overhang > -15000:
                    continue

                if not a.ordered and best_hit.self_range.end == a.start:
                    a.reverse()
                if not b.ordered and best_hit.other_range.end == b.end:
                    b.reverse()

                if best_hit.self_range.end == a.end and best_hit.other_range.end == b.start:
                    new_chain = a.marker_chain + b.marker_chain
                    if new_chain:
                        a.merge(b)
                        b.empty()
                        merged=True
                        break
                        
            if merged:
                break

        return None
    
    def merge(self, a):
        pass

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
        
    def merge(self, a):
        pass


if __name__ == '__main__':
    pass