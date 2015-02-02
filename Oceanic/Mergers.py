#!/usr/bin/env python3

from copy import deepcopy
from collections import defaultdict
from math import sqrt

class MarkerMerge():
    
    def __init__(self, this, other):
        self.joins = []

    def merge(self, a):
        pass

    def do_join(self, a, b, bridge):
        a.replace(bridge)
        b.empty()

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
    


class OverlapMerge():
    
    def __init__(self, this, other):
        
        print(this)
        if this is not other:
            print(other)
        self.hits = defaultdict(lambda:defaultdict(list))

        statement = '''select * from genome_overlaps
                       where rscaffold in ({})
                       and   qscaffold in ({})
                    '''.format(','.join('?'*len(this.scaffolds)), ','.join('?'*len(other.scaffolds)))

        for *args, rscaffold, qscaffold, hittype in this.genome.overlaps.execute(statement, this.scaffolds + other.scaffolds):
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
        def __init__(self):
            self.positions = []
            self.coverage = 0
            self.pccov = 0

        def add(self, hit):
            self.seqlen = hit.seqlen
            start, end = hit.start, hit.end
            if start > end:
                start, end = end, start
            self.positions.append((start, 'start'))
            self.positions.append((end, 'end'))

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
                if overlap[0] == 1:
                    begin = True
                if overlap[1] == self.seqlen:
                    end = True
                self.coverage += overlap[1] - overlap[0] + 1

            self.pccov = self.coverage / self.seqlen * 100
            self.contained = begin and end
            return self.coverage, self.pccov, self.contained

    def contained(self, scaffold, scftype, hits):
        container = self.Container()
        for hit in (self.Hit(h) for h in hits):
#            print("{} ({} bp) {}-{} ({} bp) hits {} ({} bp) {}-{} ({} bp)".format(
#                hit.r.scaffold, hit.r.seqlen, hit.r.start, hit.r.end, hit.r.hitlen, 
#                hit.q.scaffold, hit.q.seqlen, hit.q.start, hit.q.end, hit.q.hitlen))
            if scftype == 'r':
                container.add(hit.r)
            else:
                container.add(hit.q)
        coverage, pccov, contained = container.get_coverage()
#        print("{} ({} bp) covered {} bp, {} %, {}".format(scaffold, container.seqlen, coverage, pccov, contained))
        
        if contained and pccov > 80:
            return True

        return False


    def bridge(self, a, b):
        if (a == b or
           a.scaffold not in self.hits or
           b.scaffold not in self.hits[a.scaffold]):
            return None
        
        if self.contained(a.scaffold, 'r', self.hits[a.scaffold][b.scaffold]):
            print("{} contained in {}!".format(a.scaffold, b.scaffold))
#            a.empty()
            return True
        if self.contained(b.scaffold, 'q', self.hits[a.scaffold][b.scaffold]):
            print("{} contained in {}!".format(b.scaffold, a.scaffold))
 #           b.empty()
            return True

        return False


class PacBioMerge():
    def __init__(self, this, other):
        self.this = this
        self.other = other
        self.genome = this.genome
        self.bridges = defaultdict(lambda: defaultdict(list))
        
        self.genome_scaffolds = this.scaffolds + other.scaffolds
        
        self.hits, self.pb_lengths = self.get_hits()

        self.bridges = self.get_bridges()

    def get_hits(self):
        hits = defaultdict(lambda: defaultdict(list))
        pb_lengths = {}
        statement = 'select * from pacbio_overlaps where qscaffold in ({})'.format(','.join('?'*len(self.genome_scaffolds)))
        for args in self.genome.overlaps.execute(statement, self.genome_scaffolds):
            hit = self.Hit(args)
            pb_lengths[hit.pb.scaffold] = hit.pb.seqlen
            hits[hit.pb.scaffold][hit.g.scaffold].append(hit)
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

    def match_end(self, hit):
        return hit.g.start == 1 or hit.g.end == 1 or hit.g.start == hit.g.seqlen or hit.g.end == hit.g.seqlen

    def get_bridges(self):
        bridges = defaultdict(lambda: defaultdict(list))
        for pacbio in self.hits:
            if len(self.hits[pacbio]) <= 1: # If PacBio scaffold only hits one genome scaffold, it's useless for bridging
                continue
            for gs1 in self.hits[pacbio]:
                for gs2 in self.hits[pacbio]:
                    if gs1 == gs2:
                        continue
                    for hit1 in self.hits[pacbio][gs1]:
#                        if not self.match_end(hit1):
#                            continue
                        for hit2 in self.hits[pacbio][gs2]:
#                            if not self.match_end(hit2):
#                                continue
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
            self.gap = self.hit2.pb.start - self.hit1.pb.end + 1
            self._score = None
            self.assembly = self.get_assembly()
    
        @property
        def score(self):
            if self._score is None:
                self._score = sqrt((self.hit1.pb.hitlen * self.hit1.pcid/100) ** 2 + (self.hit2.pb.hitlen * self.hit2.pcid/100)  ** 2)
            
            return self._score

        def get_parts(self, hit):
            parts = [None, (hit.g.start, hit.g.end), None]
            if hit.g.start < hit.g.end:
                if hit.g.start != 1:
                    parts[0] = (1, hit.g.start-1)
                if hit.g.end != hit.g.seqlen:
                    parts[2] = (hit.g.end+1, hit.g.seqlen)
            elif hit.g.start > hit.g.end:
                if hit.g.end != 1:
                    parts[2] = (hit.g.end-1, 1)
                if hit.g.start != hit.g.seqlen:
                    parts[0] = (hit.g.seqlen, hit.g.start+1)
            return parts

        def get_assembly(self):
            self.g1_parts = self.get_parts(self.hit1)
            self.g2_parts = self.get_parts(self.hit2)
            self.hit1.g.lost = abs(self.g1_parts[2][0] - self.g1_parts[2][1]) + 1 if self.g1_parts[2] is not None else 0
            self.hit2.g.lost = abs(self.g2_parts[0][0] - self.g2_parts[0][1]) + 1 if self.g2_parts[0] is not None else 0
            self.total_lost = self.hit1.g.lost + self.hit2.g.lost
            return self.g1_parts, self.g2_parts

        def __repr__(self):
            out = ("{0} ({1:6d} bp) {2} - {3}\n\t{4}\n\t{5}".format(
                       self.pacbio, self.pb_len, self.hit1.pb.start, self.hit2.pb.end,
                       self.hit1, self.hit2
                   ))

            if self.hit1.pb.end < self.hit2.pb.start:
                out += "\tFill gap with {} {}-{}\t{} bp".format(self.pacbio, self.hit1.pb.end+1, self.hit2.pb.start-1, self.gap)
            elif self.hit1.pb.start < self.hit2.pb.start and self.hit1.pb.end > self.hit2.pb.end:
                out += "\tSecond hit contained in first hit"
            elif self.hit2.pb.start < self.hit1.pb.start and self.hit2.pb.end > self.hit1.pb.end:
                out += "\tFirst hit contained in second hit"
            else:
                out += "\tOverlap of {} bases".format(self.hit1.pb.end - self.hit2.pb.start + 1)
            
            out += "\n{}\n".format(self.assembly)
            

            out += "Remove {} bp from {} and {} bp from {} = {} bp".format(
                self.hit1.g.lost, self.hit1.g.scaffold, self.hit2.g.lost, self.hit2.g.scaffold, self.total_lost
            )
            return out
    
    def bridge(self, a, b):
        if a.scaffold != b.scaffold and a.scaffold in self.bridges and b.scaffold in self.bridges[a.scaffold]:
            print(a.scaffold, b.scaffold)
            for bridge in sorted(self.bridges[a.scaffold][b.scaffold], key=lambda x:x.total_lost)[:3]:
                print(bridge)

        return None
    
    def merge(self, a):
        pass

if __name__ == '__main__':
    pass