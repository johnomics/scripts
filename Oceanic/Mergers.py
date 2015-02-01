#!/usr/bin/env python3

from copy import deepcopy
from collections import defaultdict

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

        self.hits = defaultdict(lambda:defaultdict(list))

        statement = '''select * from genome_overlaps
                       where rscaffold in ({})
                       and   qscaffold in ({})
                    '''.format(','.join('?'*len(this.scaffolds)), ','.join('?'*len(other.scaffolds)))

        for *args, rscaffold, qscaffold, hittype in this.genome.overlaps.execute(statement, this.scaffolds + other.scaffolds):
            self.hits[rscaffold][qscaffold].append(args)

    class Hit:
        def __init__(self, hit):
            self.rstart, self.rend, self.qstart, self.qend, \
            self.rhitlen, self.qhitlen, self.pcid, self.rseqlen, self.qseqlen, \
            self.rpccov, self.qpccov = hit

    def merge(self, a):
        pass

    def bridge(self, a, b):
        if (a == b or
           a.scaffold not in self.hits or
           b.scaffold not in self.hits[a.scaffold]):
            return None
        
        for hit in (self.Hit(h) for h in self.hits[a.scaffold][b.scaffold]):
            pass
#            print(a.scaffold, b.scaffold, hit.rstart, hit.rend, hit.qstart, hit.qend, hit.rhitlen, hit.qhitlen)

        return None


class PacBioMerge():
    def __init__(self, this, other):
        self.genome = this.genome
        self.hits = defaultdict(lambda: defaultdict(list))
        self.pacbio_lengths = {}
        self.genome_scaffolds = this.scaffolds + other.scaffolds
        statement = 'select * from pacbio_overlaps where qscaffold in ({})'.format(','.join('?'*len(self.genome_scaffolds)))
        for args in self.genome.overlaps.execute(statement, self.genome_scaffolds):
            hit = self.Hit(args)
            self.pacbio_lengths[hit.rscaffold] = hit.rseqlen
            self.hits[hit.rscaffold][hit.qscaffold].append(hit)

        self.bridges = self.get_bridges()

    class Hit:
        def __init__(self, hit):
            self.rstart, self.rend, self.qstart, self.qend, \
            self.rhitlen, self.qhitlen, self.pcid, self.rseqlen, self.qseqlen, \
            self.rpccov, self.qpccov, self.rscaffold, self.qscaffold, self.hittype = hit
    
    def get_bridges(self):
        bridges = {}
        for pbs in self.hits:
            if len(self.hits[pbs]) <= 1:
                continue
            pacbio_order = self.order_hits(self.hits[pbs])
    
            for s1, hit1 in sorted(pacbio_order.items()):
                for s2, hit2 in sorted(pacbio_order.items()):
                    if s2 <= s1:
                        continue
                    self.add_bridge(bridges, pbs, self.pacbio_lengths[pbs], hit1, hit2)
    
        return bridges

    def order_hits(self, hits):
        pacbio_order = {}
        for gs in hits:
            for hit in sorted(hits[gs], key=lambda x:x.qhitlen, reverse=True):
                if hit.rstart > hit.rend:
                    hit.qstart, hit.qend = hit.qend, hit.qstart
                    hit.rstart, hit.rend = hit.rend, hit.rstart
                ignore = False
                for start, pbhit in pacbio_order.items():
                    if gs == pbhit['scaffold'] and hit.qstart >= pbhit['gstart'] and hit.qend <= pbhit['gend']:
                        ignore = True
                if not ignore:
                    pacbio_order[hit.rstart] = {'scaffold': gs, 'pbstart': hit.rstart, 'pbend': hit.rend,
                                                'gstart': hit.qstart, 'gend': hit.qend, 'pcid':hit.pcid}
        return pacbio_order

    
    def bridge(self, a, b):
#        if a.scaffold in self.bridges and b.scaffold in self.bridges[a.scaffold]:
#            print(self.bridges[a.scaffold][b.scaffold])

        return None
    
    def merge(self, a):
        pass

    def add_bridge(self, bridges, pbs, pblen, hit1, hit2):
        for hit in hit1, hit2:
            if hit['scaffold'] not in bridges:
                bridges[hit['scaffold']] = defaultdict(list)
        hit1seqlen = len(self.genome.sequences[hit1['scaffold']])
        hit2seqlen = len(self.genome.sequences[hit2['scaffold']])
        bridges[hit1['scaffold']][hit2['scaffold']].append(self.Bridge(pbs, pblen, hit1, hit2, hit1seqlen, hit2seqlen))
        bridges[hit2['scaffold']][hit1['scaffold']].append(self.Bridge(pbs, pblen, hit2, hit1, hit2seqlen, hit1seqlen))

    class Bridge:
        
        def __init__(self, pbs, pblen, hit1, hit2, hit1seqlen, hit2seqlen):
            self.pbs = pbs
            self.pblen = pblen
            self.hit1 = hit1
            self.hit2 = hit2
            self.hit1seqlen = hit1seqlen
            self.hit2seqlen = hit2seqlen
            self.gap = self.hit2['pbstart']-1 - self.hit1['pbend']+1
    
        def __repr__(self):
            out = ('''{0} ({1:6d} bp) {2} - {3} joins
           {4:16s} ({5:7d} bp) {6:7d} - {7:7d} ({2:7d} - {8:7d}; {9:5.2f}%) and
           {10:16s} ({11:7d} bp) {12:7d} - {13:7d} ({14:7d} - {3:7d}; {15:5.2f}%)'''.format(
                       self.pbs, self.pblen, self.hit1['pbstart'], self.hit2['pbend'],
                       self.hit1['scaffold'], self.hit1seqlen,
                       self.hit1['gstart'], self.hit1['gend'], self.hit1['pbend'], self.hit1['pcid'],
                       self.hit2['scaffold'], self.hit2seqlen,
                       self.hit2['gstart'], self.hit2['gend'], self.hit2['pbstart'], self.hit2['pcid']
                   ))
                
            if self.hit1['pbend'] < self.hit2['pbstart']:
                out += "\tFill gap with {} {}-{}".format(self.pbs, self.hit1['pbend']+1, self.hit2['pbstart']-1)
            elif self.hit1['pbstart'] < self.hit2['pbstart'] and self.hit1['pbend'] > self.hit2['pbend']:
                out += "\tSecond hit contained in first hit"
            elif self.hit2['pbstart'] < self.hit1['pbstart'] and self.hit2['pbend'] > self.hit1['pbend']:
                out += "\tFirst hit contained in second hit"
            else:
                out += "\tOverlap of {} bases".format(self.hit1['pbend'] - self.hit2['pbstart'] + 1)
            return out


def get_neighbour(group, bridges, direction):
    neighbour = None
    (thisgroup, ) = group
    if direction == 1:
        this = thisgroup.manifest[-1]
    else:
        this = thisgroup.manifest[0]
    if this.scaffold in bridges:
        for a in bridges[this.scaffold]:
            if this.scaffold != a:
                print(this.scaffold, a)
                for bridge in sorted(bridges[this.scaffold][a],key=lambda x:x.gap):
                    if not neighbour:
                        neighbour = bridge
                    print(bridge)
    return neighbour

def extend_ok(i, direction, pools):
    pacbio_bridges = get_pacbio_overlaps(pool_scaffolds, overlaps)
    
    for a in pacbio_bridges:
        for b in pacbio_bridges[a]:
            print(a, b)
            for bridge in pacbio_bridges[a][b]:
                print(bridge)
    return i
    
    neighbour = get_neighbour(pools[i].rafts, pacbio_bridges, direction)
    if not neighbour:
        return i

    print("Neighbour:")
    print(neighbour)
    (this,) = pools[i].rafts
    to_remove = set()
    rafts = [raft for raft in pools[i+direction].rafts]
    for raft in rafts:
        if raft.scaffold == neighbour.hit2['scaffold']:
            print('Neighbour:', raft)
            for (scaffold, start, direction) in raft.logs:
                if direction == 1:
                    this.append(scaffold, start, direction)
                else:
                    this.prepend(scaffold, start, direction)
            print('i:{} dir:{}'.format(i, direction))
            print("Pools are currently:")
            print(pools[i+direction].rafts)
            print("And raft is:")
            print(raft)
            pools[i+direction].rafts.remove(bl)

#    for raft in to_remove:
#        pools[i+direction].rafts.remove(bl)

    print("Pools now:")
    print(pools[i])
    print(pools[i+direction])
    print("=======")

    return i


if __name__ == '__main__':
    pass