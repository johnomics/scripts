#!/usr/bin/env python3

from copy import deepcopy

def orient(a, b):
    edge_blocks = [a.start, a.last, b.start, b.last]
    if a.start < b.start:
        a_far, a_near, b_near, b_far = sorted(edge_blocks)
    else:
        b_far, b_near, a_near, a_far = sorted(edge_blocks)
    
    a_forward = a.start <= a.last
    b_forward = b.start <= b.last
    return a_near, b_near, a_forward, b_forward

def merge_logs(this, other):
    if not this:
        return deepcopy(other)

    new = deepcopy(this) + deepcopy(other)
    new.sort(key=lambda x: x[1]) # Sort by start position

    this_dir = this[0][2]
    new = [(scaffold, start, this_dir) for scaffold, start, direction in new]
    if this_dir == -1:
        new.reverse()

    return new

def cm(a, b):
    if a.scaffold != b.scaffold:
        return None
    bridge = []
    scaffold = a.genome.blocks[a.scaffold]

    start, target, a_forward, b_forward = orient(a,b)
    this_cm = scaffold[start].cm
    direction = 1 if start < target else -1

    while True:
        if direction == 1:
            start = scaffold[start].next_block
        else:
            start = scaffold[start].prev_block
        
        if start == target:
            bridge = merge_logs(bridge, a.logs)
            bridge = merge_logs(bridge, b.logs)
            return bridge
            
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

def overlap(a, b):
    return None

def pacbio(a, b):
    return None

def get_genome_overlaps(scaffolds, overlaps):
    found = False
    statement = '''select * from genome_overlaps
                   where rscaffold in ({0}) and qscaffold in ({0})
                     and rscaffold != qscaffold'''.format(
                ','.join('?'*len(scaffolds)))
    for rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rseqlen, qseqlen, rpccov, qpccov, \
        rscaffold, qscaffold, hittype in overlaps.execute(statement, list(scaffolds)*2):
        found = True
        print(rscaffold, len(genome.sequences[rscaffold]), rstart, rend, rhitlen, rpccov,
              qscaffold, len(genome.sequences[qscaffold]), qstart, qend, qhitlen, qpccov,
              pcid, hittype, sep='\t')
    
    if not found:
        print("No genome overlaps found")


def order_hit_ends(start, end):
    if start > end:
        start, end = end, start
    return start, end


def get_pacbio_hits(scaffolds, overlaps):
    pacbio_scaffolds = {}
    pacbio_lengths = {}
    statement = 'select * from pacbio_overlaps where qscaffold in ({})'.format(','.join('?'*len(scaffolds)))
    for rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rseqlen, qseqlen, rpccov, qpccov, rscaffold, qscaffold, hittype \
      in overlaps.execute(statement, list(scaffolds)):
        gstart, gend = order_hit_ends(qstart, qend)
        hits_to_keep = {}
        for start in scaffolds[qscaffold]:
            bstart, bend = order_hit_ends(start, scaffolds[qscaffold][start])
#            if (gstart > bstart and gend < bend
#               or gstart > bend or gend < bstart):
            if (gstart > bend or gend < bstart):
                continue
            hits_to_keep[gstart] = gend

        if hits_to_keep:
            if rscaffold not in pacbio_scaffolds:
                pacbio_scaffolds[rscaffold] = defaultdict(list)
                pacbio_lengths[rscaffold] = rseqlen
            pacbio_scaffolds[rscaffold][qscaffold].append((qstart, qend, qhitlen, rstart, rend, pcid))

    return pacbio_scaffolds, pacbio_lengths

def order_pacbio_hits(pbscaffold_hits):
    pacbio_order = {}
    for gs in pbscaffold_hits:
        for (qstart, qend, qhitlen, rstart, rend, pcid) in sorted(pbscaffold_hits[gs], key=lambda x:x[2], reverse=True):
            if rstart > rend:
                qstart, qend = qend, qstart
                rstart, rend = rend, rstart
            ignore = False
            for start, hit in pacbio_order.items():
                if gs == hit['scaffold'] and qstart >= hit['gstart'] and qend <= hit['gend']:
                    ignore = True
            if not ignore:
                pacbio_order[rstart] = {'scaffold': gs, 'pbstart': rstart, 'pbend': rend, 'gstart': qstart, 'gend': qend, 'pcid':pcid}
    return pacbio_order

class Bridge:
    
    def __init__(self, pbs, pblen, hit1, hit2):
        self.pbs = pbs
        self.pblen = pblen
        self.hit1 = hit1
        self.hit2 = hit2
        self.gap = self.hit2['pbstart']-1 - self.hit1['pbend']+1

    def __repr__(self):
        out = ('''{0} ({1:6d} bp) {2} - {3} joins
       {4:16s} ({5:7d} bp) {6:7d} - {7:7d} ({2:7d} - {8:7d}; {9:5.2f}%) and
       {10:16s} ({11:7d} bp) {12:7d} - {13:7d} ({14:7d} - {3:7d}; {15:5.2f}%)'''.format(
                   self.pbs, self.pblen, self.hit1['pbstart'], self.hit2['pbend'],
                   self.hit1['scaffold'], len(genome.sequences[self.hit1['scaffold']]),
                   self.hit1['gstart'], self.hit1['gend'], self.hit1['pbend'], self.hit1['pcid'],
                   self.hit2['scaffold'], len(genome.sequences[self.hit2['scaffold']]),
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

def add_bridge(bridges, pbs, pblen, hit1, hit2):
    for hit in hit1, hit2:
        if hit['scaffold'] not in bridges:
            bridges[hit['scaffold']] = defaultdict(list)
    bridges[hit1['scaffold']][hit2['scaffold']].append(Bridge(pbs, pblen, hit1, hit2))
    bridges[hit2['scaffold']][hit1['scaffold']].append(Bridge(pbs, pblen, hit2, hit1))
    
def get_pacbio_bridges(scaffolds, pacbio_scaffolds, pacbio_lengths):
    bridges = {}
    for pbs in pacbio_scaffolds:
        if len(pacbio_scaffolds[pbs]) <= 1:
            continue
        pacbio_order = order_pacbio_hits(pacbio_scaffolds[pbs])

        for s1, hit1 in sorted(pacbio_order.items()):
            for s2, hit2 in sorted(pacbio_order.items()):
                if s2 <= s1:
                    continue
                add_bridge(bridges, pbs, pacbio_lengths[pbs], hit1, hit2)

    return bridges


def get_pacbio_overlaps(scaffolds, overlaps):
    
    bridges = {}

    pacbio_scaffolds, pacbio_lengths = get_pacbio_hits(scaffolds, overlaps)
    
    if not pacbio_scaffolds:
        print("No PacBio overlaps found")
        return bridges

    bridges = get_pacbio_bridges(scaffolds, pacbio_scaffolds, pacbio_lengths)

    return bridges

def get_pool_scaffolds(scaffolds, pool, cms=None, fillcms=False):
    for raft in pool:
        for sb in raft.manifest:
            if fillcms and sb.cm != -1 and sb.cm not in cms:
                cms[sb.cm] = 0
            if sb.cm in cms:
                if sb.scaffold not in scaffolds:
                    scaffolds[sb.scaffold] = {}
                if raft.start not in scaffolds[sb.scaffold]:
                    scaffolds[sb.scaffold][raft.start] = raft.end

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
    pool_scaffolds = {}
    cms = {}
    if direction == 1:
        print('Testing pool against next')
        print(pools[i])
        print(pools[i+direction])
    else:
        print('Testing pool against previous')
        print(pools[i+direction])
        print(pools[i])
    get_pool_scaffolds(pool_scaffolds, pools[i], cms, fillcms=True)
    get_pool_scaffolds(pool_scaffolds, pools[i+direction], cms)
    
    print(pool_scaffolds)
    genome_overlaps = get_genome_overlaps(pool_scaffolds, overlaps)
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

def extend_ok_pools(chromosome, overlaps):
    i = 0
    while i < (len(chromosome.pools)):
        thispool = chromosome.pools[i]
        if thispool.pooltype != 'ok':
            i += 1
            continue
        
        if i > 0:
            i = extend_ok(i, -1, chromosome.pools)

        if i < len(chromosome.pools)-1:
            i = extend_ok(i, 1, chromosome.pools)
        
        i += 1


def join_ok_pools(chromosome):
    for i in range(0, len(chromosome.pools)-1):
        pool_i = chromosome.pools[i]
        j=i+1
        while j<len(chromosome.pools):
            pool_j = chromosome.pools[j]

            if pool_i.pooltype == 'ok' and pool_j.pooltype == 'ok':
                if len(pool_i.rafts) != 1 or len(pool_j.rafts) != 1:
                    print("OK pool has more than one raft!")
                    sys.exit()
                (pool_i_raft,)=pool_i.rafts
                (pool_j_raft,)=pool_j.rafts

                extend(pool_i_raft, -1)
                extend(pool_j_raft, 0)
                for (scaffold, start, direction) in pool_j_raft.logs:
                    pool_i_raft.append(scaffold, start, direction)
                pool_j.rafts.pop()
                pool_j.pooltype=''
            else:
                break

            j += 1

    newpools = []
    for pool in chromosome:
        if pool.rafts:
            newpools.append(pool)
    
    chromosome.pools = newpoolss


if __name__ == '__main__':
    pass