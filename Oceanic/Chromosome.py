#!/usr/bin/env python3

from . import Pool as pl
from . import Raft as r
from . import Stats as s
from . import Mergers as merge

class Chromosome:
    def __init__(self, name, genome):
        self.name = str(name)
        self.markers = {}
        self.set_markers(genome.db)
        self.pools = []
        self.genome = genome
        self.mapped_blocks, self.mapped_blocks_length, self.placed_blocks, self.placed_blocks_length = self.set_blocks()
        
    def __repr__(self):
        return('\n'.join([self.name, repr(self.pools), repr(self.markers)]))

    def __iter__(self):
        return iter(self.pools)

    @property
    def stats(self):
        stats = s.Stats(self.name)
        stats.pool_num = len(self.pools)
        for pool in self.pools:
            gt = pool.pooltype
            for raft in pool:
                stats.log_num += len(raft.logs)
                stats.pool_stats[gt].length += raft.length
                stats.log_count[gt] += len(raft.logs)
                stats.raft_count[gt] += 1
            stats.pool_stats[gt].count += 1
        
        return stats

    def set_markers(self, db):
        prev_cm = -1
        next_cm = -1
        for cm, in db.execute("select distinct cm from scaffold_map where chromosome={} order by cm".format(self.name)):

            if cm == -1:
                self.add_marker(-1, -1, -1)
                continue

            self.add_marker(cm, prev_cm=prev_cm)
            if prev_cm != -1:
                self.update_marker(prev_cm, next_cm = cm)
        
            prev_cm = cm

    def set_blocks(self):
        mapped_blocks = mapped_blocks_length = placed_blocks = placed_blocks_length = 0

        for cm in sorted(self.markers) + [-1]:
            cm_blocks = pl.Pool(self)
            cm_block_id = 1
            statement = "select scaffold, start, end, length from scaffold_map where chromosome={} and cm={} order by scaffold, start, end".format(self.name, cm)
            for scaffold, start, end, length in self.genome.db.execute(statement):
                if scaffold not in self.genome.sequences:
                    continue
                if (scaffold, start) in self.genome.errors:
                    continue
                mapped_blocks += 1
                mapped_blocks_length += length
                self.genome.blocks[scaffold][start].add_marker(self.name,cm)

                if cm != -1:
                    placed_blocks += 1
                    placed_blocks_length += length
                cm_blocks.add(r.Raft(cm_block_id, scaffold, start, self))
                cm_block_id += 1
            if cm != -1:
                self.pools.append(cm_blocks)

        return mapped_blocks, mapped_blocks_length, placed_blocks, placed_blocks_length

    @property
    def marker_chains(self):
        
        marker_chains = {}
        
        for pool in chromosome:
            if pool.marker_chain not in marker_chains:
                marker_chains[pool.marker_chain] = {}
                marker_chains[pool.marker_chain]['count'] = 0
                marker_chains[pool.marker_chain]['type'] = ''

            marker_chains[pool.marker_chain]['type'] = ''
            marker_chains[pool.marker_chain]['count'] += 1
        
        return marker_chains

    def add_marker(self, cm, prev_cm=-1, next_cm=-1):
        self.markers[cm] = Marker(cm, prev_cm, next_cm)
        
    def update_marker(self, cm, prev_cm = -1, next_cm = -1):
        if cm not in self.markers:
            self.add_marker(cm, prev_cm, next_cm)

        self.markers[cm].update_previous(prev_cm)
        self.markers[cm].update_next(next_cm)

    def get_scaffolds(self):
        scaffolds = []
        for pool in self:
            for raft in pool:
                scaffoldlist = []
                
                for scaffold, start, direction in raft.logs:
                    scaffoldlist.append(self.genome.blocks[scaffold][start])
                scaffolds.append(scaffoldlist)
        
        return scaffolds
    
    def run_merger(self, merger):
        for pool in self.pools:
            pool.assemble(pool, merger)
        self.connect(merger)
    
    def assemble(self, genome, threads):
        self.run_merger(merge.cm)

        for pool in self.pools:
            pool.extend()
            
        self.run_merger(merge.overlap)

        self.run_merger(merge.pacbio)

        print(self.stats)

#        extend_ok_pools(self, overlaps)


    def connect(self, merger):
        p = 0
        while p < len(self.pools)-1:
            self.pools[p].assemble(self.pools[p+1], merger)
            p = self.split(p)
            if self.pools[p+1]:
                p += 1
            else:
                del self.pools[p+1]


    def split(self, p):
        ordered_rafts = [raft for raft in self.pools[p] if raft.ordered]
        if not ordered_rafts or len(ordered_rafts) == len(self.pools[p]):
            return p

        self.pools.insert(p+1, pl.Pool(self))
        for raft in ordered_rafts:
            self.pools[p+1].add(raft)
            self.pools[p].remove(raft)
        return p+1


class Marker:
    def __init__(self, cm, prev_cm=-1, next_cm=-1):
        self.cm = cm
        self.prev_cm = prev_cm
        self.next_cm = next_cm

    def __repr__(self):
        return('{}-({},{})'.format(self.cm, self.prev_cm, self.next_cm))

    def update_previous(self, prev_cm):
        if prev_cm != -1:
            self.prev_cm = prev_cm
    
    def update_next(self, next_cm):
        if next_cm != -1:
            self.next_cm = next_cm


if __name__ == '__main__':
    
    pass