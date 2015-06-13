#!/usr/bin/env python3

from itertools import chain

class Pool:
    def __init__(self, chromosome):
        self.chromosome = chromosome
        self.genome = chromosome.genome
        self.rafts = set()

    def __repr__(self):
        output = 'Chromosome {}\nType: {}\n'.format(self.chromosome.name, self.pooltype)
        for raft in self:
            output += repr(raft)
            output += '-----\n'
        output += '====='
        
        return output
    
    def __iter__(self):
        return iter(self.rafts)
    
    def __len__(self):
        return len(self.rafts)

    def add(self, raft):
        self.rafts.add(raft)

    def remove(self, raft):
        self.rafts.remove(raft)
    
    @property
    def scaffolds(self):
        return [scaffold for r in self.rafts for scaffold in r.scaffolds]
    
    @property
    def marker_chain(self):
        for raft in sorted(self.rafts, key=lambda x: x.length, reverse=True):
            return raft.marker_chain
        return ''
    
    @property
    def pooltype(self):
        if len(self.marker_chain) == 1:
            if len(self.rafts) > 1:
                return 'order'
            else:
                return 'orient'
        else:
            if len(self.rafts) > 1:
                return 'overlap'
            else:
                return 'ok'

    def cleanup(self):
        for raft in [raft for raft in self.rafts if not raft.logs]:
            self.rafts.remove(raft)

    def assemble(self, other, mergeclass, options=None):
        merger = mergeclass(self, other, options)
        for a in self.rafts:
            if not a:
                continue

            repeat = True
            seen = {}
            while repeat:
                repeat = False
                for b in other.rafts:
                    if not b or repr(a) == repr(b) or (a,b) in seen:
                        continue

                    seen[(a,b)] = 1
                    if merger.bridge(a, b):
                        repeat = True
                        break
        merger.merge()
        other.cleanup()


    def extend(self):
        for raft in self.rafts:
            raft.extend()
    
    def tie(self, hit, rope):
        for raft in self.rafts:
            knot = raft.tie(hit, rope)
            if knot:
                return knot
        else:
            return None
    
    def turn_hooks(self, end):
        for raft in self:
            raft.turn_hooks(end)

    def write(self):
        scaffolds = []
        for raft in self:
            scaffolds.append(raft.write())
        return scaffolds

if __name__ == '__main__':
    pass
