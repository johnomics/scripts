#!/usr/bin/env python3

class Pool:
    def __init__(self, chromosome):
        self.chromosome = chromosome
        self.genome = chromosome.genome
        self.rafts = set()

    def __repr__(self):
        output = 'Type: {}\n'.format(self.pooltype)
        for raft in self:
            output += repr(raft) + '\n'
            output += 'Length: {}\n'.format(raft.length)
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
    def marker_chain(self):
        for raft in self.rafts:
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

    def assemble(self, other, merger):
        for a in self.rafts:
            if not a:
                continue

            repeat = True
            while repeat:
                repeat = False
                for b in other.rafts:
                    if not b or repr(a) == repr(b):
                        continue

                    merge = merger(a, b)
                    if merge:
                        repeat = True
                        a.replace(merge)
                        b.empty()
                        break
        other.cleanup()

    def extend(self):
        for raft in self.rafts:
            raft.extend()


if __name__ == '__main__':
    pass
