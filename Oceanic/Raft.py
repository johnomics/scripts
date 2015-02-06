#!/usr/bin/env python3

from collections import defaultdict
from copy import deepcopy

class Raft:
    def __init__(self, bl_id, scaffold, start, chromosome):
        self.chromosome = chromosome
        self.genome = chromosome.genome
        self.id = bl_id
        self.logs = []
        self.manifest = []
        self.append(scaffold, start, 1)
        self.bridges = defaultdict(lambda:defaultdict(lambda:defaultdict(list)))

    def __repr__(self):
        return '\n'.join([repr(m) for m in self.manifest])

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __iter__(self):
        return iter(self.logs)
    
    def __len__(self):
        return len(self.logs)
    
    def __contains__(self, start):
        for log in self.logs:
            if log[1] == start:
                return True
        return False

    @property
    def scaffolds(self):
        return list({m.scaffold:0 for m in self.manifest}.keys())

    @property
    def scaffold(self):
        return '_'.join(sorted(self.scaffolds))
        
    @property
    def name(self):
        summarydict = defaultdict(list)
        scaffolds = []
        for m in self.manifest:
            if m.scaffold not in summarydict:
                scaffolds.append(m.scaffold)
            summarydict[m.scaffold].append((m.start,m.end))
        names = []
        for scaffold in scaffolds:
            names.append('{}_{}_{}'.format(scaffold, summarydict[scaffold][0][0], summarydict[scaffold][-1][-1]))
        name = '-'.join(names)
        return name

    @property
    def start(self):
        return self.manifest[0].start
    
    @property
    def last(self):
        return self.logs[-1][1]

    @property
    def end(self):
        return self.manifest[-1].end

    @property
    def length(self):
        length = 0
        for m in self.manifest:
            length += m.length
        return length
    
    @property
    def sequence(self):
        return sum([m.sequence for m in self.manifest], Seq("", generic_dna))

    @property
    def marker_chain(self):
        marker_chain = []
        for m in self.manifest:
            if m.cm != -1:
                marker_chain.append(m.cm)
                if len(marker_chain)>1 and marker_chain[-2] == marker_chain[-1]:
                    del marker_chain[-1]

        if not self.check_chain(marker_chain):
            return ()
#            print("To fix:", marker_chain, "\n", self, "\n")

        return tuple(marker_chain)

    @property
    def ordered(self):
        return len(self.marker_chain) > 1

    @property
    def ranges(self):
        if not self._ranges:
            self._ranges = self.set_ranges()
        return self._ranges

    def set_ranges(self):
        scaffold = start = end = None
        ranges = []
        for sb in self.manifest:
            end = sb.end
            if sb.scaffold != scaffold:
                if scaffold is not None:
                    ranges.append(Range(scaffold, start, end))
                start = sb.start
                scaffold = sb.scaffold
            else:
                end = sb.end
        ranges.append(Range(scaffold, start, end))

        return ranges
        
    def get_ranges(self, end):
        first_reversed = self.ranges[0].reversed()
        if self.ordered:
            if end == 'first':
                return [first_reversed]
            elif end == 'last':
                return [self.ranges[-1]]
        else:
            if len(self.ranges) == 1:
                return [first_reversed, self.ranges[0]]
            else:
                return [first_reversed, self.ranges[-1]]
        
    def empty(self):
        self.logs = []
        self.manifest = []
        self.update()
    
    def discard(self, c_range):
        for log in self.logs:
            self.genome.blocks[log[0]][log[1]].contained = c_range
        self.empty()
        
    def replace(self, logs):
        self.empty()
        for scaffold, start, direction in logs:
            self.append(scaffold, start, direction)
        self.update()

    def reverse(self):
        newlogs = deepcopy(self.logs)
        self.empty()
        for scaffold, start, direction in reversed(newlogs):
            self.append(scaffold, start, direction*-1)

    def append(self, scaffold, start, direction=1):
        self.logs = self.logs + [(scaffold, start, direction)]
        self.manifest = self.manifest + [SummaryBlock(scaffold, start, self.genome.blocks[scaffold][start], direction)]
        self.update()

    def prepend(self, scaffold, start, direction=1):
        self.logs = [(scaffold, start, direction)] + self.logs
        self.manifest = [SummaryBlock(scaffold, start, self.genome.blocks[scaffold][start], direction)] + self.manifest
        self.update()

    def merge(self, other):
        for scaffold, start, direction in other.logs:
            self.append(scaffold, start, direction)

    def update(self):
        self.collapse()
        self._ranges = self.set_ranges()

    def check_chain(self, chain):
        
        if len(chain) <= 1:
            return True
        
        if chain[0] < chain[1]:
            direction = 1
        else:
            direction = -1
        
        for i in range(1, len(chain)-1):
            if chain[i] < chain[i+1]:
                this_dir = 1
            else:
                this_dir = -1
            if direction != this_dir:
                return False
                break
        
        return True


    def collapse_consecutive(self):
        newsummary = []
        for i in range(0, len(self.manifest)):
            sbi = self.manifest[i]
            if len(newsummary) == 0:
                newsummary.append(sbi)
            else:
                last = newsummary[-1]
                if (last.scaffold==sbi.scaffold and
                    abs(last.end-sbi.start) == 1 and
                    last.cm == sbi.cm):
                    last.end = sbi.end
                else:
                    newsummary.append(sbi)

        self.manifest = newsummary

    def collapse_trios(self):
        if len(self.manifest) < 3:
            return

        newsummary = []
        i = 0
        while i < len(self.manifest):
            sbi = self.manifest[i]
            if i < len(self.manifest)-2:
                sbj = self.manifest[i+1]
                sbk = self.manifest[i+2]
                if (sbi.scaffold == sbj.scaffold and sbi.scaffold == sbk.scaffold and
                    abs(sbi.end - sbj.start) == 1 and abs(sbj.end - sbk.start) == 1 and
                    sbi.cm == sbk.cm             and sbj.cm == -1):
                        sbi.end = sbk.end
                        newsummary.append(sbi)
                        i += 3
                        continue

            newsummary.append(sbi)
            i += 1
        
        self.manifest = newsummary
        

    def collapse(self):
        self.collapse_consecutive()
        self.collapse_trios()

    def extend(self):
        self.extend_dir(0)
        self.extend_dir(-1)

    def extend_dir(self, item):
        first_scaffold, first_start, direction = self.logs[item]
        ext_start = first_start
        extend = 0
        starts_to_extend = []
        scaffold = self.genome.blocks[first_scaffold]
        chromosome = scaffold[first_start].chromosome
        
        while True:
            if item == 0 and direction == 1 or item == -1 and direction == -1:
                ext_start = scaffold[ext_start].prev_block
            else:
                ext_start = scaffold[ext_start].next_block
        
            # Extending to ends is OK
            if ext_start == 0:
                break
        
            # If we reach another cM, abandon extension
            if scaffold[ext_start].cm != -1 or scaffold[ext_start].chromosome != '0' and scaffold[ext_start].chromosome != chromosome:
                starts_to_extend = []
                break
        
            starts_to_extend.append(ext_start)
        
        if starts_to_extend:
            if item == 0:
                for start in starts_to_extend:
                    self.prepend(first_scaffold, start, direction)
            else:
                for start in starts_to_extend:
                    self.append(first_scaffold, start, direction)

    class RaftBridge():
        def __init__(self, self_range, self_overhang, other_range, other_overhang, pb_range, pb_overhang, bridge):
            self.self_range = self_range
            self.self_overhang = self_overhang
            self.other_range = other_range
            self.other_overhang = other_overhang
            self.pb_range = pb_range
            self.pb_overhang = pb_overhang
            self.bridge = bridge
            
        def __repr__(self):
            out = '{} {:7d} = {} ({:7d} bp) = {} {:7d}'.format(
                self.self_range, self.self_overhang, self.pb_range, self.pb_overhang, self.other_range.reversed(), self.other_overhang
            )
            (self_hit, other_hit) = (self.bridge.hit1, self.bridge.hit2) if self.bridge.hit1.g.scaffold == self.self_range.scaffold else (self.bridge.hit2, self.bridge.hit1)
            out += ' {:6.2f} {:7d} {:6.2f} {:7d} {:10.2f}'.format(
                self_hit.pcid, self_hit.g.hitlen, other_hit.pcid, other_hit.g.hitlen, self.bridge.score
            )

            return out

        @property
        def overhang(self):
            return self.self_overhang + self.other_overhang
        
        @property
        def coverage(self):
            return self.bridge.hit1.g.hitlen + self.bridge.hit2.g.hitlen

    
    def add_bridge(self, self_hit, self_range, self_range_i, other_range, other_range_i, bridge):

        if self_hit == bridge.hit1:
            self_overhang, other_overhang, pb_overhang, pb_range = bridge.get_hit_overhangs(self_range, other_range)
        else:
            other_overhang, self_overhang, pb_overhang, pb_range = bridge.get_hit_overhangs(other_range, self_range)
            pb_range.reverse()

        pb_range_i = repr(pb_range)
        cur_overhang = -10000000
        if (self_range_i in self.bridges and other_range_i in self.bridges[self_range_i]
            and pb_range_i in self.bridges[self_range_i][other_range_i]):
            cur_overhang = self.bridges[self_range_i][other_range_i][pb_range_i].overhang
        
        if cur_overhang < self_overhang + other_overhang:
            self.bridges[self_range_i][other_range_i][pb_range_i] = self.RaftBridge(self_range, self_overhang,
                                                                                    other_range, other_overhang,
                                                                                    pb_range, pb_overhang, bridge)


class SummaryBlock:
    def __init__(self, scaffold, start, block, direction):
        self.scaffold = scaffold
        if direction == 1:
            self.start = start
            self.end = block.end
        elif direction == -1:
            self.start = block.end
            self.end = start
        self.cm = block.cm

    @property
    def length(self):
        return max(self.start,self.end)-min(self.start,self.end)+1

    @property
    def sequence(self):
        if self.start < self.end:
            seq = genome.sequences[self.scaffold][self.start-1:self.end]
            seq.id = '_'.join([self.scaffold, str(self.start), str(self.end)])
        else:
            seq = genome.sequences[self.scaffold][self.end-1:self.start][::-1]
            seq.id = '_'.join([self.scaffold, str(self.start), str(self.end)])
        seq.name = seq.id
        seq.description = seq.id
        return seq

        
    def __repr__(self):
        return '{}:{}-{} ({}, {} bp)'.format(self.scaffold, str(self.start), str(self.end), str(self.cm), str(self.length))


class Range():
    def __init__(self, scaffold, start, end):
        self.scaffold = scaffold
        self.start = start
        self.end = end

    def __repr__(self):
        return "{0:16s}:{1:7d}-{2:7d}".format(self.scaffold, self.start, self.end)

    def reversed(self):
        return Range(self.scaffold, self.end, self.start)
    
    def reverse(self):
        self.start, self.end = self.end, self.start


if __name__ == '__main__':
    pass
