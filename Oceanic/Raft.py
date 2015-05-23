#!/usr/bin/env python3

from collections import defaultdict
from copy import deepcopy

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

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
        output = ''
        if self.marker_chain.chain is None:
            output += 'To fix:\n'
        if self.scaffold in self.genome.offcuts:
            output += 'Offcut to {}\n'.format(','.join(self.genome.offcuts[self.scaffold]))
        if self.name in self.genome.haplotypes:
            output += 'Haplotype {}\n'.format(self.genome.haplotypes[self.name])
        output += '\n'.join([repr(m) for m in self.manifest]) + '\n'
        output += 'Length: {}\n'.format(self.length)
        return output

    def __hash__(self):
        return self.id

    def __eq__(self, other):
        return self.id == other.id

    def __iter__(self):
        return iter(self.logs)
    
    def __len__(self):
        return len(self.logs)
    
    def __contains__(self, var):
        for log in self.logs:
            if log[0] == var or log[1] == var:
                return True
        return False

    @property
    def offcuts(self):
        if self.scaffold in self.genome.offcuts:
            return self.genome.offcuts[self.scaffold]
        return None

    @property
    def haplotype(self):
        if self.name in self.genome.haplotypes:
            return self.genome.haplotypes[self.name]
        else:
            return None

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
    def mappedlength(self):
        length = 0
        for m in self.manifest:
            if m.cm != -1:
                length += m.length
        return length

    @property
    def sequence(self):
        sequence = sum([h.sequence for h in self.hooks], Seq("", generic_dna))
        sequence.id = self.name
        sequence.description = self.name
        return sequence

    @property
    def marker_chain(self):
        mc = MarkerChain(self.manifest)
        return mc

    @property
    def ordered(self):
        return len(self.marker_chain) > 1

    @property
    def hooks(self):
        if not self._hooks:
            self._hooks = self.forge_hooks()
        return self._hooks

    def empty(self):
        self.logs = []
        self.manifest = []
        self.update()

    def discard(self):
        self.genome.refuse.append(self.summary())
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

    def get_log_start(self, pos):
        for log in self.logs:
            if pos >= log[1] and pos <= self.genome.blocks[log[0]][log[1]].end:
                return log[1]
        else:
            return None

    def merge(self, other):
        for scaffold, start, direction in other.logs:
            self.append(scaffold, start, direction)

    def trim(self, hook, overhang):
        pass

    def update(self):
        self.remove_duplicates()
        self.collapse()
        self._hooks = self.forge_hooks()

    def remove_duplicates(self):
        dupes = defaultdict(int)
        new_logs = []
        for log in self.logs:
            dupes[log] += 1
            if dupes[log] <= 1:
                new_logs.append(log)
        self.logs = new_logs
        

        dupes = defaultdict(int)
        new_manifest = []
        for sb in self.manifest:
            sbkey = '{}_{}'.format(sb.scaffold, sb.start)
            dupes[sbkey] += 1
            if dupes[sbkey] <= 1:
                new_manifest.append(sb)
        self.manifest = new_manifest

        
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
        self.update()

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

    def write(self):
        if self.genome.revised_fasta:
            SeqIO.write(self.sequence, self.genome.revised_fasta, "fasta")
            for sb in self.manifest:
                self.chromosome.revised_db.execute("insert into scaffold_map values (?,?,?,?,?,?)",
                      [self.chromosome.name, sb.cm, sb.scaffold, min(sb.start, sb.end), max(sb.start, sb.end), sb.length])
        
        return self.summary()
    
    def summary(self):
        scaffolds = []
        for scaffold, start, direction in self.logs:
            scaffolds.append(self.genome.blocks[scaffold][start])
        return scaffolds


    def forge_hooks(self):
        scaffold = start = end = None
        hooks = []
        for sb in self.manifest:
            if sb.scaffold != scaffold:
                if scaffold is not None:
                    hooks.append(Hook(self, scaffold, start, end))
                start = sb.start
                scaffold = sb.scaffold
            end = sb.end
        hooks.append(Hook(self, scaffold, start, end))

        return hooks

    def turn_hooks(self, end):
        for hook in self.hooks:
            hook.faces = []
        
        if self.ordered:
            if end == 'first':
                self.hooks[0].faces = [(self.hooks[0].end, self.hooks[0].start)]
            elif end == 'last':
                self.hooks[-1].faces = [(self.hooks[-1].start, self.hooks[-1].end)]
        else:
            if len(self.hooks) == 1:
                self.hooks[0].faces = [(self.hooks[0].start, self.hooks[0].end), (self.hooks[0].end, self.hooks[0].start)]
            else:
                self.hooks[0].faces = [(self.hooks[0].end, self.hooks[0].start)]
                self.hooks[-1].faces = [(self.hooks[-1].start, self.hooks[-1].end)]

    def hook_knots(self, end, rope_name):
        knots = {}
        for hook in self.hooks:
            if hook.faces:
                for knot in hook.knots:
                    if knot.rope.lookup[knot].scaffold == rope_name:
                        knots[knot] = 1
        return knots

    def tie(self, hit, rope):
        for hook in self.hooks:
            knot = hook.tie(hit, rope)
            if knot:
                return knot
        else:
            return None

    def get_springline(self, self_end, other):
        other_end = 'last' if self_end == 'first' else 'first'
        springline = None
        for hook in self.hooks:
            if hook.faces:
                springline = hook.get_springline(other, other_end)
        return springline



    class RaftBridge():
        def __init__(self, other, self_range, self_overhang, other_range, other_overhang, pb_range, pb_overhang, bridge):
            self.self_range = self_range
            self.self_overhang = self_overhang
            self.other = other
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

    
    def add_bridge(self, other, self_hit, self_range, other_range, bridge):

        if self_hit == bridge.hit1:
            self_overhang, other_overhang, pb_overhang, pb_range = bridge.get_hit_overhangs(self_range, other_range)
        else:
            other_overhang, self_overhang, pb_overhang, pb_range = bridge.get_hit_overhangs(other_range, self_range)
            pb_range.reverse()

        cur_overhang = -10000000
        if (self_range in self.bridges and other_range in self.bridges[self_range]
            and pb_range in self.bridges[self_range][other_range]):
            cur_overhang = self.bridges[self_range][other_range][pb_range].overhang
        
        if cur_overhang < self_overhang + other_overhang:
            self.bridges[self_range][other_range][pb_range] = self.RaftBridge(other, self_range, self_overhang,
                                                                              other_range, other_overhang,
                                                                              pb_range, pb_overhang, bridge)



class MarkerChain:
    def __init__(self, manifest=None, chain=None):
        self.chain = []
        if manifest:
            for m in manifest:
                if m.cm != -1:
                    self.chain.append(m.cm)
                    if len(self.chain)>1 and self.chain[-2] == self.chain[-1]:
                        del self.chain[-1]
        elif chain:
            self.chain = chain

        if not self.check():
            self.chain = None

    def __getitem__(self, i):
        return self.chain[i]

    def __contains__(self, key):
        for cm in self.chain:
            if cm == key:
                return True
        else:
            return False

    def __repr__(self):
        return repr(self.chain)

    def __len__(self):
        if self.chain == None:
            return 0
        else:
            return len(self.chain)

    def __add__(self, other):
        new = self.chain + other.chain
        collapsed = [new[0]]
        for i in range(1,len(new)):
            if new[i] != collapsed[-1]:
                collapsed.append(new[i])
        return MarkerChain(chain=collapsed)

    def check(self):
        
        if len(self.chain) <= 1:
            return True
        
        if self.chain[0] < self.chain[1]:
            direction = 1
        else:
            direction = -1
        
        for i in range(1, len(self.chain)-1):
            if self.chain[i] < self.chain[i+1]:
                this_dir = 1
            else:
                this_dir = -1
            if direction != this_dir:
                return False
        
        return True



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

        
    def __repr__(self):
        return '{}:{}-{} ({}, {} bp)'.format(self.scaffold, str(self.start), str(self.end), str(self.cm), str(self.length))


class Hook():
    def __init__(self, raft, scaffold, start, end, knots=None):
        self.raft = raft
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.faces = []
        self.knots = {} if knots is None else knots

    def __repr__(self):
        return "{0:16s}:{1:7d}-{2:7d} {3}".format(self.scaffold, self.start, self.end, self.faces)

    def __lt__(self, other):
        if self.scaffold < other.scaffold:
            return self
        elif self.start < other.start:
            return self
        else:
            return other

    def __iter__(self):
        return iter(self.knots)

    @property
    def length(self):
        if self.start < self.end:
            return self.end - self.start + 1
        else:
            return self.start - self.end + 1
    
    @property
    def sequence(self):
        if 'Gap' in self.scaffold:
            seq = 'N' * self.length
        elif self.start < self.end:
            seq = self.raft.genome.sequences[self.scaffold][self.start-1:self.end]
        else:
            seq = self.raft.genome.sequences[self.scaffold][self.end-1:self.start].reverse_complement()
        return seq

    @property
    def ropes(self):
        return [k.rope for k in self.knots]

    def order(self, start, end):
        if start < end:
            return start, end
        else:
            return end, start


    def tie(self, hit, rope):
        if hit.scaffold != self.scaffold:
            return None
        h_start, h_end = self.order(self.start, self.end)
        t_start, t_end = self.order(hit.start, hit.end)
        if h_start <= t_start and h_end >= t_end:
            new_knot = RaftKnot(hit.scaffold, hit.start, hit.end, self, rope)
            self.knots[new_knot] = True
            return new_knot
        else:
            return None
    
    def untie(self, knot):
        del self.knots[knot]
    
    def get_springline(self, other, other_end):
        print(self)
        springline = None
        for knot in self.knots:
            print(knot)
            springline = knot.get_springline(other, other_end)
        return springline
        
class Knot():
    def __init__(self, scaffold, start, end):
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.length = max(self.start, self.end) - min(self.start, self.end) + 1
    
    def __repr__(self):
        return '{}\t{}\t{}\t{}'.format(self.scaffold, self.start, self.end, self.length)

class RaftKnot(Knot):
    def __init__(self, scaffold, start, end, hook, rope):
        self.hook = hook
        self.rope = rope
        super().__init__(scaffold, start, end)

    def __repr__(self):
        out = super().__repr__()
        if self in self.rope.lookup:
            out += '\ttied to {}'.format(self.rope.lookup[self])
        return out

    def untie(self):
        self.hook.untie(self)
    
    def get_springline(self, other, other_end):
        springlines = self.rope.get_springlines(self, other, other_end)
        for springline in springlines:
            print(springline)
        return springlines[0] if springlines else None
        

if __name__ == '__main__':
    pass
