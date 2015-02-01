#!/usr/bin/env python3

class Raft:
    def __init__(self, bl_id, scaffold, start, chromosome):
        self.chromosome = chromosome
        self.genome = chromosome.genome
        self.id = bl_id
        self.logs = []
        self.manifest = []
        self.append(scaffold, start, 1)

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
            print("To fix:", marker_chain, "\n", self, "\n")

        return tuple(marker_chain)

    @property
    def ordered(self):
        return len(self.marker_chain) > 1


    def empty(self):
        self.logs = []
        self.manifest = []
        self.update()
        
    def replace(self, logs):
        self.empty()
        for scaffold, start, direction in logs:
            self.append(scaffold, start, direction)
        self.update()

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

    def check_chain(self, chain):
        
        if len(chain) == 1:
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





if __name__ == '__main__':
    pass
