#!/usr/bin/env python3

import sys
import argparse
import sqlite3 as db
from os.path import isfile
import operator


class Overlap:
    def __init__(self, start, end, hitlen, pccov):
        self.start = int(start)
        self.end = int(end)
        self.hitlen = int(hitlen)
        self.pccov = float(pccov)
        self.range = (self.start, self.end) if self.start < self.end else (self.end, self.start)

class Hit:
    def __init__(self, rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rpccov, qpccov, rscaffold, qscaffold, hittype):
        self.rscaffold = rscaffold
        self.qscaffold = qscaffold
        self.pcid = float(pcid)
        self.hittype = hittype

        self.r = Overlap(rstart, rend, rhitlen, rpccov)
        self.q = Overlap(qstart, qend, qhitlen, qpccov)

    def __repr__(self):
        return '{} {}-{} = {} {}-{} {}'.format(self.rscaffold, self.r.start, self.r.end, self.qscaffold, self.q.start, self.q.end, self.hittype)
        

def add_hit(hits, hit):
    if hit.rscaffold not in hits:
        hits[hit.rscaffold] = {}
    if hit.qscaffold not in hits[hit.rscaffold]:
        hits[hit.rscaffold][hit.qscaffold] = []
    hits[hit.rscaffold][hit.qscaffold].append(hit)

def collapse(lst):
    lst.sort()
    to_remove = []
    for i, hit_i in enumerate(lst):
        for j, hit_j in enumerate(lst):
            if i >= j:
                continue
            if hit_i[0] <= hit_j[0] and hit_i[1] >= hit_j[0]:
                i_end = max(hit_i[1], hit_j[1])
                lst[i] = (hit_i[0], i_end)
                to_remove.append(j)
    for i in reversed(sorted(set(to_remove))):
        del(lst[i])
    return lst

def get_prop(ranges, scaffold_length):
    ranges = collapse(ranges)
    coverage = sum([r[1]-r[0]+1 for r in ranges])
    return coverage/scaffold_length * 100, ranges

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--overlaps', type=str, required=True)
parser.add_argument('-p', '--pacbio', type=str, required=True)
parser.add_argument('-d', '--database', type=str, required=True)
args = parser.parse_args()

def load_nucmer_coords (name, file, database):
    conn = db.connect(database)
    c = conn.cursor()
    
    c.execute('drop table if exists {}'.format(name))
    
    c.execute('''CREATE TABLE {}
                 (rstart integer, rend integer, qstart integer, qend integer,
                  rhitlen integer, qhitlen integer, pcid real, rseqlen integer, qseqlen integer,
                  rpccov real, qpccov real, rscaffold text, qscaffold text, hittype text)'''.format(name))
    

    if isfile(file):
        with open(file, 'r') as o:
            for line in o:
                f = line.rstrip().split('\t')
                if len(f) < 13:
                    continue
                rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rseqlen, qseqlen, rpccov, qpccov, rscaffold, qscaffold, *rem = f
                
                if rscaffold == qscaffold:
                    continue
    
                hittype = ''
                if len(rem) == 1:
                    hittype = rem[0]
                
                c.execute("INSERT INTO {} VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)".format(name),
                          [rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rseqlen, qseqlen, rpccov, qpccov, rscaffold, qscaffold, hittype])
                
                if 'HE' in rscaffold and 'sch' in qscaffold:
                    if hittype == '[CONTAINS]':
                        hittype = '[CONTAINED]'
                    elif hittype == '[CONTAINED]':
                        hittype = '[CONTAINS]'
    
                    c.execute("INSERT INTO {} VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)".format(name),
                              [qstart, qend, rstart, rend, qhitlen, rhitlen, pcid, qseqlen, rseqlen, qpccov, rpccov, qscaffold, rscaffold, hittype])
    
    c.execute("CREATE INDEX {}_scaffolds on {} (rscaffold, qscaffold)".format(name, name))

load_nucmer_coords('pacbio_overlaps', args.pacbio, args.database)
load_nucmer_coords('genome_overlaps', args.overlaps, args.database)
sys.exit()

hits = {}
scaffold_length = {}
contained = {}
if isfile(args.overlaps):
    with open(args.overlaps, 'r') as o:
        for line in o:
            f = line.rstrip().split('\t')
            if len(f) < 13:
                continue
            rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rseqlen, qseqlen, rpccov, qpccov, rscaffold, qscaffold, *rem = f
            
            if rscaffold not in scaffold_length:
                scaffold_length[rscaffold] = int(rseqlen)
            if qscaffold not in scaffold_length:
                scaffold_length[qscaffold] = int(qseqlen)
            
            if rscaffold == qscaffold:
                continue

            hittype = ''
            if len(rem) == 1:
                hittype = rem[0]
            
            add_hit(hits, Hit(rstart, rend, qstart, qend, rhitlen, qhitlen, pcid, rpccov, qpccov, rscaffold, qscaffold, hittype))
            
            if hittype == '[CONTAINS]':
                contained[qscaffold] = rscaffold

            if 'HE' in rscaffold and 'sch' in qscaffold:
                if hittype == '[CONTAINS]':
                    hittype = '[CONTAINED]'
                elif hittype == '[CONTAINED]':
                    hittype = '[CONTAINS]'
                    contained[rscaffold] = qscaffold

                add_hit(hits, Hit(qstart, qend, rstart, rend, qhitlen, rhitlen, pcid, qpccov, rpccov, qscaffold, rscaffold, hittype))

for rscaffold in sorted(hits):
    if rscaffold in contained:
        continue
    for qscaffold in sorted(hits[rscaffold]):
        if qscaffold in contained:
            continue
        rprop, rranges = get_prop([hit.r.range for hit in hits[rscaffold][qscaffold]], scaffold_length[rscaffold])
        qprop, qranges = get_prop([hit.q.range for hit in hits[rscaffold][qscaffold]], scaffold_length[qscaffold])
        if qranges[0][0] == 1 and qranges[-1][1] == scaffold_length[qscaffold]:
            print(rscaffold, qscaffold)
            for hit in hits[rscaffold][qscaffold]:
                print(hit)
            print(rscaffold, scaffold_length[rscaffold], rprop, rranges)
            print(qscaffold, scaffold_length[qscaffold], qprop, qranges)
            print()

for rscaffold in hits:
    for qscaffold in hits[rscaffold]:
        for hit in hits[rscaffold][qscaffold]:
            c.execute("INSERT INTO overlaps VALUES (?,?,?,?,?,?,?,?,?,?,?,?)",
                      [hit.r.start, hit.r.end, hit.q.start, hit.q.end, hit.r.hitlen, hit.q.hitlen, hit.pcid, hit.r.pccov, hit.q.pccov, hit.rscaffold, hit.qscaffold, hit.hittype])

c.execute("CREATE INDEX scaffolds on overlaps (rscaffold, qscaffold)")