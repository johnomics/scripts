#!/usr/bin/python3 -u

import os
import sys
import argparse
import sqlite3 as sql
import string
from collections import defaultdict

def open_input_database(database):
    try:
        if not os.path.exists(database):
            raise IOError
        conn = sql.connect(database)
        cursor = conn.cursor()
        return conn, cursor
    except IOError:
        print("Can't open database {}".format(database))
        sys.exit(1)
    except sql.Error:
        print("SQL error")
        sys.exit(1)


def get_args():
    parser = argparse.ArgumentParser(description='''Count recombinations per individual per chromosome

        -d database
        ''')

    parser.add_argument('-d', '--database', type=str, required=True)
    return parser.parse_args()

if __name__ == '__main__':
    
    args = get_args()

    conn_in, ci = open_input_database(args.database)

    query = 'select chromosome, cm, clean, sum(length) from chromosome_map group by chromosome, cm order by chromosome, cm'
    rperi = defaultdict(list)
    previous = None
    for chromosome, cm, clean, length in ci.execute(query):
        marker = ''.join([x for x in clean if x in ['A','B','H']])
        if chromosome not in rperi:
            rperi[chromosome] = [0] * len(marker)
            previous = None
        if previous:
            for i, call in enumerate(marker):
                if previous[i] != call:
                    rperi[chromosome][i] += 1
        previous = marker

    total_rperi = [0] * len(marker)
    for chromosome in sorted(rperi):
        print(chromosome, rperi[chromosome])
        for i, rec in enumerate(rperi[chromosome]):
            total_rperi[i] += rec
            
    print('TOTAL', total_rperi)
    
    conn_in.close()
