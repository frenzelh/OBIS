#!/usr/bin/env python
#
# generate_species_list.py: Script for downloading a list of available 
# species and other taxon ranks from OBIS. By default a large range of
# possible taxon IDs (1 to 1500000) is scanned, which takes over a week
# to run due to the response time of the OBIS API.
#
# This is the first (and slowest) step in the generation of full OBIS data sets, 
# which may be skipped if it can be assumed that no relevant new species
# were added to the OBIS database.
#
# Call (with default options):
# ./generate_species_list.py > full_species_list.txt 2> species_list_errors.txt
#
# Use option -h to show a list of options.
#
# Author: Hartmut Frenzel, School of Oceanography, UW
# 
# First version: November 29, 2016

import argparse
import sys
from pyobis import taxa

def get_species_by_taxon_id(id):
    '''Try to get the name of the species with the given ID.'''
    try:
        res = taxa.taxon(id)
        if res['results'][0]['taxonID'] == id:
            if 'species' in res['results'][0]:
                print('{0:d}:{1:s}'.format(id, res['results'][0]['species']))
            elif 'taxonRank' in res['results'][0]:
                print('{0:d}:TAXON_RANK={1:s}'.format(id,
                                                      res['results'][0]['taxonRank']))
            else:
                print('{0:d}'.format(id))
        else:
            print('Taxon number {0:d} does not exist (wrong ID).'.format(id),
                  file=sys.stderr)
    except:
        print('Taxon number {0:d} does not exist.'.format(id), file=sys.stderr)
  
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # options:
    parser.add_argument("-f", "--first", default=1, type=int,
                        help="number of the first taxon scanned")
    parser.add_argument("-l", "--last", default=1500000, type=int,
                        help="number of the last taxon scanned")
    parser.add_argument("-i", "--input", default=None, type=str,
                        help="name of input file with taxon IDs")
    args = parser.parse_args()

    if args.input:
        print('reading taxon IDs from {0:s}'.format(args.input))
    else:
        for id in range(args.first,args.last+1):
            get_species_by_taxon_id(id)
