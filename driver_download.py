#!/usr/bin/env python
# 
# driver_download.py: Script for downloading raw (JSON) OBIS files.
# 
# This is the second step in the generation of the full OBIS data set.
# It typically takes several days to run due to the response time of
# the OBIS API.
#
# Call:
# driver_download.py  (input file name and output directory are currently
#                      hard-coded)
#
# Author: Hartmut Frenzel, School of Oceanography, UW
# 
# Current revision: October 11, 2019 (major revision of API)
#
# Revision:         February 19, 2019 (output to log file)
#
# First version:    March 16, 2017

from pyobis import occurrences
import argparse
import re
import simplejson

fn_species = 'full_species_list.txt' # input file
fn_log = 'driver_download.log'

def download_species(sname, f_log, dir_json):
    size = 10000 # limit set by the OBIS API (new limit Oct 2019)

    try:
        occur = occurrences.search(scientificname = sname, size=size)
    except:
        f_log.write("Failed to download data for " + sname + "\n")
        return

    this_count = len(occur['results'])
    total_count = this_count

    file_index = 1
    fn_json_base = sname.replace(' ','_').replace('(','[').replace(')',']')
    fn_json = '{0:s}/{1:s}_{2:d}.json'.format(dir_json, fn_json_base,
                                              file_index)

    with open(fn_json, 'w') as f_obj:
        f_obj.write(simplejson.dumps(occur['results'], indent=4, 
                                     sort_keys=True))

    while this_count == size:
        # we need to download the data in stages
        # get the occurrence id  of the last entry as "offset"
        last_id = occur['results'][size-1]['id']
        success = 0
        while not success:
            try:
                occur = occurrences.search(scientificname = sname, 
                                           size=size, after=last_id)
                success = 1
                this_count = len(occur['results'])
                total_count += this_count
            except ConnectionError:
                f_log.write("connection error!")
                time.sleep(60) # wait a minute before trying again
        file_index += 1
        fn_json = '%s/%s_%d.json' % (dir_json, fn_json_base, file_index)
        with open(fn_json, 'w') as f_obj:
            f_obj.write(simplejson.dumps(occur['results'], indent=4, 
                                         sort_keys=True))

    f_log.write(str(total_count) + " data points\n")
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # option:
    parser.add_argument('-f', '--first_species', default = None,
                        help='name of first species')
    parser.add_argument('-j', '--dir_json', default = 'JSON',
                        help='name of directory for JSON files')
    args = parser.parse_args()

    # create the log file
    f_log = open(fn_log, 'w') 

    # create the Regex object
    lineRegex = re.compile(r'\d+:([\w+ ()]+)')

    if args.first_species:
        do_download = 0
    else:
        do_download = 1 # default: download everything

    line_count = 0
    count = 0
    with open(fn_species) as file_object:
        for line in file_object:
            line_count += 1
            match_obj = lineRegex.search(line)
            if (match_obj):
                count += 1
                sci_name = match_obj.group(1)
                if sci_name == 'TAXON_RANK':
                    if do_download:
                        f_log.write(str(line_count) + ': not a species!\n')
                        continue
                if not do_download:
                    if sci_name == args.first_species:
                        do_download = 1
                        
                if do_download:        
                    f_log.write(str(line_count) + ':' + sci_name + ' -- ')
                    download_species(sci_name, f_log, args.dir_json)

    f_log.write(str(line_count) + " lines processed\n")
    f_log.close()

