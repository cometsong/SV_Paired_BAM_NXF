
#!/usr/bin/env python3
""" This script is correcting TRA calls from single to multi-line.
    This is done to facilitate exon searching / matching via bedtools intersect.
"""
### Usage: <script.py> inFile outFile

import sys

with open(sys.argv[1], 'r') as inFile:
    with open(sys.argv[2], 'w') as output:

        for line in inFile:
            parsed = line.rstrip().split('\t')

            parsed_01 = parsed[0]+parsed[1]
            if (parsed[0] == parsed[3]
            and parsed[1] == parsed[2]
            and parsed[4] == parsed[5]):
                output.write( '\t'.join([parsed[0],parsed[1],parsed[4],parsed_01,parsed[6]]) +'\n' )

            elif (parsed[0] != parsed[3]):
                output.write( '\t'.join([parsed[0],parsed[1],parsed[2],parsed_01,parsed[6]]) +'\n' )
                output.write( '\t'.join([parsed[3],parsed[4],parsed[5],parsed_01,parsed[6]]) +'\n' )

            else:
                print('FOUND_MIXED')
                output.write('\t'.join([line.rstrip(),'mixed_start_end',parsed_01,parsed[4]]) +'\n')
