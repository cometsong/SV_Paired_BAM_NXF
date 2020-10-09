#!/usr/bin/env python3

"""Annotate VCF file

Usage: annot_vcf_with_exon.py annotated.intersected.bedpe original.vcf
Output: annotated-vcf to stdout; errors logged to stderr

Apply 'InExon' to original VCF file.
Based on matching annotations from VCF and annotated bedpe.

"How to add custom fields to the VCF files?" idea from: https://www.biostars.org/p/209188/
"""

__author__ = str(('Mike Lloyd', 'Benjamin Leopold'))

import sys
import pysam

# Inputs:
annot_bedpe = sys.argv[1]
orig_vcf = sys.argv[2]

# Outputs:
annot_vcf = sys.stdout
error_log = sys.stderr

bedpe_dict = {}
# bedpe_alt_dict = {} #NB: see NOTE calls missing from BEDPE

with open(annot_bedpe, 'r') as bedpe:
    for line in bedpe:
        parsed = line.rstrip().split('\t')

        # alt_key = parsed[0]+str(parsed[1])

        if '.' in parsed[4]:
            cross_ref = parsed[0]+str(parsed[1])
        else:
            cross_ref = parsed[4]

        if cross_ref in bedpe_dict \
        and 'FALSE' in bedpe_dict[cross_ref] \
        and 'exon' in parsed[8]:
            bedpe_dict[cross_ref] = 'TRUE'
            # bedpe_alt_dict[alt_key] = 'TRUE'

        elif 'exon' in parsed[8]:
            bedpe_dict[cross_ref] = 'TRUE'
            # bedpe_alt_dict[alt_key] = 'TRUE'

        else:
            bedpe_dict[cross_ref] = 'FALSE'
            # bedpe_alt_dict[alt_key] = 'FALSE'

myvcf = pysam.VariantFile(orig_vcf,'r')
myvcf.header.info.add('InEXON','1','String',
                      'Is SV call contained within an exonic region or regions: TRUE, FALSE')

print(myvcf.header, end='', file=annot_vcf)

line_counter = 0

for variant in myvcf:
    if variant.id:
        try:
            if 'TRUE' in bedpe_dict[variant.id]:
                variant.info['InEXON']='TRUE'
            else:
                variant.info['InEXON']='FALSE'

        except KeyError:
            print(f"ERROR: can’t find variant.id: '{variant.id!s}' in bedpe "
                  f"file '{annot_bedpe!s}'",
                  file=error_log)
            pass
        except Exception as e:
            raise e

        else:
            print(variant, end='', file=annot_vcf)


        """
        NOTE: The reason that calls missing from BEDPE --
        INV calls are 'reciprocal' and the IDs are '123_1' and '123_2'.
        All '*_2' calls are omitted in the BEDPE calls.

        However, not all calls have _1. Some are unique ID's, and hard to parse.
        Parsing is very hard because of this, and I was not able to find a field
        or combination of files that worked.

        My script that generate the BED files for intersection accounts
        for and checks if exon is associated with INV calls on chr1 and chr2.

        The reciprocal calls are not being printed to the final file. We will
        need to convey to Franscecio The Database will have to parse the CHR2
        and END position.

        e.g.
        chr1	8506514	INV00840SUR	.	.	441	PASS	SUPP=4;SUPP_VEC=1111;SVLEN=2006;SVTYPE=INV;SVMETHOD=SURVIVOR1.0.7;CHR2=chr1;END=8508232;CIPOS=-181,1;CIEND=0,361;STRANDS=++	GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO	0/1:NA:2198:27,0:++:89,441:INV,INV:INV00000021:NA:NA:chr1_8506395-chr1_8508593,chr1_8506506-chr1_8508427	./.:NA:1921:0,9:++:.:INV:27_1:NA:NA:chr1_8506504-chr1_8508425	1/.:NA:1718:0,0:++:.:INV::.:.:chr1_8506514-chr1_8508232	./.:NA:2188:0,0:++:.,.:INV,INV:MantaINV_811_0_1_0_0_0:NA:NA:chr1_8506333-chr1_8508521,chr1_8506515-chr1_8508469
        chr1	8508425	27_2	N	N]chr1:8506504]	.	PASS	SUPP=1;SUPP_VEC=0100;SVLEN=1921;SVTYPE=INV;SVMETHOD=SURVIVOR1.0.7;CHR2=chr1;END=8506504;CIPOS=0,0;CIEND=0,0;STRANDS=++	GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO	./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN	./.:NA:1921:0,9:++:.:INV:27_2:NA:NA:chr1_8508425-chr1_8506504	./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN	./.:NaN:0:0,0:--:NaN:NaN:NaN:NAN:NAN:NAN
        """

        # try:
        #     if 'TRUE' in bedpe_alt_dict[variant.info['CHR2']+str(variant.stop)]:
        #         variant.info['InEXON']='TRUE'
        #     else:
        #         variant.info['InEXON']='FALSE'
        # except KeyError:
        #     variant.info['InEXON']='KeyError'

    else:
        """
        NOTE: When SURVIVOR does the conversion from VCF to bedpe, it adds 1 base
        to each position, so we adjust variant.pos check for numerical match.

        e.g.
        original VCF:
        chr1	4849034	.	.	N[chr8:11432035[
        bedpe:
        chr1	4849035	4849035	chr14849035	.
        """
        try:
            chr_pos = variant.chrom + str(variant.pos + 1)
            # print(f"DEBUG-> checking chr_pos: '{chr_pos}'", file=error_log)
            if 'TRUE' in bedpe_dict[chr_pos]:
                variant.info['InEXON']='TRUE'
            else:
                variant.info['InEXON']='FALSE'
        except KeyError:
            # print(f"ERROR: can’t find variant: '{variant.chrom}' at "
            #       f"'{variant.pos!s}' in the bedpe file '{anot_bedpe!s}'. Check "
            #       "if bedpe positions are 0 or +1 adjusted relative to original "
            #       "vcf positions.”, file=error_log)
            # print(f'............ERROR: {chr_pos}')
                  #, file=error_log)
            print(f'ERROR: unable to find value: "{variant.chrom}" at postition '
                  f'"{variant.pos!s}" in the bedpe file "{anot_bedpe!s}".\n\t'
                  'Check if bedpe positions are 0 or +1 adjusted relative to '
                  'original vcf positions.', file=error_log)
            pass
        except Exception as e:
            raise e
        else:
            print(variant, end='', file=annot_vcf)
