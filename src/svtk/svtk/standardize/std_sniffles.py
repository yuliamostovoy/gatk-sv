# -*- coding: utf-8 -*-
#
"""
std_sniffles.py

Standardize Sniffles2 records.
"""
from collections import deque
from .standardize import VCFStandardizer, parse_bnd_pos, parse_bnd_strands
from svtk.utils import is_smaller_chrom


@VCFStandardizer.register('sniffles')
class SnifflesStandardizer(VCFStandardizer):
    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize Sniffles2 record.

        1) add CHR2 for all entries
        2) add END for BNDs
        3) change STRAND to STRANDS for BNDs and change + to ++ / - to --
        4) add SVLEN for BNDs
        5) for non-BNDs, make sure SVLEN is positive (originally negative for DELs)
        """

        std_rec.info['SVTYPE'] = raw_rec.info['SVTYPE']

        if std_rec.info['SVTYPE'] == 'BND':
            # ADD STRANDS
            if raw_rec.info['STRAND'] == '+':
                std_rec.info['STRANDS'] = '++'
            elif raw_rec.info['STRAND'] == '-':
                std_rec.info['STRANDS'] = '--'
            elif raw_rec.info['STRAND'] == '+-':
                std_rec.info['STRANDS'] = raw_rec.info['STRAND']
            else:
                sys.stderr.write("Error, unknown strand detected: %s\n" % raw_rec.info['STRAND'])
                sys.exit(1)

            # ADD SVLEN
            if 'SVLEN' not in raw_rec.info:
                std_rec.info['SVLEN'] = -1

            # Parse CHR2 and END
            chr2, end = parse_bnd_pos(std_rec.alts[0])
            end = end + 1 #there are end==0 so it seems to be 0-based but should be 1-based (and throws error if swapping end/pos when end==0)

            if not is_smaller_chrom(std_rec.chrom, chr2):
                # swap chr2/chrom, pos/end, and reverse strandedness
                std_rec.pos, end = end, std_rec.pos
                std_rec.chrom, chr2 = chr2, std_rec.chrom
                std_rec.info['STRANDS'] = std_rec.info['STRANDS'][::-1]
            elif std_rec.chrom==chr2 and end < std_rec.pos:
                # same chr but mates still need to be flipped
                std_rec.pos, end = end, std_rec.pos
                std_rec.info['STRANDS'] = std_rec.info['STRANDS'][::-1]
        else:
                chr2, end = raw_rec.chrom, raw_rec.stop
                std_rec.info['SVLEN'] = raw_rec.info['SVLEN']
                if std_rec.info['SVLEN'] < 0:
                    std_rec.info['SVLEN'] = std_rec.info['SVLEN'] * -1

        std_rec.info['CHR2'] = chr2
        std_rec.stop = end

        std_rec.info['ALGORITHMS'] = ['sniffles2']

        return std_rec

    def standardize_alts(self, std_rec, raw_rec):
        """
        Standardize ALT.

        When the full ref/alt sequence is specified for deletions or
        insertions, replace with N and <SVTYPE>
        """

        # Set reference to null
        stop = std_rec.stop
        std_rec.ref = 'N'
        std_rec.stop = stop

        # Format BND ALT
        #std_rec = super().standardize_alts(std_rec, raw_rec)

        # Format remainder of SVTYPES
        svtype = std_rec.info['SVTYPE']
        simple_alt = '<{0}>'.format(svtype)
        #if svtype != 'BND':
        #    std_rec.alts = (simple_alt, )
        std_rec.alts = (simple_alt, )

        return std_rec


    def standardize_format(self, std_rec, raw_rec):
        """
        rename DR and DV to RR and SR (num ref-supporting and variant-supporting reads)

        Note: self.std_sample_names order must match raw_rec.samples
        """

        source = std_rec.info['ALGORITHMS'][0]

        # Copy GT and GQ and rename SR and RR
        for sample, std_sample in zip(raw_rec.samples, self.std_sample_names):
            std_rec.samples[std_sample]['GT'] = raw_rec.samples[sample]['GT']
            if 'GQ' not in raw_rec.samples[sample] and 'PL' in raw_rec.samples[sample]: # format written by LRcaller
                std_rec.samples[std_sample]['GQ'] = sorted(list(raw_rec.samples[sample]['PL']))[1] # the GQ is the second-largest PL
                RR, SR, total = raw_rec.samples[sample]['AD']
                std_rec.samples[std_sample]['SR'] = SR
                std_rec.samples[std_sample]['RR'] = RR

            else: # original sniffles format
                std_rec.samples[std_sample]['GQ'] = raw_rec.samples[sample]['GQ']
                std_rec.samples[std_sample]['RR'] = raw_rec.samples[sample]['DR']
                std_rec.samples[std_sample]['SR'] = raw_rec.samples[sample]['DV']


        return std_rec
