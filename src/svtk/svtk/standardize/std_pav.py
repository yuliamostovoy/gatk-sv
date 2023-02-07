# -*- coding: utf-8 -*-
#
"""
std_pav.py

Standardize PAV records.
"""
import sys
from collections import deque
from .standardize import VCFStandardizer, parse_bnd_pos, parse_bnd_strands
from svtk.utils import is_smaller_chrom

@VCFStandardizer.register('pav')
class PAVStandardizer(VCFStandardizer):
    def standardize_records(self):
        """
        Filter PAV VCF

        remove SNVs
        """

        for record in self.filter_raw_vcf():
            if record.info['SVTYPE']=="SNV":
                continue

            yield self.standardize_record(record)

    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize PAV record.

        1) add CHR2 for all entries
        5) for non-BNDs, make sure SVLEN is positive (originally negative for DELs)
        """

        std_rec.info['SVTYPE'] = raw_rec.info['SVTYPE']

        chr2, end = raw_rec.chrom, raw_rec.stop
        std_rec.info['SVLEN'] = raw_rec.info['SVLEN']
        if std_rec.info['SVLEN'] < 0:
            std_rec.info['SVLEN'] = std_rec.info['SVLEN'] * -1

        std_rec.info['CHR2'] = chr2
        std_rec.stop = end

        std_rec.info['ALGORITHMS'] = ['pav']

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

        # Format SVTYPES
        svtype = std_rec.info['SVTYPE']
        simple_alt = '<{0}>'.format(svtype)
        std_rec.alts = (simple_alt, )

        return std_rec


    def standardize_format(self, std_rec, raw_rec):
        """
        Add GT if missing
        Add GQ

        Note: self.std_sample_names order must match raw_rec.samples
        """

        source = std_rec.info['ALGORITHMS'][0]
        null_GTs = [(0, 0), (None, None), (0, ), (None, )]

        # Add per-sample genotypes and copy CN, AD, and DP if they exist
        for sample, std_sample in zip(raw_rec.samples, self.std_sample_names):
            if 'GT' not in raw_rec.samples[sample]:
               raw_rec.samples[sample]['GT'] = (None, None)
            gt = raw_rec.samples[sample]['GT']
            if self.call_null_sites:
                if gt == (None, None):
                    gt = (0, 1)
                if gt == (None,):
                    gt = (1,)

            # standardizing GTs
            gt = list(gt)
            if 1 in gt: # if GT includes a 1, make sure it's not "./1" or "1/." because that causes issues downstream
                for x in [0,1]:
                    if gt[x]==None:
                        gt[x] = 0
            if None not in gt: # don't try to sort [None, None]
                gt.sort() # so (1,0) becomes (0,1)
            gt = tuple(gt)

            std_rec.samples[std_sample]['GT'] = gt

            std_rec.samples[std_sample][source] = 1
            std_rec.samples[std_sample]['GQ'] = 20

        return std_rec
