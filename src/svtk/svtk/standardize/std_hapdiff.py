# -*- coding: utf-8 -*-
#
"""
std_hapdiff.py

Standardize hapdiff records.
"""
import sys
from collections import deque
from .standardize import VCFStandardizer, parse_bnd_pos, parse_bnd_strands
from svtk.utils import is_smaller_chrom


@VCFStandardizer.register('hapdiff')
class HapdiffStandardizer(VCFStandardizer):
    def standardize_records(self):
        """
        Filter hapdiff VCF

        Mated events are not marked with SECONDARY tag; must filter manually.
        """

        partner_chr_pos_prev = deque()
        self_chr_pos_prev = deque()
        for record in self.filter_raw_vcf():
            # Filter mates of BND calls
            if record.info['SVTYPE'] == "BND":
                self_chr_pos = '%s_%d' % (record.chrom, record.pos)
                partner_chr_pos = '%s_%d' % (parse_bnd_pos(record.alts[0]))
                if self_chr_pos in partner_chr_pos_prev and partner_chr_pos in self_chr_pos_prev:
                    continue

                # Track info of observed records
                self_chr_pos_prev.append(self_chr_pos)
                partner_chr_pos_prev.append(partner_chr_pos)

            yield self.standardize_record(record)

    def standardize_info(self, std_rec, raw_rec):
        """
        Standardize hapdiff record.

        1) add CHR2 for all entries
        2) add END for BNDs
        3) add STRANDS for BNDs
        4) add SVLEN for BNDs and INVs
        5) for non-BNDs, make sure SVLEN is positive (originally negative for DELs)
        """

        std_rec.info['SVTYPE'] = raw_rec.info['SVTYPE']
        if raw_rec.info['SVTYPE'] == 'DUP:TANDEM':
            std_rec.info['SVTYPE'] = 'DUP'

        if std_rec.info['SVTYPE'] == 'INV':
            std_rec.info['SVLEN'] = abs(std_rec.stop-std_rec.pos)
            chr2, end = raw_rec.chrom, raw_rec.stop
        elif std_rec.info['SVTYPE'] == 'BND':
            # ADD STRANDS
            std_rec.info['STRANDS'] = parse_bnd_strands(raw_rec.alts[0])

            # ADD SVLEN
            std_rec.info['SVLEN'] = -1

            # Parse CHR2 and END
            chr2, end = parse_bnd_pos(std_rec.alts[0])

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

        std_rec.info['ALGORITHMS'] = ['hapdiff']

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
        Add GT if missing
        Add GQ

        Note: self.std_sample_names order must match raw_rec.samples
        """

        source = std_rec.info['ALGORITHMS'][0]
        null_GTs = [(0, 0), (None, None), (0, ), (None, )]

        # Add per-sample genotypes and copy CN if it exists
        for sample, std_sample in zip(raw_rec.samples, self.std_sample_names):
            if 'GT' not in raw_rec.samples[sample]:
               raw_rec.samples[sample]['GT'] = (None, None)
            gt = raw_rec.samples[sample]['GT']
            if self.call_null_sites:
                if gt == (None, None):
                    gt = (0, 1)
                if gt == (None,):
                    gt = (1,)
            std_rec.samples[std_sample]['GT'] = gt

#            sys.stderr.write("%s\t%s\t%s\n" % (sample, raw_rec.id, raw_rec.samples[sample]['AD']))
            std_rec.samples[std_sample][source] = 1
            std_rec.samples[std_sample]['GQ'] = 1

        return std_rec
