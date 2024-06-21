#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from Read import Read
from Window import Window

class Improve:

    def __init__(self):
        pass

    def get_reads_info_improve_1(self):

        windowObject=Window()
        readsObject=Read()

        result_list = []
        for read in reads.reads.values():
            read.microsatellites = {ms_id: self.microsatellites[ms_id] for ms_id in read.support_microsatellites}

            self.this_ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, start=self.align_start,
                                                                      end=self.align_end).upper()
            self.this_ref_list = list(self.this_ref_str)
            sub_read_str = []
            sub_read_quals = []

            read_pos = 0
            for cigartuple in self.cigartuples:
                if cigartuple[0] in [0, 7, 8]:
                    match_read = list(self.this_read_str[read_pos:read_pos + cigartuple[1]])
                    match_quals = list(self.this_read_quals[read_pos:read_pos + cigartuple[1]])
                    sub_read_str.extend(match_read)
                    sub_read_quals.extend(match_quals)
                    read_pos += cigartuple[1]
                elif cigartuple[0] in [1, 4, 5]:
                    if cigartuple[0] == 1:
                        if len(sub_read_str) == 0:
                            continue
                        else:
                            sub_read_str[-1] += self.this_read_str[
                                                read_pos:read_pos + cigartuple[1]]
                            if self.tech == "contig": continue
                            sub_read_quals[-1] += self.this_read_quals[read_pos:read_pos + cigartuple[1]]
                    elif cigartuple[0] == 5:
                        continue
                    read_pos += cigartuple[1]
                elif cigartuple[0] in [2, 3]:
                    sub_read_str.extend([""] * cigartuple[1])
                    if self.tech == "contig": continue
                    sub_read_quals.extend([""] * cigartuple[1])
                else:
                    return -1
            self.this_read_list = sub_read_str
            self.this_quals_list = sub_read_quals
            self.this_read_quals = ""

            self.microsatellites_num = len(self.microsatellites)
            read_muts = {}
            repeat_lengths = {}

            for ms_id, ms in self.microsatellites.items():
                ms_start = ms.start
                ms_end = ms.end
                ms_start_pre = ms.start_pre
                ms_end_suf = ms.end_suf

                if compar_software == 'minimap2' or compar_software == 'winnowmap2' or compar_software == 'ngmlr':
                    ms_motif = ms.repeat_unit
                    str_varia = "".join(self.this_read_list[
                                ms_start - 1 - self.align_start:ms_end - self.align_start +1])

                    if len(str_varia) >= 2:
                        if str_varia[0] != ms_motif[0]:
                            str_varia = str_varia[1:]
                        if str_varia[-1] != ms_motif[-1]:
                            str_varia = str_varia[:-1]

                    query_repeat_length=len(str_varia)

                elif compar_software == 'lra':
                    query_repeat_length = len(
                        "".join(self.this_read_list[ms_start + 3 - self.align_start:ms_end - self.align_start+0]))
                else:
                    print('compar_software is wrong!')

                prefix = self.this_read_list[ms_start_pre - self.align_start:ms_start - self.align_start]
                suffix = self.this_read_list[ms_end - self.align_start:ms_end_suf - self.align_start]
                ms_content = self.this_read_list[ms_start - self.align_start:ms_end - self.align_start]
                mismatches = []
                deletions = []
                insertions = []
                pos_based_info = {}
                pos_deletion = []
                ref_pos = ms_start_pre - 2
                for pot in range(ms_start_pre - self.align_start - 1, ms_end_suf - self.align_start):
                    if pot >= len(self.this_read_list):
                        print('Read.py:pot',pot,' len',len(self.this_read_list))
                        continue

                    ref_pos += 1
                    this_read_base = self.this_read_list[pot]
                    this_ref_base = self.this_ref_list[pot]
                    if this_read_base == this_ref_base:
                        continue
                    else:
                        if ref_pos < ms_start_pre:
                            band = [1]
                        elif ref_pos < ms_start:
                            band = [2]
                        elif ref_pos < ms_end:
                            band = [3]
                        elif ref_pos < ms_end_suf:
                            band = [4]
                        else:
                            band = [5]
                        this_read_base_len = len(this_read_base)
                        this_ref_base_len = len(this_ref_base)
                        if this_read_base_len == this_ref_base_len:
                            mismatches.append([ref_pos, this_read_base, band])
                            pos_based_info[ref_pos] = ["SNV", this_ref_base, this_read_base, band]
                        else:
                            if this_read_base_len < this_ref_base_len:
                                deletions.append([ref_pos, this_read_base, band])
                                pos_deletion.append([ref_pos, this_ref_base, band[0]])
                            else:
                                if ref_pos == ms_start - 1 and this_read_base_len > ms.repeat_unit_len:
                                    band = [3]
                                insertions.append([ref_pos, this_read_base, band])
                                if band[0] != 3:
                                    pos_based_info[ref_pos] = ["INS", "", this_read_base[1:], band]

                if len(pos_deletion) > 0:
                    deletion_start = pos_deletion[0][0]
                    deletion_len = 1
                    del_str = pos_deletion[0][1]
                    band = [pos_deletion[0][2]]
                    for i in range(1, len(pos_deletion)):
                        if pos_deletion[i - 1][0] + 1 == pos_deletion[i][0]:
                            deletion_len += 1
                            del_str += pos_deletion[i][1]
                            band.append(pos_deletion[i][2])
                        else:
                            pos_based_info[deletion_start] = ["DEL", del_str, "", set(band)]

                            deletion_len = 1
                            deletion_start = pos_deletion[i][0]
                            del_str = pos_deletion[i][1]
                            band = [pos_deletion[i][2]]
                    if len(set(band)) > 1 or list(set(band))[0] != 3:
                        pos_based_info[deletion_start] = ["DEL", del_str, "", set(band)]
                read_str = self.this_read_list[ms_start_pre - 1 - self.align_start:ms_end_suf + 2 - self.align_start]
                read_muts[ms_id] = Read_Mutation(self.read_id, repeat_length=query_repeat_length, strand=self.strand,
                                                 hap=self.hap,
                                                 mismatches=mismatches, deletions=deletions, insertions=insertions,
                                                 prefix=prefix, suffix=suffix, ms_content=ms_content,
                                                 pos_based_info=pos_based_info,
                                                 read_str=read_str
                                                 )
                repeat_lengths[ms_id] = query_repeat_length

            self.repeat_lengths = repeat_lengths
            self.mut_info = read_muts

            result_list.append(read)
        self.reads = {read.read_id: read for read in result_list}



