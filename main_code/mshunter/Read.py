#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from units import *
import pysam
import numpy as np
from collections import *
import math

class Read:

    def __init__(self, read_id, alignment, reference, chrom, tech=""):
        self.chrom = chrom
        self.read_name = alignment.query_name
        self.align_start = alignment.reference_start
        self.align_end = alignment.reference_end
        self.this_read_str = alignment.query_sequence.upper()
        self.tech = tech
        if tech == "contig":
            self.this_read_quals = []
        else:
            self.this_read_quals = "".join([chr(i + 33) for i in alignment.query_qualities])
        self.strand = False if alignment.is_reverse else True
        self.this_read_list = []
        self.this_quals_list = []
        self.this_ref_str = ""
        self.this_ref_list = []
        self.read_id = read_id
        self.reference = reference
        self.support_microsatellites = []
        self.microsatellites = {}
        self.direction = False if alignment.is_reverse else True
        self.hap = int(alignment.get_tag("HP")) if alignment.has_tag("HP") else 0
        self.cigartuples = alignment.cigartuples
        self.mut_info = {}
        self.repeat_lengths = {}

        self.cigars_list=[]

        self.repeat_lengths_GmmHmm = {}

        self.hmmRes_draw={}

    def get_read_str(self):

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
                        sub_read_str[-1] += self.this_read_str[read_pos:read_pos + cigartuple[1]]

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

    def get_read_str_2(self):

        self.this_ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, start=self.align_start,
                                                                  end=self.align_end).upper()
        self.this_ref_list = list(self.this_ref_str)
        sub_read_str = []
        sub_read_quals = []

        read_pos = 0

        cigars_list=[]
        for cigartuple in self.cigartuples:

            if cigartuple[0] in [0, 7, 8]:

                match_read = list(self.this_read_str[read_pos:read_pos + cigartuple[1]])
                match_quals = list(self.this_read_quals[read_pos:read_pos + cigartuple[1]])
                sub_read_str.extend(match_read)

                sub_read_quals.extend(match_quals)
                read_pos += cigartuple[1]

                cigars_list.extend(['M'] *cigartuple[1])

            elif cigartuple[0] in [1, 4, 5]:

                if cigartuple[0] == 1:
                    if len(sub_read_str) == 0:
                        continue
                    else:
                        sub_read_str[-1] += self.this_read_str[read_pos:read_pos + cigartuple[1]]

                        cigars_list[-1]+=''.join(['I']*cigartuple[1])

                        if self.tech == "contig": continue
                        sub_read_quals[-1] += self.this_read_quals[read_pos:read_pos + cigartuple[1]]

                elif cigartuple[0] == 5 or cigartuple[0] == 4:
                    continue

                read_pos += cigartuple[1]

            elif cigartuple[0] in [2, 3]:
                sub_read_str.extend([""] * cigartuple[1])

                cigars_list.extend(['D'] *cigartuple[1])

                if self.tech == "contig": continue
                sub_read_quals.extend([""] * cigartuple[1])
            else:
                return -1
        self.this_read_list = sub_read_str
        self.this_quals_list = sub_read_quals
        self.this_read_quals = ""

        self.cigars_list=cigars_list

    def get_read_str_speedUp(self):
        self.this_ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, start=self.align_start,
                                                                  end=self.align_end).upper()
        self.this_ref_list = list(self.this_ref_str)
        sub_read_str = []
        sub_read_quals = []

        read_pos = 0

        cigartuples_array=np.array(self.cigartuples)
        for cigartuple in cigartuples_array:

            if cigartuple[0] in [0, 7, 8]:

                sub_read_str.extend(list(self.this_read_str[read_pos:read_pos + cigartuple[1]]))

                sub_read_quals.extend(list(self.this_read_quals[read_pos:read_pos + cigartuple[1]]))
                read_pos += cigartuple[1]

            elif cigartuple[0] in [1, 4, 5]:

                if cigartuple[0] == 1:
                    if len(sub_read_str) == 0:
                        continue
                    else:
                        sub_read_str[-1] += self.this_read_str[read_pos:read_pos + cigartuple[1]]

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

    def get_repeat_length(self, ms_start, ms_end):
        query_repeat_length = len(
            "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
        return query_repeat_length

    def get_repeat_length_all_ms(self):
        self.microsatellites_num = len(self.microsatellites)
        repeat_lengths = {}
        for ms_id, ms in self.microsatellites.items():
            ms_start = ms.start
            ms_end = ms.end
            query_repeat_length = len(
                "".join(self.this_read_list[ms_start - 1 - self.align_start:ms_end - self.align_start - 1]))
            repeat_lengths[ms_id] = query_repeat_length
        self.repeat_lengths = repeat_lengths

    def get_quals(self, q_start, q_end):
        quals_list = self.this_quals_list[q_start - self.align_start - 1:q_end - self.align_start - 1]
        return quals_list, self.strand

    def get_ms_info_one_read(self,compar_software='minimap2'):
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
                flankLen=1
                str_varia = "".join(self.this_read_list[
                            ms_start - flankLen - self.align_start:ms_end - self.align_start +flankLen])

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

    def predict_hmm(self,str_varia,model,flankLen,ms_motif):

        obs_list_temp = []
        for i in str_varia[:flankLen]:
            if i == 'A':
                obs_list_temp.append(0)
            elif i == 'T':
                obs_list_temp.append(1)
            elif i == 'C':
                obs_list_temp.append(2)
            elif i == 'G':
                obs_list_temp.append(3)
        for i in range(0, len(str_varia[flankLen:-flankLen]), len(ms_motif)):
            if (i + len(ms_motif) - 1) >= len(str_varia[flankLen:-flankLen]): break
            str = str_varia[flankLen:-flankLen][i]

            if str == 'A':
                obs_list_temp.append(0)
            elif str == 'T':
                obs_list_temp.append(1)
            elif str == 'C':
                obs_list_temp.append(2)
            elif str == 'G':
                obs_list_temp.append(3)
        for i in str_varia[-flankLen:]:
            if i == 'A':
                obs_list_temp.append(0)
            elif i == 'T':
                obs_list_temp.append(1)
            elif i == 'C':
                obs_list_temp.append(2)
            elif i == 'G':
                obs_list_temp.append(3)
        obs_list = np.array([obs_list_temp]).T
        predict_list = model.predict(obs_list, lengths=None)
        predict_count = Counter(predict_list)
        dict_correct = {}
        dict_correct['match'] = predict_count.most_common()[0][0]
        if len(predict_count.most_common()) >= 2:
            dict_correct['ins'] = predict_count.most_common()[1][0]
        else:
            dict_correct['ins'] = 'no ins'
        if len(predict_count.most_common()) >= 3:
            dict_correct['del'] = predict_count.most_common()[2][0]
        else:
            dict_correct['del'] = 'no del'

        count = 0
        temp_list = predict_list[flankLen:]
        temp_list = temp_list[:-flankLen]
        for i in temp_list:
            if i == dict_correct['match']:
                count += 1 * len(ms_motif)
            elif i == dict_correct['del']:
                count -= 1 * len(ms_motif)
            elif i == dict_correct['ins']:
                count += 1 * len(ms_motif)
        for i in predict_list[flankLen:]:
            if i == dict_correct['del']:
                count -= 1
            elif i == dict_correct['ins']:
                count += 1
        for i in predict_list[:-flankLen]:
            if i == dict_correct['del']:
                count -= 1
            elif i == dict_correct['ins']:
                count += 1

        return count

    def predict_hmm_2(self, str_varia, model, flankLen, ms_motif):

        obs_list_temp = []
        for i in str_varia:
            if i == 'A':
                obs_list_temp.append(0)
            elif i == 'T':
                obs_list_temp.append(1)
            elif i == 'C':
                obs_list_temp.append(2)
            elif i == 'G':
                obs_list_temp.append(3)

        obs_list = np.array([obs_list_temp]).T

        predict_list = model.predict(obs_list, lengths=None)
        predict_count = Counter(predict_list)
        dict_correct = {}
        dict_correct['match'] = predict_count.most_common()[0][0]
        if len(predict_count.most_common()) >= 2:
            dict_correct['ins'] = predict_count.most_common()[1][0]
        else:
            dict_correct['ins'] = 'no ins'
        if len(predict_count.most_common()) >= 3:
            dict_correct['del'] = predict_count.most_common()[2][0]
        else:
            dict_correct['del'] = 'no del'

        count = 0
        for i in predict_list:
            if i == dict_correct['match']:
                count += 1 * len(ms_motif)
            elif i == dict_correct['del']:
                count -= 1 * len(ms_motif)
            elif i == dict_correct['ins']:
                count += 1 * len(ms_motif)

        return count

    def predict_hmm_3(self, str_varia, model, ms_motif):

        obs_list_temp = []
        for i in str_varia:
            if i == 'A':
                obs_list_temp.append(1)
            elif i == 'T':
                obs_list_temp.append(2)
            elif i == 'C':
                obs_list_temp.append(3)
            elif i == 'G':
                obs_list_temp.append(4)

        obs_list = np.array([obs_list_temp]).T
        predict_list = model.predict(obs_list, lengths=None)
        count=0
        for nums in predict_list:
            if nums==0 or nums==1 or nums==2 or nums==3:
                count-=1
            elif nums==4 or nums==5 or nums==6 or nums==7:
                count+=1
            else:
                count

        return count

    def predict_hmm_4(self, str_varia, model, ms_motif):

        obs_list_temp = []
        for i in str_varia:
            if i == 'A':
                obs_list_temp.append(1)
            elif i == 'T':
                obs_list_temp.append(2)
            elif i == 'C':
                obs_list_temp.append(3)
            elif i == 'G':
                obs_list_temp.append(4)
        obs_list = np.array([obs_list_temp]).T

        predict_list = model.predict(obs_list)

        predict_count = Counter(predict_list)
        dict_correct = {}

        dict_correct['match'] = predict_count.most_common()[0][0]
        if len(predict_count.most_common()) >= 2:
            dict_correct['ins'] = predict_count.most_common()[1][0]
        else:
            dict_correct['ins'] = 'no ins'
        if len(predict_count.most_common()) >= 3:
            dict_correct['del'] = predict_count.most_common()[2][0]
        else:
            dict_correct['del'] = 'no del'

        count = 0
        for i in predict_list:
            if i == dict_correct['match']:
                count += 1
            elif i == dict_correct['del']:
                count -= 1
            elif i == dict_correct['ins']:
                count += 2

        return count

    def predict_hmm_5(self, str_varia, model,model_gmm, ms_motif):

        obs_list_temp = []
        for i in str_varia:
            if i == 'A':
                obs_list_temp.append(1)
            elif i == 'T':
                obs_list_temp.append(2)
            elif i == 'C':
                obs_list_temp.append(3)
            elif i == 'G':
                obs_list_temp.append(4)
        obs_list = np.array([obs_list_temp]).T
        predict_list = model.predict(obs_list, lengths=None)
        predict_count = Counter(predict_list)
        dict_correct = {}

        dict_correct['match'] = predict_count.most_common()[0][0]
        if len(predict_count.most_common()) >= 2:
            dict_correct['ins'] = predict_count.most_common()[1][0]
        else:
            dict_correct['ins'] = 'no ins'
        if len(predict_count.most_common()) >= 3:
            dict_correct['del'] = predict_count.most_common()[2][0]
        else:
            dict_correct['del'] = 'no del'

        predict_list_corrected=[]
        for i in predict_list:
            if i == dict_correct['match']:
                predict_list_corrected.append(1)
            elif i == dict_correct['del']:
                predict_list_corrected.append(3)
            elif i == dict_correct['ins']:
                predict_list_corrected.append(2)

        count=0
        temp_MI=''
        for i in range(len(predict_list_corrected)):
            if predict_list_corrected[i]==1:
                count+=1
                temp_MI='M'
            elif predict_list_corrected[i]==2:
                count+=2
                temp_MI = 'I'
            elif predict_list_corrected[i]==3:
                if temp_MI=='':
                    temp_MI='M'
                temp_key=str_varia[i]+temp_MI
                count-=math.ceil(model_gmm[temp_key])*1

        count=round(count,0)

        return count

    def predict_hmm_6(self, str_varia, model, model_gmm, ms_motif,ms_id):

        obs_list_temp = []
        for i in str_varia:
            if i == 'A':
                obs_list_temp.append(1)
            elif i == 'T':
                obs_list_temp.append(2)
            elif i == 'C':
                obs_list_temp.append(3)
            elif i == 'G':
                obs_list_temp.append(4)
        obs_list = np.array([obs_list_temp]).T
        predict_list = model.predict(obs_list, lengths=None)

        dict_ATCG={'A':[],'T':[],'C':[],'G':[]}
        for i in range(len(predict_list)):
            if predict_list[i] not in dict_ATCG[str_varia[i]]:
                dict_ATCG[str_varia[i]].append(predict_list[i])

        dict_preList_count={}
        for i in predict_list:
            if i not in dict_preList_count.keys():
                dict_preList_count[i]=1
            else:
                dict_preList_count[i]+=1

        dict_ATCG_count={'A':{},'T':{},'C':{},'G':{}}
        for key,value in dict_ATCG.items():
            for i in value:
                dict_ATCG_count[key][i]=dict_preList_count[i]

        dict_ATCG_sorted={'A':[],'T':[],'C':[],'G':[]}
        for keys,values in dict_ATCG_count.items():
            temp_dict_sorted=sorted(dict_ATCG_count[keys].items(), key=lambda x: x[1], reverse=True)
            for i in range(len(temp_dict_sorted)):

                dict_ATCG_sorted[keys].append(int(temp_dict_sorted[i][0]))

        del_num_list,m_ins_num_list=[],[]
        for key,value in dict_ATCG_sorted.items():
            if len(value)==0:continue
            elif len(value)==1:
                m_ins_num_list.append(value[0])
            elif len(value)==2:
                m_ins_num_list.append(value[0])
                m_ins_num_list.append(value[1])
            elif len(value)==3:
                m_ins_num_list.append(value[0])
                m_ins_num_list.append(value[1])
                del_num_list.append(value[2])
            else:
                for i in range(len(value)):
                    if i==len(value)-1:
                        del_num_list.append(value[i])
                    else:
                        m_ins_num_list.append(value[i])
        count=0
        for i in predict_list:
            if i in del_num_list:
                count+=0
                # count+=1
            elif i in m_ins_num_list:
                count+=1

        return count

    def predict_hmm_7(self, str_varia, model, model_gmm, ms_motif):

        obs_list_temp = []
        for i in str_varia:
            if i == 'A':
                obs_list_temp.append(1)
            elif i == 'T':
                obs_list_temp.append(2)
            elif i == 'C':
                obs_list_temp.append(3)
            elif i == 'G':
                obs_list_temp.append(4)
        obs_list = np.array([obs_list_temp]).T
        predict_list = model.predict(obs_list, lengths=None)

        dict_ATCG = {'A': [], 'T': [], 'C': [], 'G': []}
        for i in range(len(predict_list)):
            if predict_list[i] not in dict_ATCG[str_varia[i]]:
                dict_ATCG[str_varia[i]].append(predict_list[i])

        dict_preList_count = {}
        for i in predict_list:
            if i not in dict_preList_count.keys():
                dict_preList_count[i] = 1
            else:
                dict_preList_count[i] += 1

        dict_ATCG_count = {'A': {}, 'T': {}, 'C': {}, 'G': {}}
        for key, value in dict_ATCG.items():
            for i in value:
                dict_ATCG_count[key][i] = dict_preList_count[i]

        dict_ATCG_sorted = {'A': [], 'T': [], 'C': [], 'G': []}
        for keys, values in dict_ATCG_count.items():
            temp_dict_sorted = sorted(dict_ATCG_count[keys].items(), key=lambda x: x[1], reverse=True)
            for i in range(len(temp_dict_sorted)):
                dict_ATCG_sorted[keys].append(int(temp_dict_sorted[i][0]))

        del_num_list, m_num_list, ins_num_list = [], [], []
        for key, value in dict_ATCG_sorted.items():
            if len(value) == 0:
                continue
            elif len(value) == 1:
                m_num_list.append(value[0])
            elif len(value) == 2:
                m_num_list.append(value[0])
                ins_num_list.append(value[1])
            else:
                m_num_list.append(value[0])
                ins_num_list.append(value[1])
                del_num_list.append(value[2])

        count = 0
        temp_MI=''
        for i in range(len(predict_list)):
            if predict_list[i] in del_num_list:
                if temp_MI=='':
                    temp_MI='M'
                temp_key = str_varia[i] + temp_MI
                count -= math.ceil(model_gmm[temp_key])*1
            elif predict_list[i] in m_num_list:
                count += 1
                temp_MI = 'M'
            elif predict_list[i] in ins_num_list:
                count += 2
                temp_MI = 'I'

        return count

    def predict_hmm_8(self, str_varia, model,model_gmm, ms_motif):

        if 'N' in str_varia:
            return len(str_varia)

        obs_list_temp = []
        for i in str_varia:
            if i == 'A':
                obs_list_temp.append(1)
            elif i == 'T':
                obs_list_temp.append(2)
            elif i == 'C':
                obs_list_temp.append(3)
            elif i == 'G':
                obs_list_temp.append(4)
        obs_list = np.array([obs_list_temp]).T
        predict_list = model.predict(obs_list, lengths=None)
        predict_count = Counter(predict_list)
        dict_correct = {}

        dict_correct['match'] = predict_count.most_common()[0][0]
        if len(predict_count.most_common()) >= 2:
            dict_correct['ins'] = predict_count.most_common()[1][0]
        else:
            dict_correct['ins'] = 'no ins'
        if len(predict_count.most_common()) >= 3:
            dict_correct['del'] = predict_count.most_common()[2][0]
        else:
            dict_correct['del'] = 'no del'

        predict_list_corrected=[]
        for i in predict_list:
            if i == dict_correct['match']:
                predict_list_corrected.append(1)  # 匹配
            elif i == dict_correct['del']:
                predict_list_corrected.append(3)
            elif i == dict_correct['ins']:
                predict_list_corrected.append(2)

        count=0
        temp_MI=''
        for i in range(len(predict_list_corrected)):
            if predict_list_corrected[i]==1:
                count+=1
                temp_MI='M'
            elif predict_list_corrected[i]==2:
                count+=1
                temp_MI = 'I'
            elif predict_list_corrected[i]==3:
                if temp_MI=='':
                    temp_MI='M'
                temp_key=str_varia[i]+temp_MI
                if model_gmm[temp_key]!=[]:
                    count+=model_gmm[temp_key]
                else:
                    count+=0

        return count

    def predict_hmm_9(self, str_varia, model, ms_motif):

        obs_list_temp = []
        for i in str_varia:
            if i == 'A':
                obs_list_temp.append(1)
            elif i == 'T':
                obs_list_temp.append(2)
            elif i == 'C':
                obs_list_temp.append(3)
            elif i == 'G':
                obs_list_temp.append(4)
        obs_list = np.array([obs_list_temp]).T
        predict_list = model.predict(obs_list, lengths=None)
        predict_count = Counter(predict_list)
        dict_correct = {}

        dict_correct['match'] = predict_count.most_common()[0][0]
        if len(predict_count.most_common()) >= 2:
            dict_correct['ins'] = predict_count.most_common()[1][0]
        else:
            dict_correct['ins'] = 'no ins'
        if len(predict_count.most_common()) >= 3:
            dict_correct['del'] = predict_count.most_common()[2][0]
        else:
            dict_correct['del'] = 'no del'

        count=0
        for i in predict_list:
            if i == dict_correct['match']:
                count += 1
            elif i == dict_correct['del']:
                count += 0
            elif i == dict_correct['ins']:
                count += 1

        # count=0
        # for i in predict_list:
        #     if i == dict_correct['match']:
        #         count += 1
        #     elif i == dict_correct['del']:
        #         count += 1
        #     # elif i==2:count+=2
        #     elif i == dict_correct['ins']:
        #         # count += 1 * len(ms_motif)
        #         count += 1

        return count

    def predict_hmm_9_draw(self, str_varia, model, ms_motif):

        obs_list_temp = []
        for i in str_varia:
            if i == 'A':
                obs_list_temp.append(1)
            elif i == 'T':
                obs_list_temp.append(2)
            elif i == 'C':
                obs_list_temp.append(3)
            elif i == 'G':
                obs_list_temp.append(4)
        obs_list = np.array([obs_list_temp]).T
        predict_list = model.predict(obs_list, lengths=None)
        predict_count = Counter(predict_list)
        dict_correct = {}

        dict_correct['match'] = predict_count.most_common()[0][0]
        if len(predict_count.most_common()) >= 2:
            dict_correct['ins'] = predict_count.most_common()[1][0]
        else:
            dict_correct['ins'] = 'no ins'
        if len(predict_count.most_common()) >= 3:
            dict_correct['del'] = predict_count.most_common()[2][0]
        else:
            dict_correct['del'] = 'no del'

        count = 0
        predictRes=[]
        for i in predict_list:
            if i == dict_correct['match']:
                count += 1
                predictRes.append('M')
            elif i == dict_correct['del']:
                count += 0
                predictRes.append('D')
            elif i == dict_correct['ins']:
                count += 1
                predictRes.append('I')

        # count=0
        # for i in predict_list:
        #     if i == dict_correct['match']:
        #         count += 1
        #     elif i == dict_correct['del']:
        #         count += 1
        #     # elif i==2:count+=2
        #     elif i == dict_correct['ins']:
        #         # count += 1 * len(ms_motif)
        #         count += 1

        return count,predictRes

    def get_ms_info_one_read_improved(self,model,compar_software='minimap2'):

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
                flankLen=5
                str_varia = "".join(self.this_read_list[
                            ms_start - flankLen - self.align_start:ms_end - self.align_start +flankLen])
                if len(str_varia)==0:continue

                query_repeat_length=self.predict_hmm_2(str_varia,model,flankLen,ms_motif)

            elif compar_software == 'lra':
                query_repeat_length = len(
                    "".join(self.this_read_list[ms_start + 3 - self.align_start:ms_end - self.align_start+0]))
            else:
                print('compare_software is wrong!')

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

    def get_ms_info_one_read_improved_2(self,model, model_gmm,compar_software='minimap2'):

        self.microsatellites_num = len(self.microsatellites)
        read_muts = {}
        repeat_lengths = {}

        repeat_lengths_GmmHmm = {}

        hmmRes_draw = {}

        for ms_id, ms in self.microsatellites.items():
            ms_start = ms.start
            ms_end = int(ms.end)
            ms_start_pre = ms.start_pre
            ms_end_suf = int(ms.end_suf)

            ms_motif = ms.repeat_unit

            flankLen=1
            str_varia = "".join(self.this_read_list[
                        ms_start - flankLen - self.align_start:ms_end - self.align_start +flankLen])

            if len(str_varia) >= 2:
                if str_varia[0] != ms_motif[0]:
                    str_varia = str_varia[1:]
                if str_varia[-1] != ms_motif[-1]:
                    str_varia = str_varia[:-1]
            if len(str_varia) == 0: continue

            if len(set(str_varia)) == 1 and 'N' in str_varia:continue

            # query_repeat_length=self.predict_hmm(str_varia,model,flankLen,ms_motif)
            # query_repeat_length=self.predict_hmm_2(str_varia,model,flankLen,ms_motif)
            # query_repeat_length=self.predict_hmm_3(str_varia,model,ms_motif)
            # query_repeat_length=self.predict_hmm_4(str_varia,model,ms_motif)
            # query_repeat_length=int(self.predict_hmm_5(str_varia,model, model_gmm,ms_motif))
            # query_repeat_length=int(self.predict_hmm_6(str_varia,model, model_gmm,ms_motif))
            # query_repeat_length=self.predict_hmm_6(str_varia,model, model_gmm,ms_motif,ms_id)
            # query_repeat_length = int(self.predict_hmm_7(str_varia, model, model_gmm, ms_motif))
            # if query_repeat_length!=len(str_varia):
            #     print(query_repeat_length,'**',len(str_varia))

            # query_repeat_length=len(str_varia)
            # query_repeat_length_GmmHMM=self.predict_hmm_4(str_varia,model,ms_motif)

            query_repeat_length_GmmHMM=11

            # query_repeat_length=self.predict_hmm_5(str_varia,model,model_gmm,ms_motif)
            # query_repeat_length_GmmHMM=int(self.predict_hmm_5(str_varia,model, model_gmm,ms_motif))

            # query_repeat_length = len(str_varia)
            # query_repeat_length_GmmHMM = int(self.predict_hmm_5(str_varia, model, model_gmm, ms_motif))

            # query_repeat_length = int(self.predict_hmm_7(str_varia, model, model_gmm, ms_motif))
            # query_repeat_length_GmmHMM = int(self.predict_hmm_5(str_varia, model, model_gmm, ms_motif))

            # query_repeat_length=len(str_varia)
            # query_repeat_length_GmmHMM = int(self.predict_hmm_7(str_varia, model, model_gmm, ms_motif))

            # query_repeat_length=len(str_varia)
            # query_repeat_length_GmmHMM=int(self.predict_hmm_6(str_varia,model, model_gmm,ms_motif))

            # query_repeat_length=len(str_varia)
            # query_repeat_length_GmmHMM=int(self.predict_hmm_5(str_varia,model, model_gmm,ms_motif))
            # query_repeat_length_GmmHMM=round(self.predict_hmm_5(str_varia,model, model_gmm,ms_motif),0)
            # query_repeat_length_GmmHMM=math.ceil(self.predict_hmm_5(str_varia,model, model_gmm,ms_motif))

            # query_repeat_length = round(self.predict_hmm_8(str_varia,model, model_gmm,ms_motif),0)
            # query_repeat_length = self.predict_hmm_8(str_varia,model, model_gmm,ms_motif)
            # query_repeat_length_GmmHMM = int(self.predict_hmm_6(str_varia, model, model_gmm, ms_motif))
            # if query_repeat_length>(int(ms_end)-int(ms_start)):
            #     query_repeat_length=math.ceil(query_repeat_length)
            # elif query_repeat_length<(int(ms_end)-int(ms_start)):
            #     query_repeat_length=int(query_repeat_length)

            predictRes=[]

            if model==None:
                query_repeat_length=len(str_varia)
            else:
                query_repeat_length = self.predict_hmm_9(str_varia,model,ms_motif)

            # query_repeat_length=len(str_varia)

            if query_repeat_length == -len(str_varia):query_repeat_length=-query_repeat_length
            if query_repeat_length_GmmHMM==-len(str_varia):query_repeat_length_GmmHMM=-query_repeat_length_GmmHMM

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

            repeat_lengths_GmmHmm[ms_id] = query_repeat_length_GmmHMM

            predictRes=self.cigars_list[ms_start - flankLen - self.align_start:ms_end - self.align_start + flankLen]

            if ms_id not in hmmRes_draw.keys():
                hmmRes_draw[ms_id]={}

            if 'ms_start' not in list(hmmRes_draw[ms_id].keys()):
                hmmRes_draw[ms_id]={'ms_start':ms_start,'ms_end':ms_end,'motif':ms_motif,'predict_res':[predictRes]}
            else:
                hmmRes_draw[ms_id]['predict_res'].append(predictRes)

        self.repeat_lengths = repeat_lengths
        self.mut_info = read_muts
        self.repeat_lengths_GmmHmm=repeat_lengths_GmmHmm
        self.hmmRes_draw=hmmRes_draw

    def get_ms_info_one_read_improved_3(self, compar_software='minimap2'):

        self.microsatellites_num = len(self.microsatellites)
        read_muts = {}
        repeat_lengths = {}

        for ms_id, ms in self.microsatellites.items():
            ms_start = ms.start
            ms_end = ms.end
            ms_start_pre = ms.start_pre
            ms_end_suf = ms.end_suf

            ms_motif = ms.repeat_unit
            flankLen = 1
            str_varia = "".join(self.this_read_list[
                                ms_start - flankLen - self.align_start:ms_end - self.align_start + flankLen])
            if len(str_varia) >= 2:
                if str_varia[0] != ms_motif[0]:
                    str_varia = str_varia[1:]
                if str_varia[-1] != ms_motif[-1]:
                    str_varia = str_varia[:-1]

            query_repeat_length = len(str_varia)

            if query_repeat_length == -len(str_varia): query_repeat_length = -query_repeat_length

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
                    print('Read.py:pot', pot, ' len', len(self.this_read_list))
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

    def get_all_str_varia(self,flankLen_left=1,flankLen_right=0):

        out_list,out_cigar_list,out_varia_num_list=[],[],[]
        for ms_id, ms in self.microsatellites.items():
            ms_start = ms.start
            ms_end = int(ms.end)
            str_varia = "".join(self.this_read_list[
                        ms_start - flankLen_left - self.align_start:ms_end - self.align_start +flankLen_right])
            str_cigar="".join(self.cigars_list[
                        ms_start - flankLen_left - self.align_start:ms_end - self.align_start +flankLen_right])
            out_list.append(str_varia)
            out_cigar_list.append(str_cigar)

            temp=[]
            for strr in str_varia:
                if strr=='A':
                    temp.append(1)
                elif strr=='T':
                    temp.append(2)
                elif strr=='C':
                    temp.append(3)
                elif strr=='G':
                    temp.append(4)
            out_varia_num_list.append(temp)
            return out_list,out_cigar_list,out_varia_num_list


