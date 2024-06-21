#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from global_dict import *
import numpy as np
import pysam
from units import *
from Read import Read
from gmm import get_repeat_gmm
from collections import Counter
from sklearn.cluster import KMeans
from sklearn import mixture
import math

import logging
mpl_logger = logging.getLogger('matplotlib')
mpl_logger.setLevel(logging.ERROR)

import matplotlib.patches as patches
import matplotlib.pyplot as plt
import os


class Microsatellite:

    def __init__(self, ms_info, only_simple=True):
        self.chrom = ms_info["chr"]
        self.start = ms_info["pos"]
        self.prefix = ms_info["prefix"]
        self.suffix = ms_info["suffix"]
        self.ms_id = self.chrom + "_" + str(self.start)
        self.repeat_times = ms_info["repeatTimes"]
        self.repeat_unit = ms_info["motif"]
        self.repeat_unit_len = ms_info["motifLen"]
        self.reference = ms_info["reference"]
        self.repeat_len = self.repeat_times * self.repeat_unit_len
        self.end = self.start + self.repeat_len
        self.start_pre = self.start - ms_info["prefix_len"]
        self.end_suf = self.end + ms_info["suffix_len"]
        self.suffix_len = ms_info["suffix_len"]
        self.prefix_len = ms_info["prefix_len"]
        self.minimum_support_reads = 0
        self.min_allele_fraction = 0
        self.sequencing_error = get_value("paras")["sequencing_error"]
        self.reads_info = {}
        self.muts = {}
        self.length_dis_reads = {}
        self.depth = 0
        self.check = True
        self.check_status = []
        self.ms_error = []
        self.ref_all = pysam.FastaFile(self.reference).fetch(self.chrom, self.start_pre - 1, self.end_suf + 1)
        self.ref_str = ""
        self.ref_str_ms = ""
        self.ref_str_mut = ""
        self.alt_str = ""
        self.alt_ms = (".",)
        self.alt = (".",)
        self.dis_stat = True
        self.mut_start = self.start
        self.mut_end = self.end
        self.ms_dis_info = {}
        self.ms_dis = {}
        self.ms_dis_forward = {}
        self.ms_dis_reversed = {}
        self.ms_dis_hap1 = {}
        self.ms_dis_hap2 = {}
        self.ms_dis_hap0 = {}
        self.support_reads = 0
        self.support_hap0 = 0
        self.support_hap1 = 0
        self.support_hap2 = 0
        self.support_forward = 0
        self.support_reversed = 0
        self.query_repeat_length = 0
        self.deletions = {}
        self.ms_mutation = False
        self.insertions = {}
        self.mismatches = {}
        self.reads_phased = True
        self.model_stat = True
        self.format_GT_ms = (0, 0)
        self.format_AL_ms = "/".join(["0", "0"])
        self.format_DP_ms = "/".join(["0", "0", "0"])
        self.format_QL_ms = "/".join(["0", "0", "0"])
        self.hap1_ms_mut = False
        self.hap2_ms_mut = False
        self.hap1_ms_repeat = 0
        self.hap2_ms_repeat = 0

        self.qual_ms = 0
        self.qual_ms_hap1 = 0
        self.qual_ms_hap2 = 0
        self.phased = False
        self.report_micro = True
        self.only_simple = only_simple
        self.minimum_support_reads = get_value("paras")["minimum_support_reads"]
        self.min_allele_fraction = get_value("paras")["min_allele_fraction"]
        self.maximum_distance_of_two_complex_events = get_value("paras")["maximum_distance_of_two_complex_events"]

        self.kTemp = 0

        self.support_reads_totalNum = 0
        self.support_reads_mutNum = [0,0]

        self.GmmHmm_res={}

        self.hmmRes_draw={}

        if only_simple:
            self.report_indel = False
            self.report_snv = False
            self.report_complex = False
        else:
            self.report_indel = True
            self.report_snv = True
            self.report_complex = True
            self.format_GT = (0, 0)
            self.format_AL = "/".join(["0", "0"])
            self.format_DP = "/".join(["0", "0", "0"])
            self.format_QL = "/".join(["0", "0", "0"])

    def set_reads_info(self, reads_info):
        self.reads_info = reads_info
        self.depth = len(self.reads_info)
        self.minimum_support_reads = max(round(self.depth * self.min_allele_fraction), self.minimum_support_reads)

    def set_muts_info(self, muts):
        self.muts = muts
        self.depth = len(self.muts)
        self.minimum_support_reads = max(round(self.depth * self.min_allele_fraction), self.minimum_support_reads)

    def set_len_info(self,ms_id,gmmhmm_info):
        gmm=mixture.GaussianMixture(n_components=1,random_state=0)
        gmm.fit(np.array(gmmhmm_info).reshape(-1,1))
        self.GmmHmm_res[ms_id]=gmm.means_[0][0]

    def set_read_dis_info(self, reads_info):
        if len(reads_info) <= 0:
            return
        dis = {}
        dis_strand = {True: {}, False: {}}
        dis_hap = {0: {}, 1: {}, 2: {}}

        for read_id, ms_info in reads_info.items():
            repeat_length, strand, hap = ms_info

            if repeat_length not in dis:
                dis[repeat_length] = 0
            dis[repeat_length] += 1
            if repeat_length not in dis_hap[hap]:
                dis_hap[hap][repeat_length] = 0
            dis_hap[hap][repeat_length] += 1
            if repeat_length not in dis_strand[strand]:
                dis_strand[strand][repeat_length] = 0
            dis_strand[strand][repeat_length] += 1

        new_dis = {}
        depth = sum(dis.values())
        errors = []
        for rp, times in dis.items():
            if times > self.sequencing_error * depth:
                new_dis[rp] = times
            else:
                errors.append(rp)
        for rp in errors:
            if rp in dis_hap[0]:
                dis_hap[0].pop(rp)
            if rp in dis_hap[1]:
                dis_hap[1].pop(rp)
            if rp in dis_hap[2]:
                dis_hap[2].pop(rp)
            if rp in dis_strand[True]:
                dis_strand[True].pop(rp)
            if rp in dis_strand[False]:
                dis_strand[False].pop(rp)

        self.ms_error = errors
        self.dis_stat = True if self.depth > 0 else False
        self.ms_dis = dis

        for values in dis.values():
            self.support_reads_totalNum += values

        self.ms_dis_hap0 = dis_hap[0]
        self.ms_dis_hap1 = dis_hap[1]
        self.ms_dis_hap2 = dis_hap[2]
        self.ms_dis_forward = dis_strand[True]
        self.ms_dis_reversed = dis_strand[False]
        self.query_repeat_length = get_max_support_index(dis)
        dis_hap0_num = sum(dis_hap[0].values())
        dis_hap1_num = sum(dis_hap[1].values())
        dis_hap2_num = sum(dis_hap[2].values())
        self.support_hap0 = dis_hap0_num
        self.support_hap1 = dis_hap1_num
        self.support_hap2 = dis_hap2_num

        self.depth = dis_hap0_num + dis_hap1_num + dis_hap2_num
        self.support_reads = dis_hap0_num + dis_hap1_num + dis_hap2_num

        if abs(dis_hap1_num - dis_hap2_num) > self.depth * 0.4:
            self.reads_phased = False
        elif dis_hap0_num > self.depth * 0.5:
            self.reads_phased = False
        elif dis_hap1_num < 2 or dis_hap2_num < 2:
            self.reads_phased = False
        else:
            self.reads_phased = True

    def get_dis(self):
        samfile = pysam.Samfile(get_value("paras")["input"])
        reads = {}
        for alignment in samfile.fetch(self.chrom, self.start_pre - 1, self.end_suf + 1):
            if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary or alignment.is_supplementary:
                continue
            if alignment.reference_start > self.start_pre - 1 or alignment.reference_end < self.end_suf + 1:
                continue
            if len(alignment.query_sequence) < 2:
                continue
            read_ms = Read(read_id="", alignment=alignment, reference=self.reference, chrom=self.chrom)
            read_ms.get_read_str()
            q_repeat_length = read_ms.get_repeat_length(self.start, self.end)
            if q_repeat_length not in reads:
                reads[q_repeat_length] = 1
            else:
                reads[q_repeat_length] += 1
        return reads

    def get_dis_qual(self):
        samfile = pysam.Samfile(get_value("paras")["input"])
        reads = {}
        quals = {}
        prefix_forward = []
        suffix_forward = []
        ms_forward = []
        prefix_reversed = []
        suffix_reversed = []
        ms_reversed = []
        num_forward = 0
        num_reversed = 0
        for alignment in samfile.fetch(self.chrom, self.start_pre - 1, self.end_suf + 1):
            if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary or alignment.is_supplementary:
                continue
            if alignment.reference_start > self.start_pre - 1 or alignment.reference_end < self.end_suf + 1:
                continue
            if len(alignment.query_sequence) < 2:
                continue
            read_ms = Read(read_id="", alignment=alignment, reference=self.reference, chrom=self.chrom)
            read_ms.get_read_str()
            q_repeat_length = read_ms.get_repeat_length(self.start, self.end)
            if q_repeat_length not in reads:
                reads[q_repeat_length] = 1
            else:
                reads[q_repeat_length] += 1
            prefix_qual = read_ms.get_quals(self.start_pre, self.start)
            suffix_qual = read_ms.get_quals(self.end, self.end_suf)
            ms_qual = read_ms.get_quals(self.start, self.end)
            if prefix_qual[1]:
                num_forward += 1
                prefix_forward.append(np.array(list(map(str2int, prefix_qual[0]))))
                suffix_forward.append(np.array(list(map(str2int, suffix_qual[0]))))
                ms_forward.append(np.array(list(map(str2int, ms_qual[0]))))
            else:
                num_reversed += 1
                prefix_reversed.append(np.array(list(map(str2int, prefix_qual[0]))))
                suffix_reversed.append(np.array(list(map(str2int, suffix_qual[0]))))
                ms_reversed.append(np.array(list(map(str2int, ms_qual[0]))))
        quals["num_forward"] = num_forward
        quals["prefix_forward"] = np.array(prefix_forward)
        quals["suffix_forward"] = np.array(suffix_forward)
        quals["ms_forward"] = np.array(ms_forward)

        quals["num_reversed"] = num_reversed
        quals["prefix_reversed"] = np.array(prefix_reversed)
        quals["suffix_reversed"] = np.array(suffix_reversed)
        quals["ms_reversed"] = np.array(ms_reversed)
        return reads, quals

    def build_patterns(self):
        if self.reads_phased:
            hap1 = []
            hap2 = []
            for read_id, read_info in self.muts.items():
                if read_info.hap == 1 and read_info.other_mutation:
                    hap1.append(read_info)
                elif read_info.hap == 2 and read_info.other_mutation:
                    hap2.append(read_info)
            if len(hap1) >= self.minimum_support_reads:
                mutation_events_hap1 = self.get_mut_one_hap(hap1, self.hap1_ms_mut)
            if len(hap2) >= self.minimum_support_reads:
                mutation_events_hap2 = self.get_mut_one_hap(hap2, self.hap2_ms_mut)

        else:
            pass

    def deletion_merge(self):
        for read_id, read_info in self.muts.items():
            deletions_len = len(read_info.deletions)
            deletions = []

            if deletions_len == 0:
                pass
            elif deletions_len == 1:
                deletions = [[read_info.deletions[0][0], 1, read_info.deletions[0][-1]]]
            else:
                band = []
                current_id = read_info.deletions[0][0]
                deletion_index = {current_id: 1}
                for pos_id in range(1, deletions_len):
                    band.extend(read_info.deletions[pos_id][2])
                    if read_info.deletions[pos_id][0] == read_info.deletions[pos_id - 1][0] + 1:
                        deletion_index[current_id] += 1
                    else:
                        current_id = read_info.deletions[pos_id][0]
                        deletion_index[current_id] = 1
                for pos in sorted(deletion_index.keys()):
                    deletions.append([pos, deletion_index[pos], list(set(band))])
            self.muts[read_id].deletions = deletions

    def get_alt(self, alt_dict, offset):

        alt_list = list(self.ref_str)
        for pos, info in alt_dict:
            alt_list[pos - offset] = info
        return "".join(alt_list)

    def call_init(self):
        if self.depth == 0:
            self.check = False
            self.check_status.append("No_read_covered")
            self.dis_stat = False
            self.ref_str_ms = pysam.FastaFile(self.reference).fetch(self.chrom, self.start, self.end + 1)
            self.alt_ms = "."
            self.report_micro = False
            self.report_indel = False
            self.report_snv = False
            self.report_complex = False
            return False
        else:
            self.report_micro = True
            if not self.only_simple:
                self.report_indel = False
                self.report_snv = False
                self.report_complex = False
            else:
                self.report_indel = True
                self.report_snv = True
                self.report_complex = True

            return True

    def pre_cluster(self,mode='GaussianMixture',n_clusters = 2):

        ms_dis_list = [[item[0]] * int(item[1]) for item in self.ms_dis.items()]
        ms_dis_list = sum(ms_dis_list, [])
        X = np.array(ms_dis_list).reshape(-1, 1)

        if len(X) == 1:
            X = X.tolist()
            X.append(X[0])
            X = np.array(X)

        if mode == 'kmeans':
            k_means = KMeans(n_clusters=n_clusters, random_state=10)
            k_means.fit(X)
            y_predict = k_means.predict(X)
            cluster_center = sum(k_means.cluster_centers_.tolist(),[])
            cluster_res = k_means.labels_.tolist()
        elif mode == 'GaussianMixture':
            clst = mixture.GaussianMixture(n_components=n_clusters, max_iter=100,
                                           covariance_type='full')
            clst.fit(X)
            predicted_labels = clst.predict(X)
            cluster_center = sum(clst.means_.tolist(), [])
            cluster_res = predicted_labels.tolist()

        ms_dis_dict = {}
        for i in range(n_clusters):
            ms_dis_dict[round(cluster_center[i])] = cluster_res.count(i)
        self.ms_dis = ms_dis_dict

    def call_micro(self):

        self.ref_str_ms = pysam.FastaFile(self.reference).fetch(self.chrom, self.mut_start, self.mut_end + 1)
        if not self.call_init():
            return

        if self.reads_phased:
            dis_dict_temp = {}
            for key,value in self.ms_dis_hap1.items():
                dis_dict_temp[round(key/self.repeat_unit_len,2)] = value
            self.ms_dis_hap1 = dis_dict_temp
            dis_dict_temp = {}
            for key,value in self.ms_dis_hap2.items():
                dis_dict_temp[round(key/self.repeat_unit_len,2)] = value
            self.ms_dis_hap2 = dis_dict_temp

            res1 = get_repeat_gmm(self.ms_dis_hap1, target=1)
            res2 = get_repeat_gmm(self.ms_dis_hap2, target=1)

            hap1_repeat_length = res1["genotype"][0]
            hap2_repeat_length = res2["genotype"][0]
            self.qual_ms_hap1 = np.round(res1["qual"], )
            self.qual_ms_hap2 = np.round(res2["qual"], 2)
            self.qual_ms = np.round(np.mean([res1["qual"], res2["qual"]]), 2)

            if hap1_repeat_length in self.ms_dis_hap1.keys() \
                    and hap2_repeat_length in self.ms_dis_hap2.keys():
                self.support_reads_mutNum = [self.ms_dis_hap1[hap1_repeat_length],self.ms_dis_hap2[hap2_repeat_length]]
        else:
            dis_dict_temp = {}
            for key,value in self.ms_dis.items():
                dis_dict_temp[round(key/self.repeat_unit_len,2)] = value
            self.ms_dis = dis_dict_temp

            res = get_repeat_gmm(self.ms_dis, target=2)
            hap1_repeat_length, hap2_repeat_length = res["genotype"]
            self.qual_ms = np.round(res["qual"], 2)

            if hap1_repeat_length in self.ms_dis.keys() and hap2_repeat_length in self.ms_dis.keys():
                self.support_reads_mutNum = [self.ms_dis[hap1_repeat_length],self.ms_dis[hap2_repeat_length]]

        if hap1_repeat_length == hap2_repeat_length:
            if hap1_repeat_length == self.repeat_len:
                self.format_GT_ms = (0, 0)

            else:
                self.format_GT_ms = (1, 1)
                self.alt_ms = (str(hap1_repeat_length) + "[" + self.repeat_unit + "]",)

                reads_num=len(self.hmmRes_draw[self.ms_id]['predict_res'])
                chrName=self.ms_id.split('_')[0]
                temp1_str=chrName+':'+str(self.start)+'-'+str(self.end)+' '+self.repeat_unit

                fig, ax = plt.subplots()
                bars1 = ax.barh(np.arange(reads_num), width=self.end-self.start, color='skyblue')
                ax.set_xticks([])
                ax.set_yticks([])
                ax.spines['top'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['left'].set_visible(False)
                plt.text(0, 0, temp1_str, fontsize=8, color='b', transform=ax.transAxes)
                plt.text(0.01, 1, str(self.start), fontsize=8, color='black', transform=ax.transAxes)
                plt.text(0.85, 1, str(self.end), fontsize=8, color='black', transform=ax.transAxes)

                temp1_num=0
                for bar in bars1:
                    bar_x, bar_y = bar.get_xy()
                    bar_width = bar.get_width()
                    bar_height = bar.get_height()

                    predict_res=self.hmmRes_draw[self.ms_id]['predict_res'][temp1_num]

                    if predict_res==[]:continue

                    for ii in range(len(predict_res)):
                        indel = predict_res[ii]
                        if 'I' in indel:
                            nums = indel.count('I')
                            gap = 1
                            rectangle = patches.Rectangle((bar_x + ii, bar_y), gap, bar_height * 1,
                                                          facecolor='r', linestyle='None')

                            center_x = rectangle.get_x() + 0.25 * rectangle.get_width()
                            center_y = rectangle.get_y() + 0.5 * rectangle.get_height()
                            text_fontsize = min(rectangle.get_width(), bar_height) * 12 * 0.5
                            plt.text(center_x, center_y, str(nums), color='w', fontsize=text_fontsize, weight='bold')

                            ax.add_patch(rectangle)
                        elif indel == 'D':
                            gap = 1
                            rectangle = patches.Rectangle((bar_x + ii, bar_y), gap, bar_height * 1,
                                                          facecolor='w', linestyle='None')
                            ax.add_patch(rectangle)
                            line = patches.ConnectionPatch((bar_x + ii, bar_y + bar_height * 0.5),
                                                           (bar_x + ii + gap, bar_y + bar_height * 0.5),
                                                           coordsA='data', coordsB='data', linewidth=bar_height * 0.4,
                                                           color='skyblue', arrowstyle='-')
                            ax.add_patch(line)

                    temp1_num+=1

                ifVisual = get_value("paras")['mut_information_visual']
                outPath=get_value("paras")["output"]
                if not os.path.exists(str(outPath)+str('res_picture/')):
                    os.makedirs(str(outPath)+str('res_picture/'))

                if ifVisual=='True':
                    plt.savefig(str(outPath)+str('res_picture/%s.jpg' % (self.ms_id)),dpi=300)
                plt.close()
        else:
            if self.repeat_len in [hap1_repeat_length, hap2_repeat_length]:
                if self.reads_phased:
                    if hap1_repeat_length == self.repeat_len:
                        self.format_GT_ms = (0, 1)
                        self.alt_ms = (str(hap2_repeat_length) + "[" + self.repeat_unit + "]",)
                    else:
                        self.alt_ms = (str(hap1_repeat_length) + "[" + self.repeat_unit + "]",)
                        self.format_GT_ms = (1, 0)
                else:
                    self.format_GT_ms = (0, 1)
                    if hap1_repeat_length == self.repeat_len:
                        self.alt_ms = (str(hap2_repeat_length) + "[" + self.repeat_unit + "]",)
                    else:
                        self.alt_ms = (str(hap1_repeat_length) + "[" + self.repeat_unit + "]",)
            else:
                if self.reads_phased:
                    if hap1_repeat_length < hap2_repeat_length:
                        self.alt_ms = (str(hap1_repeat_length) + "[" + self.repeat_unit + "]",
                                       str(hap2_repeat_length) + "[" + self.repeat_unit + "]")
                        self.format_GT_ms = (1, 2)
                    else:
                        self.alt_ms = (str(hap2_repeat_length) + "[" + self.repeat_unit + "]",
                                       str(hap1_repeat_length) + "[" + self.repeat_unit + "]")
                        self.format_GT_ms = (2, 1)
                else:
                    self.format_GT_ms = (1, 2)
                    if hap1_repeat_length < hap2_repeat_length:
                        self.alt_ms = (str(hap1_repeat_length) + "[" + self.repeat_unit + "]",
                                       str(hap2_repeat_length) + "[" + self.repeat_unit + "]")
                    else:
                        self.alt_ms = (str(hap2_repeat_length) + "[" + self.repeat_unit + "]",
                                       str(hap1_repeat_length) + "[" + self.repeat_unit + "]")

            reads_num=len(self.hmmRes_draw[self.ms_id]['predict_res'])
            chrName=self.ms_id.split('_')[0]
            temp1_str=chrName+':'+str(self.start)+'-'+str(self.end)+' '+self.repeat_unit

            fig, ax = plt.subplots()
            bars1 = ax.barh(np.arange(reads_num), width=self.end-self.start, color='skyblue')
            ax.set_xticks([])
            ax.set_yticks([])
            ax.spines['top'].set_visible(False)
            ax.spines['bottom'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            plt.text(0, 0, temp1_str, fontsize=8, color='b', transform=ax.transAxes)
            plt.text(0.01, 1, str(self.start), fontsize=8, color='black', transform=ax.transAxes)
            plt.text(0.85, 1, str(self.end), fontsize=8, color='black', transform=ax.transAxes)

            temp1_num=0
            for bar in bars1:
                bar_x, bar_y = bar.get_xy()
                bar_width = bar.get_width()
                bar_height = bar.get_height()

                predict_res=self.hmmRes_draw[self.ms_id]['predict_res'][temp1_num]

                if predict_res==[]:continue

                for ii in range(len(predict_res)):
                    indel = predict_res[ii]
                    if 'I' in indel:
                        nums = indel.count('I')
                        gap = 1
                        rectangle = patches.Rectangle((bar_x + ii, bar_y), gap, bar_height * 1,
                                                      facecolor='r',linestyle='None')

                        center_x = rectangle.get_x() + 0.25 * rectangle.get_width()
                        center_y = rectangle.get_y() + 0.5 * rectangle.get_height()
                        text_fontsize = min(rectangle.get_width(), bar_height) * 12 *0.5
                        plt.text(center_x, center_y, str(nums), color='w', fontsize=text_fontsize, weight='bold')

                        ax.add_patch(rectangle)
                    elif indel=='D':
                        gap = 1
                        rectangle = patches.Rectangle((bar_x + ii, bar_y), gap, bar_height * 1,
                                                      facecolor='w', linestyle='None')
                        ax.add_patch(rectangle)
                        line = patches.ConnectionPatch((bar_x + ii, bar_y + bar_height * 0.5),
                                                       (bar_x + ii + gap, bar_y + bar_height * 0.5),
                                                       coordsA='data', coordsB='data', linewidth= bar_height*0.4,
                                                       color='skyblue', arrowstyle='-')
                        ax.add_patch(line)

                temp1_num+=1

            ifVisual = get_value("paras")['mut_information_visual']
            outPath = get_value("paras")["output"]
            if not os.path.exists(str(outPath) + str('res_picture/')):
                os.makedirs(str(outPath) + str('res_picture/'))

            if ifVisual == 'True':
                plt.savefig(str(outPath) + str('res_picture/%s.jpg' % (self.ms_id)),dpi=300)
            plt.close()

        if self.reads_phased:
            self.format_AL_ms = "/".join(list(map(str, [hap1_repeat_length, hap2_repeat_length])))
            self.format_DP_ms = "/".join(list(map(str, [self.support_hap0, self.support_hap1, self.support_hap2])))
            self.format_QL_ms = "/".join(list(map(str, [self.qual_ms, self.qual_ms_hap1, self.qual_ms_hap2])))
            if hap1_repeat_length != self.repeat_len: self.hap1_ms_mut = True
            if hap2_repeat_length != self.repeat_len: self.hap2_ms_mut = True

        else:
            if hap1_repeat_length <= hap2_repeat_length:
                self.format_AL_ms = "/".join(list(map(str, [hap1_repeat_length, hap2_repeat_length])))
                self.format_DP_ms = "/".join(
                    list(map(str, [self.support_hap0, self.support_hap1, self.support_hap2])))
                self.format_QL_ms = "/".join(list(map(str, [self.qual_ms, self.qual_ms_hap1, self.qual_ms_hap2])))
            else:
                self.format_AL_ms = "/".join(list(map(str, [hap2_repeat_length, hap1_repeat_length])))
                self.format_DP_ms = "/".join(
                    list(map(str, [self.support_hap0, self.support_hap2, self.support_hap1])))
                self.format_QL_ms = "/".join(list(map(str, [self.qual_ms, self.qual_ms_hap2, self.qual_ms_hap1])))

        if sum(self.format_GT_ms) > 0:
            self.ms_mutation = True
        else:
            self.ms_mutation = False

    def call_micro_and_other(self):
        self.call_micro()
        self.deletion_merge()
        self.build_patterns()

    def one_hap_genotype(self):
        if self.depth == 0:
            self.check = False
            self.check_status.append("No_read_covered")
            self.dis_stat = False
            self.ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, self.mut_start, self.mut_end + 1)
            return
        self.deletion_merge()
        mismatches = {}
        deletions = {}
        insertions = {}
        ms_dis = {}
        mut_start = self.start
        mut_end = self.end - 1
        for read_id, read_info in self.muts.items():
            for mut in read_info.mismatches:
                mut_start = min(mut_start, mut[0])
                mut_end = max(mut_end, mut[0])
                if mut[0] not in mismatches:
                    mismatches[mut[0]] = {}
                if mut[1] not in mismatches[mut[0]]:
                    mismatches[mut[0]][mut[1]] = 1
                else:
                    mismatches[mut[0]][mut[1]] += 1
            for mut in read_info.insertions:
                mut_start = min(mut_start, mut[0])
                mut_end = max(mut_end, mut[0])
                if mut[0] not in insertions:
                    insertions[mut[0]] = {}
                if mut[1] not in insertions[mut[0]]:
                    insertions[mut[0]][mut[1]] = 1
                else:
                    insertions[mut[0]][mut[1]] += 1
            for mut in read_info.deletions:
                if mut[0] not in deletions:
                    deletions[mut[0]] = {}
                if mut[1] not in deletions[mut[0]]:
                    deletions[mut[0]][mut[1]] = 1
                else:
                    deletions[mut[0]][mut[1]] += 1
                mut_start = min(mut_start, mut[0])
                mut_end = max(mut_end, mut[0] + mut[1])

            if read_info.repeat_length not in ms_dis:
                ms_dis[read_info.repeat_length] = 1
            else:
                ms_dis[read_info.repeat_length] += 1
        self.mut_start = mut_start
        self.mut_end = mut_end
        self.ms_dis = ms_dis
        self.query_repeat_length = get_max_support_index(ms_dis)
        for mut in deletions.values():
            if len(mut) > 1:
                self.check = False
                self.check_status.append("More_alleles_in_MS")
        for mut in insertions.values():
            if len(mut) > 1:
                self.check = False
                self.check_status.append("More_alleles_in_MS")
        for mut in mismatches.values():
            if len(mut) > 1:
                self.check = False
                self.check_status.append("More_alleles_in_MS")

        alt_list = []
        for mut_pos, info in mismatches.items():
            gt_list = list(info.keys())
            alt_list.append([mut_pos, gt_list[0]])
            if mut_pos < self.start:
                self.mut.var_pre[mut_pos] = ["SNV", gt_list]
            elif mut_pos < self.end:
                self.mut.var_ms[mut_pos] = ["SNV", gt_list]
            else:
                self.mut.var_suf[mut_pos] = ["SNV", gt_list]
        for mut_pos, info in insertions.items():
            gt_list = list(info.keys())
            alt_list.append([mut_pos, gt_list[0]])
            if mut_pos < self.start - 1:
                self.mut.var_pre[mut_pos] = ["INS", gt_list]
            elif mut_pos < self.end:
                self.mut.var_ms[mut_pos] = ["INS", gt_list]
            else:
                self.mut.var_suf[mut_pos] = ["INS", gt_list]
        for mut_pos, info in deletions.items():
            gt_list = list(info.keys())
            for del_pos in range(mut_pos, mut_pos + gt_list[0]):
                alt_list.append([del_pos, ""])
            del_end = mut_pos + gt_list[0] - 1
            if del_end < self.start:
                self.mut.var_pre[mut_pos] = ["DEL", gt_list]
            elif mut_pos >= self.end:
                self.mut.var_suf[mut_pos] = ["DEL", gt_list]
            else:
                if del_end < self.end or mut_pos >= self.start:
                    self.mut.var_suf[mut_pos] = ["DEL", gt_list]
                else:
                    self.check = False
                    self.check_status.append("DEL_Span")
        self.mut.compute()
        self.ref_str = pysam.FastaFile(self.reference).fetch(self.chrom, self.mut_start - 1, self.mut_end + 1)
        self.alt_str = "." if self.mut.var_type == "None" else self.get_alt(alt_list, offset=self.mut_start - 1)

    def get_mut_one_hap(self, hap2, hap2_ms_mut):
        pass


