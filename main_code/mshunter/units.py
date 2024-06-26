#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np
import time
import pandas as pd
from global_dict import *


class Read_Mutation:

    def __init__(self, read_id, repeat_length, strand, mismatches, deletions, insertions, hap, prefix="", suffix="",
                 ms_content="", pos_based_info={}, read_str=""):
        self.read_id = read_id
        self.repeat_length = repeat_length
        self.strand = strand
        self.mismatches = mismatches
        self.deletions = deletions
        self.insertions = insertions
        self.hap = hap
        self.prefix = prefix
        self.suffix = suffix
        self.ms_content = ms_content
        self.pos_based_info = pos_based_info
        self.del_span = "False"
        self.before = False
        self.after = False
        self.ms = False
        self.read_str = read_str
        self.other_mutation = True if len(pos_based_info) > 0 else False

    def get_mut_by_read(self):

        mut_type = set([info[0] for pos, info in self.pos_based_info.items()])
        mut_start = [pos for pos, info in self.pos_based_info.items()]
        mut_end = []
        for pos, info in self.pos_based_info.items():
            if info[0] == "SNV" or info[0] == "INS":
                mut_end.append(pos)
            else:
                mut_end.append(pos + len(info[1]))
        if len(mut_start) > 0:
            return min(mut_start), max(mut_end), mut_type
        else:
            return -1, -1, mut_type


class PatternCluster:
    reads_id = []
    pattern = {}
    pattern_num = 0
    read_num = 0
    no_mutation = False
    pattern_id = 0
    forward = []
    reversed = []
    hap0 = []
    hap1 = []
    hap2 = []

    def __init__(self, pattern_id):
        self.pattern_id = pattern_id

    def contain(self, pattern={}):
        match_num = 0
        for pos, ann in pattern.items():
            if pos in self.pattern:
                if self.pattern[pos] == ann:
                    match_num += 1

        match_ratio = (match_num / self.pattern_num + match_num / len(pattern)) * 0.5
        return True if match_ratio == 1 else False

    def init_patttern(self, pattern, read_id, hap, strand):
        self.pattern = pattern
        self.reads_id = [read_id]
        self.pattern_num = len(pattern)
        self.read_num = 1
        if self.pattern_num == 0:
            self.no_mutation = True
        if hap == 0: self.hap0.append(read_id)
        if hap == 1: self.hap1.append(read_id)
        if hap == 2: self.hap2.append(read_id)
        if strand:
            self.forward.append(read_id)
        else:
            self.reversed.append(read_id)

        pass

    def add_reads(self, read_id, hap, strand):
        self.reads_id.append(read_id)
        self.read_num += 1
        if hap == 0: self.hap0.append(read_id)
        if hap == 1: self.hap1.append(read_id)
        if hap == 2: self.hap2.append(read_id)
        if strand:
            self.forward.append(read_id)
        else:
            self.reversed.append(read_id)


class Mutation:

    def __init__(self):
        self.var_pre = {}
        self.var_suf = {}
        self.var_ms = {}
        self.type_pre = []
        self.type_suf = []
        self.type_ms = []
        self.mut_start = 0
        self.mut_end = 0
        self.alt = []
        self.ms = 0
        self.pre = 0
        self.suf = 0
        self.var_detail = []
        self.var_type = "None"
        self.var_type_detail = ""
        self.var_detail_str = ""

    def compute(self, ):
        for pos in sorted(list(self.var_pre.keys())):
            var_type = self.var_pre[pos][0]
            var_content = self.var_pre[pos][1][0]
            self.type_pre.append(self.var_pre[pos][0])
            self.alt.append([pos, self.var_pre[pos][1][0]])
            self.var_detail.append([str(pos), var_type, str(var_content)])

        for pos in sorted(list(self.var_ms.keys())):
            var_type = self.var_ms[pos][0]
            var_content = self.var_ms[pos][1][0]
            self.type_ms.append(var_type)
            self.alt.append([pos, var_content])
            self.var_detail.append([str(pos), var_type, str(var_content)])

        for pos in sorted(list(self.var_suf.keys())):
            var_type = self.var_suf[pos][0]
            var_content = self.var_suf[pos][1][0]
            self.type_suf.append(self.var_suf[pos][0])
            self.alt.append([pos, self.var_suf[pos][1]])
            self.var_detail.append([str(pos), var_type, str(var_content)])

        self.ms = len(self.type_ms)
        self.pre = len(self.type_pre)
        self.suf = len(self.type_suf)
        var_list = self.type_pre + self.type_ms + self.type_suf
        var_num = len(var_list)
        if var_num > 1:
            self.var_type = "Complex"
        elif var_num == 0:
            self.var_type = "None"
        else:
            self.var_type = var_list[0]
        self.var_type_detail = "|".join(
            [":".join(var_type) for var_type in [self.type_pre, self.type_ms, self.type_suf]])
        self.var_detail_str = "|".join([":".join(var) for var in self.var_detail])

    def show(self):
        if len(self.var_ms) > 0:
            print(self.var_ms)
        if len(self.var_pre) > 0:
            print(self.var_pre)
        if len(self.var_suf) > 0:
            print(self.var_suf)


def remove_zero_dict(dict):
    newdict = {}
    for key, value in dict.items():
        if value > 0.000001:
            newdict[key] = value
    return newdict


def getDisdistance(dict1, dict2):
    dictkey = list(dict1.keys()) + list(dict2.keys())
    for key in dictkey:
        if key not in dict1:
            dict1[key] = 0
        if key not in dict2:
            dict2[key] = 0
    sum = 0
    for key in dictkey:
        err = dict1[key] - dict2[key]
        sum += (err * err)
    return round(np.sqrt(sum), 6)


def get_disdistance(dict1, dict2):
    dictkey = list(dict1.keys()) + list(dict2.keys())
    for key in dictkey:
        if key not in dict1:
            dict1[key] = 0
        if key not in dict2:
            dict2[key] = 0
    sum = 0
    for key in dictkey:
        err = dict1[key] - dict2[key]
        sum += (err * err)
    return round(np.sqrt(sum), 6)


def getDisdistance2(dict1, dict2):
    dictkey = list(dict1.keys()) + list(dict2.keys())
    list1 = []
    list2 = []
    for key in dictkey:
        if key not in dict1:
            list1.append(0)
        else:
            list1.append(dict1[key])
        if key not in dict2:
            list2.append(0)
        else:
            list2.append(dict2[key])
    return np.linalg.norm(np.array(list1) - np.array(list2))


def getDisdistance3(dict1, dict2):
    dictkey = list(dict1.keys()) + list(dict2.keys())
    list1 = []
    list2 = []
    for key in dictkey:
        if key not in dict1:
            list1.append(0)
        else:
            list1.append(dict1[key])
        if key not in dict2:
            list2.append(0)
        else:
            list2.append(dict2[key])
    return np.sqrt(np.sum(np.square(np.array(list1) - np.array(list2))))


def load_microsatellites(args):

    logger.info("Loading microsatellite file...")
    ms = args["microsatellite"]
    separator = args["microsatellite_region_format"]
    if separator == "comma":
        df_microSatellites = pd.read_csv(ms, index_col=0)
    elif separator == "msisensor_scan":
        df_microSatellites = pd.read_table(ms, header=0)
        columns = df_microSatellites.columns
        if "chromosome" in columns:
            df_microSatellites.rename(columns={"chromosome": "chr"}, inplace=True)
        if "location" in columns:
            df_microSatellites.rename(columns={"location": "pos"}, inplace=True)
        if "repeat_unit_bases" in columns:
            df_microSatellites.rename(columns={"repeat_unit_bases": "motif"}, inplace=True)
        if "repeat_unit_length" in columns:
            df_microSatellites.rename(columns={"repeat_unit_length": "motifLen"}, inplace=True)
        if "repeat_times" in columns:
            df_microSatellites.rename(columns={"repeat_times": "repeatTimes"}, inplace=True)
        if "left_flank_bases" in columns:
            df_microSatellites.rename(columns={"left_flank_bases": "prefix"}, inplace=True)
        if "right_flank_bases" in columns:
            df_microSatellites.rename(columns={"right_flank_bases": "suffix"}, inplace=True)
    elif separator == "space":
        df_microSatellites = pd.read_table(ms, header=0, sep=" ")
    chromList = get_value("chrom_list")
    df_microSatellites = df_microSatellites[df_microSatellites['chr'].isin(chromList)]
    if args["only_homopolymer"]:
        df_microSatellites = df_microSatellites[df_microSatellites['motifLen'] == 1]

    repeatRange = args["ranges_of_repeat_times"]
    repeatUnitList = sorted(repeatRange.keys())

    df_microsatellite_pass = pd.DataFrame()
    for ul in repeatUnitList:
        minr = repeatRange[ul]["min"]
        maxr = repeatRange[ul]["max"]
        df_microsatellite_pass = pd.concat(
            [df_microsatellite_pass, df_microSatellites[(df_microSatellites["motifLen"] == ul) &
                                                        (df_microSatellites["repeatTimes"] >= minr) &
                                                        (df_microSatellites["repeatTimes"] <= maxr)
                                                        ]])

    logger.info("There are total " + str(len(df_microsatellite_pass)) + " microsatellites(bed文件位点总数).")
    set_value("ms_number", len(df_microsatellite_pass))
    set_value("motif_list", set(df_microsatellite_pass["motif"]))
    return df_microsatellite_pass

def get_max_support_index(input_dict):

    m = max(input_dict.keys(), key=(lambda x: input_dict[x]))
    return m

def dis_stats(dis, values=["mean", "std"]):
    repeat_length_list = []
    for repeat_len, times in dis.items():
        repeat_length_list.extend([repeat_len] * times)
    res_list = []
    for value in values:
        if value == "mean":
            this_value = np.mean(repeat_length_list)
        elif value == "std":
            this_value = np.std(repeat_length_list)
        elif value == "max":
            this_value = np.max(repeat_length_list)
        elif value == "min":
            this_value = np.min(repeat_length_list)
        else:
            this_value = None
        res_list.append(this_value)
    return res_list

def str2int(item):
    if len(item) > 0:
        return np.nanmean([ord(i) - 33 for i in item])
    else:
        return np.nan

def int2str(qual_int):
    if np.isnan(qual_int):
        return chr(33)
    else:
        return chr(int(qual_int) + 33)

def dis_sum(dict_list):
    keylist = []
    for item in dict_list:
        for key in item:
            if key not in keylist:
                keylist.append(key)
    res_dict = {}
    for key in keylist:
        res_dict[key] = 0
    for item in dict_list:
        for key in item:
            res_dict[key] += item[key]
    return res_dict

def change_dim_for_pool_map(input, threads):
    item_num = 0
    for win in input:
        item_num += len(win)
    item_per_win = item_num // threads + 1
    item_index = 0
    output = []
    win_sub = []
    for win in input:
        for item in win:
            item_index += 1
            win_sub.append(win_sub)
            if item_index % item_per_win == 0:
                output.append(win_sub)
                win_sub = []
    if len(win_sub) > 0:
        output.append(win_sub)
    return output

if __name__ == "__main__":
    a = get_max_support_index({1: 5, 6: 40, 3: 2})
    print(get_max_support_index({13: 8, 14: 11, 12: 1}))

