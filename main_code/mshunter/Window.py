#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
sys.path.append('/home/zyx/miniconda3/envs/mshunter_env/lib/python3.7/site-packages')

import multiprocessing
from global_dict import *
from Microsatellite import Microsatellite
from Read import Read
import pysam
import pandas as pd
import numpy as np
from hmmlearn import hmm
from collections import *
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
import os
import pandas as pd
import pickle

class Window:

    def __init__(self, ms_info_list, tech=""):
        self.tech = tech
        self.contig = ms_info_list[0]["chr"]
        self.paras = get_value("paras")
        self.ms_list = ms_info_list
        self.bam_path = self.paras["input"]
        self.win_start = ms_info_list[0]["pos"] - self.paras["prefix_len"]
        self.win_end = ms_info_list[-1]["pos"] + ms_info_list[-1]["repeatTimes"] * ms_info_list[-1]["motifLen"] + \
                       self.paras["suffix_len"]
        self.reads = {}
        self.reads_num = 0
        self.microsatellites = {}
        self.microsatellites_id = [it["chr"] + "_" + str(it["pos"]) for it in ms_info_list]
        self.vcf_recs = []

        self.minimum_support_reads_ratio = get_value("paras")["minimum_support_reads_ratio"]

    def init_microsatellites(self, only_simple=True):

        microsatellites = []
        for ms in self.ms_list:
            microsatellites.append(Microsatellite(ms, only_simple=only_simple))
        self.microsatellites = {ms_info.ms_id: ms_info for ms_info in microsatellites}

    def init_reads(self):

        reads = {}
        sam_file = pysam.AlignmentFile(self.paras["input"], mode="rb", reference_filename=self.paras["reference"])

        for ms_id, ms_info in self.microsatellites.items():
            for alignment in sam_file.fetch(ms_info.chrom, ms_info.start_pre, ms_info.end_suf ):

                if alignment.reference_start > ms_info.end_suf:
                    break

                if alignment.is_unmapped or alignment.is_duplicate or alignment.is_secondary :
                    continue
                k_para = 100

                if alignment.reference_start > ms_info.start_pre - k_para or alignment.reference_end < ms_info.end_suf + k_para:
                    continue
                if len(alignment.query_sequence) < 2:
                    continue

                read_id = alignment.query_name + "_" + str(alignment.reference_start)

                if alignment.query_qualities == None:
                    continue

                if read_id not in reads:
                    reads[read_id] = Read(read_id=read_id,
                                          chrom=self.contig,
                                          alignment=alignment,
                                          reference=self.paras["reference"],
                                          tech=self.tech)
                if ms_info.ms_id not in reads[read_id].support_microsatellites:
                    reads[read_id].support_microsatellites.append(ms_info.ms_id)
        self.reads = reads
        self.reads_num = len(self.reads)

    def get_reads_dis(self):

        result_list = []
        for read in self.reads.values():
            read.microsatellites = {ms_id: self.microsatellites[ms_id] for ms_id in read.support_microsatellites}
            read.get_read_str()
            read.get_repeat_length_all_ms()
            result_list.append(read)
        self.reads = {read.read_id: read for read in result_list}

    def get_reads_info_improveHMM_3(self, only_simple=True):

        result_list = []
        for read in self.reads.values():
            read.microsatellites = {ms_id: self.microsatellites[ms_id] for ms_id in read.support_microsatellites}

            read.get_read_str()

            if only_simple:
                read.get_repeat_length_all_ms()
            else:
                read.get_ms_info_one_read(compar_software=self.paras['compare_software'])

            result_list.append(read)
        self.reads = {read.read_id: read for read in result_list}

    def get_reads_info_improveHMM(self, only_simple=True):

        start_pi,trainsition_probability,model=self.train_HMM()

        result_list = []
        for read in self.reads.values():
            read.microsatellites = {ms_id: self.microsatellites[ms_id] for ms_id in read.support_microsatellites}

            read.get_read_str()

            if only_simple:
                read.get_repeat_length_all_ms()
            else:
                read.get_ms_info_one_read_improved(model,compar_software=self.paras['compare_software'])

            result_list.append(read)
        self.reads = {read.read_id: read for read in result_list}

    def get_reads_info_improveHMM_2(self, only_simple=True):

        result_list = []
        for read in self.reads.values():
            read.microsatellites = {ms_id: self.microsatellites[ms_id] for ms_id in read.support_microsatellites}

            read.get_read_str_2()
            result_list.append(read)

        str_varia_list,str_cigar_list,str_varia_toNum=[],[],[]
        for read in self.reads.values():
            str_varia_list.extend(read.get_all_str_varia(flankLen_left=1,flankLen_right=0)[0])
            str_cigar_list.extend(read.get_all_str_varia(flankLen_left=1,flankLen_right=0)[1])
            str_varia_toNum.extend(read.get_all_str_varia(flankLen_left=0,flankLen_right=0)[2])

        # model=self.train_HMM_2()
        # model=self.train_HMM_3(str_varia_list,str_cigar_list,str_varia_toNum,flankLen_left=1,flankLen_right=0)
        model=self.train_HMM_4(str_varia_toNum,stateNum=3)
        # model=self.train_HMM_4_useModel(str_varia_toNum,stateNum=3)
        # model=self.train_HMM_4_saveModel(str_varia_toNum,stateNum=3)
        # model=self.train_HMM_4(str_varia_toNum,stateNum=12)

        # min_train_dataSet=100
        #
        # X_len=[]
        # for i in str_varia_toNum:X_len.append(len(i))
        # X_len_y=int(list(dict(Counter(X_len)).keys())[0])
        # str_varia_toNum_len=[]
        # for i in str_varia_toNum:
        #     if len(i)==X_len_y:
        #         str_varia_toNum_len.append(i)
        # X = np.array(str_varia_toNum_len)
        # print(len(X))
        # if len(X)<min_train_dataSet:
        #     model = self.train_HMM_4_useModel(str_varia_toNum, stateNum=3)
        # else:
        #     model = self.train_HMM_4(str_varia_toNum, stateNum=3)

        # min_saveModel_dataSet=100
        # min_train_dataSet=100
        # X_len=[]
        # for i in str_varia_toNum:X_len.append(len(i))
        # # if len(X_len)==0:
        # #     with open('/mnt/d/science/simulateData/remakeData/STRtool/MSHunter/model_hmm/model1_12.pkl', 'rb') as f:
        # #         model = pickle.load(f)
        # # else:
        # X_len_y=int(list(dict(Counter(X_len)).keys())[0])
        # str_varia_toNum_len=[]
        # for i in str_varia_toNum:
        #     if len(i)==X_len_y:
        #         str_varia_toNum_len.append(i)
        # X = np.array(str_varia_toNum_len)
        # print(len(X))
        # if len(X)>min_saveModel_dataSet:
        #     model = self.train_HMM_4_saveModel_2(str_varia_toNum,stateNum=3)
        # elif len(X)>min_train_dataSet:
        #     model = self.train_HMM_4(str_varia_toNum, stateNum=3)
        # else:
        #     model = self.train_HMM_4_useModel(str_varia_toNum, stateNum=3)

        model_gmm=self.train_GMM(str_varia_list,str_cigar_list,flankLen_left=1,flankLen_right=0)
        # model_gmm=self.train_GMM_2(str_varia_list,str_cigar_list,flankLen_left=1,flankLen_right=0)
        # model_gmm={}
        # print(model_gmm)

        for read in self.reads.values():
            # read.get_ms_info_one_read_improved_2(model, compar_software=self.paras['compare_software'])
            read.get_ms_info_one_read_improved_2(model, model_gmm,compar_software=self.paras['compare_software'])

        self.reads = {read.read_id: read for read in result_list}

    def get_reads_info_improveHMM_4(self, only_simple=True):

        result_list = []
        for read in self.reads.values():
            read.microsatellites = {ms_id: self.microsatellites[ms_id] for ms_id in read.support_microsatellites}

            read.get_read_str()

            if only_simple:
                read.get_repeat_length_all_ms()
            else:
                read.get_ms_info_one_read_improved_3(compar_software=self.paras['compare_software'])

        # for read in self.reads.values():
        #     read.get_ms_info_one_read_improved_3(compar_software=self.paras['compare_software'])

            result_list.append(read)

        self.reads = {read.read_id: read for read in result_list}

    def get_reads_info(self, only_simple=True):

        result_list = []
        for read in self.reads.values():
            read.microsatellites = {ms_id: self.microsatellites[ms_id] for ms_id in read.support_microsatellites}

            read.get_read_str()

            if only_simple:
                read.get_repeat_length_all_ms()
            else:
                read.get_ms_info_one_read(compar_software=self.paras['compare_software'])

            result_list.append(read)
        self.reads = {read.read_id: read for read in result_list}

    def train_HMM(self):

        X=[]
        for ms_id, ms_info in self.microsatellites.items():
            X_temp = []
            seq = ms_info.prefix + ms_info.repeat_unit * ms_info.repeat_times + ms_info.suffix
            if len(seq) != 20: continue
            for i in range(len(seq)):
                if seq[i] == 'A':
                    X_temp.append(0)
                elif seq[i] == 'T':
                    X_temp.append(1)
                elif seq[i] == 'C':
                    X_temp.append(2)
                elif seq[i] == 'G':
                    X_temp.append(3)
            X.append(X_temp)
        print('训练数据量：', len(X))

        X = np.array(X)

        states = ["match", "deletion", "insertion"]
        n_states = len(states)

        model = hmm.GaussianHMM(n_components=n_states, n_iter=1000, tol=0.01, covariance_type='diag', min_covar=0,
                                startprob_prior=1, means_prior=0, means_weight=0, random_state=10)
        model.fit(X)

        return model.startprob_, model.transmat_, model

    def train_HMM_2(self,flankLen=0):

        reads_observe,cigars_state=[],[]
        for read in self.reads.values():
            for strr in read.this_read_list:
                if strr=='A':
                    reads_observe.append(0)
                elif strr=='T':
                    reads_observe.append(1)
                elif strr=='C':
                    reads_observe.append(2)
                elif strr=='G':
                    reads_observe.append(3)
            print(len(read.this_ref_str),len(read.cigars_list),len(read.this_read_list))
            for i in range(len(read.cigars_list)):
                nums=read.cigars_list[i]
                if (nums==2 or nums==3) and read.this_ref_str[i]=='A':
                    cigars_state.append(0)
                elif (nums==2 or nums==3) and read.this_ref_str[i]=='T':
                    cigars_state.append(1)
                elif (nums==2 or nums==3) and read.this_ref_str[i]=='C':
                    cigars_state.append(2)
                elif (nums==2 or nums==3) and read.this_ref_str[i]=='G':
                    cigars_state.append(3)
                elif (nums==1) and read.this_ref_str[i]=='A':
                    cigars_state.append(4)
                elif (nums==1) and read.this_ref_str[i]=='T':
                    cigars_state.append(5)
                elif (nums==1) and read.this_ref_str[i]=='C':
                    cigars_state.append(6)
                elif (nums==1) and read.this_ref_str[i]=='G':
                    cigars_state.append(7)
                elif (nums==0 or nums==7 or nums==8) and read.this_ref_str[i]=='A':
                    cigars_state.append(8)
                elif (nums==0 or nums==7 or nums==8) and read.this_ref_str[i]=='T':
                    cigars_state.append(9)
                elif (nums==0 or nums==7 or nums==8) and read.this_ref_str[i]=='C':
                    cigars_state.append(10)
                elif (nums==0 or nums==7 or nums==8) and read.this_ref_str[i]=='G':
                    cigars_state.append(11)

        model=seqlearn.hmm.MultinomialHMM(decode='viterbi', alpha=0.01)
        model.fit(reads_observe,cigars_state,lengths=None)

        return model

    def train_HMM_3(self, str_varia_list,str_cigar_list,str_varia_toNum_list,flankLen_left=1,flankLen_right=0):

        state_in_list=[]
        dict_ATCG,dict_MI={'A':1,'T':2,'C':3,'G':4},{'M':1,'I':2}
        for i in range(len(str_varia_list)):
            str_varia=str_varia_list[i]
            str_cigar=str_cigar_list[i]
            temp=[]
            temp_D_count,temp_D_MI,temp_last_ATCG=0,'',''
            for j in range(1,len(str_varia)):
                if str_cigar[j]=='M' or str_cigar[j]=='I':
                    if j < len(str_cigar) - 1 and str_cigar[j + 1] == 'D':
                        temp_D_MI = str_cigar[j]
                        temp_last_ATCG=str_varia[j]
                    else:
                        MI_toNum=list(dict_MI.values())[list(dict_MI.keys()).index(str_cigar[j])]
                        ATCG_toNum=list(dict_ATCG.values())[list(dict_ATCG.keys()).index(str_varia[j])]
                        temp.append(int(str(MI_toNum)+str(ATCG_toNum)))
                elif str_cigar[j] == 'D':
                    temp_D_count += 1
                    if j==1:
                        if str_cigar[j-1]=='D':
                            temp_D_MI ='M'
                        else:
                            temp_D_MI=str_cigar[j-1]
                        temp_last_ATCG=str_varia[j-1]

                    if (j < len(str_cigar) - 1 and str_cigar[j + 1] != 'D') or\
                        (j == len(str_cigar) - 1):
                        MI_toNum=list(dict_MI.values())[list(dict_MI.keys()).index(temp_D_MI)]
                        ATCG_toNum=list(dict_ATCG.values())[list(dict_ATCG.keys()).index(temp_last_ATCG)]
                        temp.append(int(str(temp_D_count)+'1'+str(MI_toNum)+str(ATCG_toNum)))
                        temp_D_count, temp_D_MI, temp_last_ATCG = 0, '', ''
            state_in_list.append(temp)

        model = pg.HiddenMarkovModel.from_samples(distribution=pg.MultivariateGaussianDistribution, n_components=3,
                                               X=str_varia_toNum_list,
                                               labels=state_in_list,
                                               algorithm='labeled')

        return model

    def train_HMM_4(self, str_varia_toNum,stateNum):

        if str_varia_toNum==[]:
            return None

        X_len=[]
        for i in str_varia_toNum:
            if len(i)>0:
                X_len.append(len(i))
        X_len_y=int(list(dict(Counter(X_len)).keys())[0])
        str_varia_toNum_len=[]
        for i in str_varia_toNum:
            if len(i)==X_len_y:
                str_varia_toNum_len.append(i)
        X = np.array(str_varia_toNum_len)
        print('训练数据量：', len(X))

        states_3 = ["match", "deletion", "insertion"]
        states_12 = ["matchA","matchT","matchC","matchG", "deletionA","deletionT","deletionC","deletionG",
                     "insertionA","insertionT","insertionC","insertionG"]
        n_states = stateNum

        model = hmm.CategoricalHMM(n_components=n_states, n_iter=1000, tol=0.1, algorithm='viterbi', random_state=42,startprob_prior=1, transmat_prior=1, emissionprob_prior=1)
        model.fit(X)

        # print('score:',model.score(X))

        return model

    def train_HMM_4_saveModel(self, str_varia_toNum,stateNum):

        if str_varia_toNum==[]:
            return None

        X_len=[]
        for i in str_varia_toNum:X_len.append(len(i))
        X_len_y=int(list(dict(Counter(X_len)).keys())[0])
        str_varia_toNum_len=[]
        for i in str_varia_toNum:
            if len(i)==X_len_y:
                str_varia_toNum_len.append(i)

        X = np.array(str_varia_toNum_len)
        print('训练数据量：', len(X))

        states_3 = ["match", "deletion", "insertion"]
        states_12 = ["matchA","matchT","matchC","matchG", "deletionA","deletionT","deletionC","deletionG",
                     "insertionA","insertionT","insertionC","insertionG"]
        n_states = stateNum

        model = hmm.CategoricalHMM(n_components=n_states, n_iter=1000, tol=0.01, algorithm='viterbi', random_state=20,startprob_prior=1, transmat_prior=1, emissionprob_prior=1)
        model.fit(X)
        # print('score:',model.score(X))

        with open('/mnt/d/science/simulateData/remakeData/STRtool/MSHunter/model_hmm/model1.pkl','wb') as f:
            pickle.dump(model,f)

        return model

    def train_HMM_4_saveModel_2(self, str_varia_toNum, stateNum):

        if str_varia_toNum == []:
            return None

        X_len = []
        for i in str_varia_toNum:
            if len(i)==0:continue
            X_len.append(len(i))
        X_len_y = int(list(dict(Counter(X_len)).keys())[0])
        str_varia_toNum_len = []
        for i in str_varia_toNum:
            if len(i) == X_len_y:
                str_varia_toNum_len.append(i)

        X = np.array(str_varia_toNum_len)
        print('训练数据量：', len(X))

        states_3 = ["match", "deletion", "insertion"]
        states_12 = ["matchA", "matchT", "matchC", "matchG", "deletionA", "deletionT", "deletionC", "deletionG",
                     "insertionA", "insertionT", "insertionC", "insertionG"]
        n_states = stateNum

        model = hmm.CategoricalHMM(n_components=n_states, n_iter=1000, tol=0.01, algorithm='viterbi',
                                   random_state=42, startprob_prior=1, transmat_prior=1,
                                   emissionprob_prior=1)
        model.fit(X)

        save_path='/mnt/d/science/simulateData/remakeData/STRtool/MSHunter/model_hmm/model1_12.pkl'
        if os.path.isfile(save_path):
            pass
        else:
            with open(save_path, 'wb') as f:
                pickle.dump(model, f)

        return model

    def train_HMM_4_useModel(self, str_varia_toNum,stateNum):

        states_3 = ["match", "deletion", "insertion"]
        states_12 = ["matchA","matchT","matchC","matchG", "deletionA","deletionT","deletionC","deletionG",
                     "insertionA","insertionT","insertionC","insertionG"]
        n_states = stateNum

        with open('/mnt/d/science/simulateData/remakeData/STRtool/MSHunter/model_hmm/model1_12.pkl','rb') as f:
            model=pickle.load(f)

        return model

    def train_GMM(self,str_varia_list,str_cigar_list,flankLen_left=1,flankLen_right=0):

        return {'AM':[],'AI':[],'TM':[],'TI':[],'CM':[],'CI':[],'GM':[],'GI':[]}
        pass

        dict_temp={'AM':[],'AI':[],'TM':[],'TI':[],'CM':[],'CI':[],'GM':[],'GI':[]}
        for i in range(len(str_varia_list)):

            str_varia=str_varia_list[i]
            str_cigar=str_cigar_list[i]

            if 'N' in str_varia:
                continue

            temp_D_count = 0
            index_v=flankLen_left
            for index_c in range(flankLen_left,len(str_cigar)):
                if str_cigar[flankLen_left-1]=='D':
                    break

                if str_cigar[index_c]!='D':
                    index_v += 1
                    if index_v>=len(str_varia): break

                    if str_cigar[index_c-1]!='D':
                        temp_MI=str_cigar[index_c]
                        temp_ATCG=str_varia[index_v]
                    else:
                        temp_str=temp_ATCG+temp_MI
                        if temp_D_count>0:
                            dict_temp[temp_str].append(temp_D_count)
                        temp_D_count=0
                        temp_MI=str_cigar[index_c]
                        temp_ATCG=str_varia[index_v]
                else:
                    if index_c==len(str_cigar)-1:
                        temp_str = temp_ATCG + temp_MI
                        if temp_D_count > 0:
                            dict_temp[temp_str].append(temp_D_count)
                    elif str_cigar[index_c-1]!='D':
                        temp_MI = str_cigar[index_c-1]
                        temp_ATCG = str_varia[index_v-1]
                        temp_D_count+=1
                    else:
                        temp_D_count += 1

        dict_gmmModel={'AM':[],'AI':[],'TM':[],'TI':[],'CM':[],'CI':[],'GM':[],'GI':[]}
        for key,value in dict_temp.items():
            if len(value)==0:
                value=[1,1] #至
            elif len(value)==1:
                value.append(value[0])

            gmm = GaussianMixture(n_components=1)
            gmm.fit(np.array(value).reshape(-1,1))
            # print(gmm.means_)
            dict_gmmModel[key] =round(gmm.means_[0][0],2)

        return dict_gmmModel

    def train_GMM_2(self, str_varia_list, str_cigar_list, flankLen_left=1, flankLen_right=0):

        dict_temp = {'A': [], 'T': [],'C': [],'G': []}
        for i in range(len(str_varia_list)):
            str_varia = str_varia_list[i]
            str_cigar = str_cigar_list[i]

            temp_D_count = 0
            index_v = flankLen_left
            for index_c in range(flankLen_left, len(str_cigar)):
                if str_cigar[flankLen_left - 1] == 'D':
                    break

                if str_cigar[index_c] != 'D':
                    index_v += 1
                    if index_v >= len(str_varia): break

                    if str_cigar[index_c - 1] != 'D':
                        temp_ATCG = str_varia[index_v]
                    else:
                        temp_str = temp_ATCG
                        if temp_D_count > 0:
                            dict_temp[temp_str].append(temp_D_count)
                        temp_D_count = 0
                        temp_ATCG = str_varia[index_v]
                else:
                    if index_c == len(str_cigar) - 1:
                        temp_str = temp_ATCG
                        if temp_D_count > 0:
                            dict_temp[temp_str].append(temp_D_count)
                    elif str_cigar[index_c - 1] != 'D':
                        temp_ATCG = str_varia[index_v - 1]
                        temp_D_count += 1
                    else:
                        temp_D_count += 1

        dict_gmmModel =  {'A': [], 'T': [],'C': [],'G': []}
        for key, value in dict_temp.items():
            if len(value) == 0:
                value = [1, 1]

            gmm = GaussianMixture(n_components=1)
            gmm.fit(np.array(value).reshape(-1, 1))
            # print(gmm.means_)
            dict_gmmModel[key] = round(gmm.means_[0][0], 2)

        return dict_gmmModel

    def merge_reads_repeat_length_distribution(self):
        microsatellites_dict = {ms_id: {} for ms_id in self.microsatellites}
        for read_id, read in self.reads.items():
            stand = read.strand
            hap = read.hap
            for ms_id, ms_read_repeat_length in read.repeat_lengths.items():
                microsatellites_dict[ms_id][read_id] = [ms_read_repeat_length, stand, hap]
        self.reads = {}
        for ms_id, ms_read_repeat_length_info in microsatellites_dict.items():
            self.microsatellites[ms_id].set_read_dis_info(ms_read_repeat_length_info)

    def merge_reads_info(self):
        microsatellites_dict = {ms_id: {} for ms_id in self.microsatellites}
        for read_id, read in self.reads.items():
            for ms_id in read.microsatellites:
                microsatellites_dict[ms_id][read_id] = read
        self.reads = {}
        for ms_id, reads_info in microsatellites_dict.items():
            self.microsatellites[ms_id].set_reads_info(reads_info)

    def merge_muts_info(self, only_simple=True):
        microsatellites_dict_dis = {ms_id: {} for ms_id in self.microsatellites}
        microsatellites_dict_mut = {ms_id: {} for ms_id in self.microsatellites}

        microsatellites_dict_dis_GmmHmm={ms_id: [] for ms_id in self.microsatellites}

        predictRes={}

        for read_id, read in self.reads.items():
            stand = read.strand
            hap = read.hap
            for ms_id, ms_read_repeat_length in read.repeat_lengths.items():
                microsatellites_dict_dis[ms_id][read_id] = [ms_read_repeat_length, stand, hap]

            for ms_id, ms_read_repeat_length in read.repeat_lengths_GmmHmm.items():
                microsatellites_dict_dis_GmmHmm[ms_id].append(ms_read_repeat_length)

            if not only_simple:
                for ms_id, ms_read_mut in read.mut_info.items():
                    microsatellites_dict_mut[ms_id][read_id] = ms_read_mut

            for ms_id, hmm_item in read.hmmRes_draw.items():
                if ms_id not in predictRes.keys():
                    predictRes[ms_id]=hmm_item
                else:
                    for aa in hmm_item['predict_res']:
                        predictRes[ms_id]['predict_res'].append(aa)

        self.reads = {}

        for ms_id, ms_read_repeat_length_info in microsatellites_dict_dis_GmmHmm.items():
            if len(ms_read_repeat_length_info)==1:
                ms_read_repeat_length_info.append(ms_read_repeat_length_info[0])
            elif len(ms_read_repeat_length_info)==0:
                self.microsatellites[ms_id].GmmHmm_res=0
                continue

            self.microsatellites[ms_id].set_len_info(ms_id,ms_read_repeat_length_info)

        for ms_id, hmm_item in predictRes.items():
            self.microsatellites[ms_id].hmmRes_draw[ms_id]=hmm_item

        for ms_id, ms_read_repeat_length_info in microsatellites_dict_dis.items():
            self.microsatellites[ms_id].set_read_dis_info(ms_read_repeat_length_info)
            if not only_simple:
                self.microsatellites[ms_id].set_muts_info(microsatellites_dict_mut[ms_id])

    def genotype_one_microsatellite_ccs(self, microsatellite):

        return microsatellite

    def call_variants(self):
        microsatellites = []

        for microsatellite in self.microsatellites.values():

            if self.paras["only_simple"]:
                microsatellite.call_init()
                microsatellite.call_micro()
            else:
                microsatellite.call_init()
                microsatellite.call_micro_and_other()

            if microsatellite.support_reads_totalNum == 0:
                continue
            if ((microsatellite.support_reads_mutNum[0] / microsatellite.support_reads_totalNum) < self.minimum_support_reads_ratio) and\
                ((microsatellite.support_reads_mutNum[1] / microsatellite.support_reads_totalNum) < self.minimum_support_reads_ratio):
                continue

            microsatellites.append(microsatellite)

        self.microsatellites = {ms.ms_id: ms for ms in microsatellites}

    def write_to_vcf_call_variants_complex(self, file_output):
        recs = []
        for ms_id in self.microsatellites_id:
            ms = self.microsatellites[ms_id]
            print(get_value("case"))
            print(ms.format_GT)
            print("ALT", ms.alt_str, ms.alt)
            print(ms.ref_str)

            vcfrec = file_output.new_record()
            vcfrec.contig = ms.chrom
            vcfrec.stop = ms.start + ms.repeat_times * ms.repeat_unit_len
            vcfrec.pos = ms.start
            vcfrec.ref = ms.ref_str
            vcfrec.alts = ms.alt
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            vcfrec.info["Quality"] = "|".join(map(str, [ms.qual_ms, ms.qual_ms_hap1, ms.qual_ms_hap2]))

            vcfrec.samples[get_value("case")]["GT"] = ms.format_GT
            vcfrec.samples[get_value("case")]["DP"] = ms.format_DP
            vcfrec.samples[get_value("case")]["QL"] = ms.format_QL
            vcfrec.samples[get_value("case")]["AL"] = ms.format_AL
            vcfrec.samples[get_value("case")].phased = ms.reads_phased

            recs.append(vcfrec)
        return recs

    def write_to_vcf_call_variants_micro(self, file_output):
        recs = []
        for ms_id in self.microsatellites_id:

            if ms_id not in self.microsatellites.keys():continue

            ms = self.microsatellites[ms_id]

            vcfrec = file_output.new_record()
            vcfrec.contig = ms.chrom
            vcfrec.stop = ms.start + ms.repeat_times * ms.repeat_unit_len
            vcfrec.pos = ms.start
            vcfrec.ref=str(ms.repeat_times)+"["+ms.repeat_unit+"]"
            if ms.report_micro:
                vcfrec.alts = ms.alt_ms
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = int(ms.end)
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = int(ms.repeat_len)
            vcfrec.info["query_repeat_length"] = int(ms.query_repeat_length)
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            vcfrec.info["Quality"] = "|".join(map(str, [ms.qual_ms, ms.qual_ms_hap1, ms.qual_ms_hap2]))
            if ms.report_micro:
                vcfrec.samples[get_value("case")]["GT"] = ms.format_GT_ms
                vcfrec.samples[get_value("case")]["DP"] = ms.format_DP_ms
                vcfrec.samples[get_value("case")]["QL"] = ms.format_QL_ms
                vcfrec.samples[get_value("case")]["AL"] = ms.format_AL_ms
                vcfrec.samples[get_value("case")].phased = ms.reads_phased
            recs.append(vcfrec)
        return recs

    def write_to_vcf_call_variants_indel(self, file_output):
        recs = []
        for ms_id in self.microsatellites_id:
            ms = self.microsatellites[ms_id]
            print(get_value("case"))
            print(ms.format_GT)
            print("ALT", ms.alt_str, ms.alt)
            print(ms.ref_str)
            vcfrec = file_output.new_record()
            vcfrec.contig = ms.chrom
            vcfrec.stop = ms.start + ms.repeat_times * ms.repeat_unit_len
            vcfrec.pos = ms.start
            vcfrec.ref = ms.ref_str
            vcfrec.alts = ms.alt
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            vcfrec.info["Quality"] = "|".join(map(str, [ms.qual_ms, ms.qual_ms_hap1, ms.qual_ms_hap2]))

            vcfrec.samples[get_value("case")]["GT"] = ms.format_GT
            vcfrec.samples[get_value("case")]["DP"] = ms.format_DP
            vcfrec.samples[get_value("case")]["QL"] = ms.format_QL
            vcfrec.samples[get_value("case")]["AL"] = ms.format_AL
            vcfrec.samples[get_value("case")].phased = ms.reads_phased

            recs.append(vcfrec)
        return recs

    def write_to_vcf_call_variants_snv(self, file_output):
        recs = []
        for ms_id in self.microsatellites_id:
            ms = self.microsatellites[ms_id]
            print(get_value("case"))
            print(ms.format_GT)
            print("ALT", ms.alt_str, ms.alt)
            print(ms.ref_str)

            vcfrec = file_output.new_record()
            vcfrec.contig = ms.chrom
            vcfrec.stop = ms.start + ms.repeat_times * ms.repeat_unit_len
            vcfrec.pos = ms.start
            vcfrec.ref = ms.ref_str
            vcfrec.alts = ms.alt
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            vcfrec.info["Quality"] = "|".join(map(str, [ms.qual_ms, ms.qual_ms_hap1, ms.qual_ms_hap2]))

            vcfrec.samples[get_value("case")]["GT"] = ms.format_GT
            vcfrec.samples[get_value("case")]["DP"] = ms.format_DP
            vcfrec.samples[get_value("case")]["QL"] = ms.format_QL
            vcfrec.samples[get_value("case")]["AL"] = ms.format_AL
            vcfrec.samples[get_value("case")].phased = ms.reads_phased

            recs.append(vcfrec)
        return recs

    def write_to_vcf_call_variants(self, file_output):
        recs = []
        for ms_id in self.microsatellites_id:
            ms = self.microsatellites[ms_id]
            print(get_value("case"))
            print(ms.format_GT)
            print("ALT", ms.alt_str, ms.alt)
            print(ms.ref_str)

            vcfrec = file_output.new_record()
            vcfrec.contig = ms.chrom
            vcfrec.stop = ms.start + ms.repeat_times * ms.repeat_unit_len
            vcfrec.pos = ms.start
            vcfrec.ref = ms.ref_str
            vcfrec.alts = ms.alt
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            vcfrec.info["Quality"] = "|".join(map(str, [ms.qual_ms, ms.qual_ms_hap1, ms.qual_ms_hap2]))

            vcfrec.samples[get_value("case")]["GT"] = ms.format_GT
            vcfrec.samples[get_value("case")]["DP"] = ms.format_DP
            vcfrec.samples[get_value("case")]["QL"] = ms.format_QL
            vcfrec.samples[get_value("case")]["AL"] = ms.format_AL
            vcfrec.samples[get_value("case")].phased = ms.reads_phased
            recs.append(vcfrec)
        return recs

    def write_to_vcf_pre_stat(self, file_output):
        recs = []
        for ms_id in self.microsatellites_id:
            ms = self.microsatellites[ms_id]
            vcfrec = file_output.new_record()
            vcfrec.contig = ms.chrom
            vcfrec.pos = ms.start
            vcfrec.ref = "."
            vcfrec.alts = (ms.alt_str,) if ms.alt_str != "" else (".",)
            vcfrec.id = ms.ms_id
            vcfrec.stop = ms.end
            vcfrec.info["ms_start"] = ms.start
            vcfrec.info["ms_end"] = ms.end
            vcfrec.info["motif"] = ms.repeat_unit
            vcfrec.info["repeat_times"] = ms.repeat_times
            vcfrec.info["motif_len"] = ms.repeat_unit_len
            vcfrec.info["ref_repeat_length"] = ms.repeat_len
            vcfrec.info["query_repeat_length"] = ms.query_repeat_length
            vcfrec.info["depth"] = ms.depth
            vcfrec.info["dis_stat"] = str(ms.dis_stat)
            vcfrec.info["dis"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis.items()]
            )
            vcfrec.info["dis_hap0"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap0.items()]
            )
            vcfrec.info["dis_hap1"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap1.items()]
            )
            vcfrec.info["dis_hap2"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_hap2.items()]
            )
            vcfrec.info["dis_forward"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_forward.items()]
            )
            vcfrec.info["dis_reversed"] = "|".join(
                [str(key) + ":" + str(value) for key, value in ms.ms_dis_reversed.items()]
            )
            recs.append(vcfrec)
        return recs

    def run_window_pre_stat(self):
        """
        For ngs/ccs pre_stat
        Returns:
        """
        self.init_microsatellites()
        self.init_reads()
        self.get_reads_dis()
        self.merge_reads_repeat_length_distribution()

    def run_window_call_variant(self):
        """
        For ngs/ccs variant calling
        Returns:

        """
        self.init_microsatellites(only_simple=self.paras["only_simple"])
        self.init_reads()
        self.get_reads_info_improveHMM_2(only_simple=self.paras["only_simple"])
        self.merge_muts_info(only_simple=self.paras["only_simple"])
        self.call_variants()
