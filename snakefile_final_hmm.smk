import os.path
import os
import copy
import numpy
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import get_clearSite
import get_clearSite_pj
from remakeData import *
from accuracy import *
from draw_all import *
from Site_screen import *
from evaluating_indicator import *

NUM_PBSIM = ['%04d' % i for i in range(1, 2,1)]

MUTlEN = [10]

PBSIM_DEPTH = 30
PBSIM_ACCURACY = 1
PBSIM_LEN = 15000
PBSIM_LEN_SD=1800
sreenSite_STRLen=[0,10000]

motifLen = 'all'
# motifLen = [1,2,3,4,5,6,'6+']
colors = ['gold','orange','y','r','g','c','b','m','purple','tan','sage','lime','teal','cyan','navy','black']

wildcard_constraints:
    mut_len="|".join([f'{i}' for i in MUTlEN]),
    num = "|".join([f'{i}' for i in NUM_PBSIM]),

#******************hg38/chr1:前0-11000000bp
# #输入：hg38/chr1:前0-11000000bp
# COMPAR_SOFTWARE ='minimap2' #用户输入的比对软件的名称——'minimap2'/'winnowmap2'/'ngmlr'/'lra'
# FA_INITIAL = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/chr1_0_11000000.fa' #改造前的、原始的FA文件
# VNTR_BED = '/mnt/d/science/simulateData/remakeData/clearSite/GRCh38_VNTR_chr1.sort.bed' #用来筛选位点的VNTR区域
# SR_bed = '/mnt/d/science/simulateData/remakeData/clearSite/simplerepeat_chr1.bed'
# SD_BED = '/mnt/d/science/simulateData/remakeData/clearSite/seg_dup_chr1.bed'
# FILE_ComplexBed = [VNTR_BED,SD_BED,SR_bed]
# # PBSIM1_MODEL = "/mnt/d/pbsim/data/model_qc_ccs" #pbsim规则的model——版本1
# PBSIM3_MODEL = '/mnt/d/science/simulateData/remakeData/STRtool/pbsim3/data/ERRHMM-SEQUEL.model' #pbsim规则的model——版本3
# MSHunter_PYTHON = '/home/zyx/miniconda3/envs/snakemake_env/bin/python3.10' #MSHunter的python解释器
# MSHunter_MAIN = '/mnt/d/science/simulateData/remakeData/STRtool/MSHunter/mshunter/main.py' #MSHunter的main文件
#
# #输出:hg38/chr1:前0-11000000bp
# IMAGE_PATH = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/heatmap_STRLEN.jpg' #最终输出的热图路径
# INDICATOR_IMAGE_PATH = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/precision-recall_{mut_len}.jpg' #precision-recall曲线
# IMAGE_ERROR_PATH = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/error_STRLEN.jpg' #error曲线图
# IMAGE_ERROR_PATH_ALL = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/error_all_STRLEN.jpg' #error曲线图
# IMAGE_Histogram_PATH_REAL = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/Histogram_real_{mut_len}.jpg'
# IMAGE_Histogram_PATH_PREDICT = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/Histogram_predict_{mut_len}.jpg'
# IMAGE_ERROR_MEAN_PATH = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/error_mean_STRLEN.jpg' #不同变异长度的error均值曲线
# IMAGE_PATH_REAL_ALL = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/Histogram_real_all.jpg' #柱状图：不同变异长度总的STR位点数量
# IMAGE_PATH_PREDICT_ALL = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/Histogram_predict_all.jpg'
# DF_EXCEL_ERROR = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/excel_error_all_STRLEN.xlsx' #用来画不同motif len的mean error的表格
# DF_EXCEL_HISTOGRAM = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/excel_histogram_all_{motifLen}.xlsx' #用来画不同motif len的str num的表格
# LIST_1 = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/hg38Chr1_list.list' #getList规则：msi-pro得到的list文件
# LIST_2 = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/hg38Chr1_clear.list' #clearSite规则：对VNTR/SD/SR筛选后的list文件
# LIST_3 = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/hg38Chr1_clear_sreen.list' #sreenSite规则：根据位点长度，减少位点数后的list文件
# LIST_4 = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/hg38Chr1_clear_sreen_motifLen.list' #sreenSite_motifLen规则
# REMAKED_FA = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/chr1_after.{mut_len}.fa' #改造后的fa文件(纯合)
# REMAKED_TXT = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/remake_chr1.{mut_len}.txt' #具体改造情况(纯合)
# PBSIM_FQ_PART = '/mnt/d/science/simulateData/remakeData/STRtool/sd_{mut_len}_{num}.fastq' #纯合/合并后的杂合路径都是这个——expand
# PBSIM_REF_PART = '/mnt/d/science/simulateData/remakeData/STRtool/sd_{mut_len}_{num}.ref' #expand
# PBSIM_REF_FAI = '/mnt/d/science/simulateData/remakeData/STRtool/sd_{mut_len}_{num}.ref.fai'
# #内存占用太多，换一个地方
# PBSIM_FQ = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/sd_{mut_len}_all.fastq' #纯合/合并后的杂合路径都是这个——expand
# PBSIM_REF = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/sd_{mut_len}_all.ref' #expand
# MINIMAP2_INDEX = "/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/mini2_ref.mmi" #minimap2_index规则
# MINIMAP2_SAM = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/mini2_{mut_len}.sam'
# MINIMAP2_BAM = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/mini2_{mut_len}.bam'
# #内存占用太多，换一个地方
# MINIMAP2_BAM_SORTED = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/mini2_{mut_len}.sorted.bam'
# MINIMAP2_BAM_SORTED_BAI = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/mini2_{mut_len}.sorted.bam.bai'
# VCF_GZ = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/output_test_{mut_len}/mini2_{mut_len}_micro.vcf.gz'
# VCF = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3/output_test_{mut_len}/mshunter_micro.vcf'

#******************hg38/chr1
#输入：hg38/chr1
COMPAR_SOFTWARE ='lra' #用户输入的比对软件的名称——'minimap2'/'winnowmap2'/'ngmlr'/'lra'
# FA_INITIAL = '/mnt/d/science/simulateData/remakeData/realDataToMshunter/chr1.fa' #改造前的、原始的FA文件
FA_INITIAL = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/chr1_1000000.fa' #改造前的、原始的FA文件
VNTR_BED = '/mnt/d/science/simulateData/remakeData/clearSite/GRCh38_VNTR_chr1.sort.bed' #用来筛选位点的VNTR区域
SR_bed = '/mnt/d/science/simulateData/remakeData/clearSite/simplerepeat_chr1.bed'
SD_BED = '/mnt/d/science/simulateData/remakeData/clearSite/seg_dup_chr1.bed'
FILE_ComplexBed = [VNTR_BED,SD_BED,SR_bed]
# PBSIM1_MODEL = "/mnt/d/pbsim/data/model_qc_ccs" #pbsim规则的model——版本1
PBSIM3_MODEL = '/mnt/d/science/simulateData/remakeData/STRtool/pbsim3/data/ERRHMM-SEQUEL.model' #pbsim规则的model——版本3
# PBSIM3_MODEL = '/mnt/d/science/simulateData/remakeData/STRtool/pbsim3/data/ERRHMM-RSII.model' #pbsim规则的model——版本3
# PBSIM3_MODEL = '/mnt/d/science/simulateData/remakeData/STRtool/pbsim3/data/QSHMM-RSII.model' #pbsim规则的model——版本3
MSHunter_PYTHON = '/home/zyx/miniconda3/envs/snakemake_env/bin/python3.10' #MSHunter的python解释器
MSHunter_MAIN = '/mnt/d/science/simulateData/remakeData/STRtool/MSHunter/mshunter/main.py' #MSHunter的main文件

#输出:hg38/chr1
IMAGE_PATH = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/heatmap_{mut_len}.jpg.jpg' #最终输出的热图路径
INDICATOR_IMAGE_PATH = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/precision-recall_{mut_len}.jpg' #precision-recall曲线
IMAGE_ERROR_PATH = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/error_{mut_len}.jpg' #error曲线图
IMAGE_ERROR_PATH_ALL = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/error_all_STRLEN.jpg' #error曲线图
IMAGE_Histogram_PATH_REAL = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/Histogram_real_{mut_len}.jpg'
IMAGE_Histogram_PATH_PREDICT = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/Histogram_predict_{mut_len}.jpg'
IMAGE_ERROR_MEAN_PATH = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/error_mean_STRLEN.jpg' #不同变异长度的error均值曲线
IMAGE_PATH_REAL_ALL = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/Histogram_real_all.jpg' #柱状图：不同变异长度总的STR位点数量
IMAGE_PATH_PREDICT_ALL = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/Histogram_predict_all.jpg'
DF_EXCEL_ERROR = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/excel_error_all_STRLEN.xlsx' #用来画不同motif len的mean error的表格
DF_EXCEL_HISTOGRAM = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/excel_histogram_all_{motifLen}.xlsx' #用来画不同motif len的str num的表格
# LIST_1 = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/hg38Chr1_list.list' #getList规则：msi-pro得到的list文件 —— 之前一直用的这个，这个没有motifLen=6
LIST_1 = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/hg38Chr1_list_have6.list' #getList规则：msi-pro得到的list文件 —— 之前一直用的这个，这个没有motifLen=6
LIST_2 = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/hg38Chr1_clear.list' #clearSite规则：对VNTR/SD/SR筛选后的list文件
LIST_3 = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/hg38Chr1_clear_sreen.list' #sreenSite规则：根据位点长度，减少位点数后的list文件
LIST_4 = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/hg38Chr1_clear_sreen_motifLen.list' #sreenSite_motifLen规则
REMAKED_FA = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/chr1_after.{mut_len}.fa' #改造后的fa文件(纯合)
REMAKED_TXT = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/remake_chr1.{mut_len}.txt' #具体改造情况(纯合)

#原来的
PBSIM_FQ_PART = '/mnt/d/science/simulateData/remakeData/STRtool/sd_{mut_len}_{num}.fastq' #纯合/合并后的杂合路径都是这个——expand
PBSIM_REF_PART = '/mnt/d/science/simulateData/remakeData/STRtool/sd_{mut_len}_{num}.ref' #expand
PBSIM_REF_FAI = '/mnt/d/science/simulateData/remakeData/STRtool/sd_{mut_len}_{num}.ref.fai'
#内存占用太多，换一个地方
PBSIM_FQ = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/sd_{mut_len}_all.fastq' #纯合/合并后的杂合路径都是这个——expand
PBSIM_REF = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/sd_{mut_len}_all.ref' #expand
MINIMAP2_INDEX = "/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/mini2_ref.mmi" #minimap2_index规则
MINIMAP2_SAM = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/mini2_{mut_len}.sam'
MINIMAP2_BAM = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/mini2_{mut_len}.bam'
#内存占用太多，换一个地方
MINIMAP2_BAM_SORTED = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/mini2_{mut_len}.sorted.bam'
MINIMAP2_BAM_SORTED_BAI = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/mini2_{mut_len}.sorted.bam.bai'
VCF_GZ = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/output_test_{mut_len}/mini2_{mut_len}_micro.vcf.gz'
VCF = '/mnt/d/science/simulateData/final/补充实验/不同比对软件/LRA/output_test_{mut_len}/mshunter_micro.vcf'

#****全基因组的输出
# IMAGE_PATH = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/heatmap_{mut_len}.jpg.jpg' #最终输出的热图路径
# INDICATOR_IMAGE_PATH = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/precision-recall_{mut_len}.jpg' #precision-recall曲线
# IMAGE_ERROR_PATH = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/error_{mut_len}.jpg' #error曲线图
# IMAGE_ERROR_PATH_ALL = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/error_all_STRLEN.jpg' #error曲线图
# IMAGE_Histogram_PATH_REAL = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/Histogram_real_{mut_len}.jpg'
# IMAGE_Histogram_PATH_PREDICT = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/Histogram_predict_{mut_len}.jpg'
# IMAGE_ERROR_MEAN_PATH = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/error_mean_STRLEN.jpg' #不同变异长度的error均值曲线
# IMAGE_PATH_REAL_ALL = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/Histogram_real_all.jpg' #柱状图：不同变异长度总的STR位点数量
# IMAGE_PATH_PREDICT_ALL = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/Histogram_predict_all.jpg'
# DF_EXCEL_ERROR = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/excel_error_all_STRLEN.xlsx' #用来画不同motif len的mean error的表格
# DF_EXCEL_HISTOGRAM = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/excel_histogram_all_{motifLen}.xlsx' #用来画不同motif len的str num的表格
# # LIST_1 = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/hg38Chr1_list.list' #getList规则：msi-pro得到的list文件
# # LIST_2 = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/hg38Chr1_clear.list' #clearSite规则：对VNTR/SD/SR筛选后的list文件
# # LIST_3 = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/hg38Chr1_clear_sreen.list' #sreenSite规则：根据位点长度，减少位点数后的list文件
# # LIST_4 = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/hg38Chr1_clear_sreen_motifLen.list' #sreenSite_motifLen规则
# # REMAKED_FA = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/chr1_after.{mut_len}.fa' #改造后的fa文件(纯合)
# REMAKED_TXT = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/remake.{mut_len}.txt' #具体改造情况(纯合)
#
# #原来的
# # PBSIM_FQ_PART = '/mnt/d/science/simulateData/remakeData/STRtool/sd_{mut_len}_{num}.fastq' #纯合/合并后的杂合路径都是这个——expand
# # PBSIM_REF_PART = '/mnt/d/science/simulateData/remakeData/STRtool/sd_{mut_len}_{num}.ref' #expand
# # PBSIM_REF_FAI = '/mnt/d/science/simulateData/remakeData/STRtool/sd_{mut_len}_{num}.ref.fai'
# # #内存占用太多，换一个地方
# # PBSIM_FQ = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/sd_{mut_len}_all.fastq' #纯合/合并后的杂合路径都是这个——expand
# # PBSIM_REF = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/sd_{mut_len}_all.ref' #expand
# # MINIMAP2_INDEX = "/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/mini2_ref.mmi" #minimap2_index规则
# # MINIMAP2_SAM = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/mini2_{mut_len}.sam'
# # MINIMAP2_BAM = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/mini2_{mut_len}.bam'
# # #内存占用太多，换一个地方
# # MINIMAP2_BAM_SORTED = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/mini2_{mut_len}.sorted.bam'
# # MINIMAP2_BAM_SORTED_BAI = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/mini2_{mut_len}.sorted.bam.bai'
# # VCF_GZ = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/output_test_{mut_len}/mini2_{mut_len}_micro.vcf.gz'
# VCF = '/mnt/d/science/simulateData/final/hmm/HMM_pbsim3_2/全基因组/ms_res_micro.vcf'

rule all:
    input:
        image_path = expand(IMAGE_PATH,mut_len=MUTlEN),
        image_error_path= expand(IMAGE_ERROR_PATH,mut_len=MUTlEN),
        image_error_path_all= IMAGE_ERROR_PATH_ALL,
        # image_path_real = expand(IMAGE_Histogram_PATH_REAL,mut_len=MUTlEN),
        # image_path_predict= expand(IMAGE_Histogram_PATH_PREDICT,mut_len=MUTlEN),
        image_errorMean_path=IMAGE_ERROR_MEAN_PATH,
        # image_path_real_all= IMAGE_PATH_REAL_ALL,
        # image_path_predict_all= IMAGE_PATH_PREDICT_ALL,
        df_excel_error= expand(DF_EXCEL_ERROR,motifLen=motifLen),
        # df_excel_histogram= expand(DF_EXCEL_HISTOGRAM,motifLen=motifLen)

        image_path_real= expand(IMAGE_Histogram_PATH_REAL,mut_len=MUTlEN),
        image_path_predict= expand(IMAGE_Histogram_PATH_PREDICT,mut_len=MUTlEN),
        # image_path_real_all= IMAGE_PATH_REAL_ALL,
        # image_path_predict_all= IMAGE_PATH_PREDICT_ALL,
        df_excel_histogram= expand(DF_EXCEL_HISTOGRAM,motifLen=motifLen)


rule getList:
    input:
        fa = FA_INITIAL
    output:
        list = LIST_1
    shell:
        'msisensor-pro scan -d {input.fa} -o {output.list} -s 6 -m 1000'

rule clearSite:
    input:
        file_MSIbed=LIST_1,
        file_VNTRbed=VNTR_BED,
        file_SDbed=SD_BED,
        file_SRbed=SR_bed,
        file_complexBed = FILE_ComplexBed
    output:
        list = LIST_2
    run:
        get_clearSite_pj.getClearSite(input.file_MSIbed,input.file_complexBed,output.list)

rule sreenSite:
    input:
        file_bed=LIST_2
    output:
        file_bed_out = LIST_3
    run:
        site_screen_len(input.file_bed,output.file_bed_out,STRLen = sreenSite_STRLen)

rule sreenSite_motifLen:
    input:
        file_bed=LIST_3
    output:
        file_bed_out = LIST_4
    run:
        site_screen_motifLen(input.file_bed,output.file_bed_out,motifLen)

rule remakeBigData:
    input:
        file_bed = LIST_4,
        file_fa=FA_INITIAL
    output:
        file_out_fa = expand(REMAKED_FA,mut_len=MUTlEN),
        file_out_remake = expand(REMAKED_TXT,mut_len=MUTlEN),
    run:
        for i in range(len(MUTlEN)):
            this_mut_len = MUTlEN[i]
            file_out_fa_temp = output.file_out_fa[i]
            file_out_remake_temp = output.file_out_remake[i]

            str_remake_gap=100

            if this_mut_len<0:
                print('ERROR:The variaLen can not be negative! ',this_mut_len)
                remakeData_len_homo_mulChr(f"{input.file_fa}",f"{input.file_bed}",this_mut_len,f"{file_out_fa_temp}",f"{file_out_remake_temp}",method_num=1,str_remake_gap=str_remake_gap)
            elif this_mut_len>0:
                print('The variaLen being remaked is ',this_mut_len)
                remakeData_len_homo_mulChr(f"{input.file_fa}",f"{input.file_bed}",this_mut_len,f"{file_out_fa_temp}",f"{file_out_remake_temp}",method_num=5,str_remake_gap=str_remake_gap)
            else:
                print('This mutation length is zero!')
                remakeData_len_homo_mulChr(f"{input.file_fa}",f"{input.file_bed}",this_mut_len,f"{file_out_fa_temp}",f"{file_out_remake_temp}",method_num=5,str_remake_gap=str_remake_gap)


rule pbsim:
    input:
        model=PBSIM3_MODEL,
        file_out_remakeFa=expand(REMAKED_FA,mut_len=MUTlEN)
    output:
        pbsim_fq=expand(PBSIM_FQ_PART,mut_len=MUTlEN,num=NUM_PBSIM),
        pbsim_ref=expand(PBSIM_REF_PART,mut_len=MUTlEN,num=NUM_PBSIM),
    run:
        for i in range(len(MUTlEN)):
            this_mut_len = MUTlEN[i]
            input_file_out_remakeFa = input.file_out_remakeFa[i]

            shell('pbsim --strategy wgs --method errhmm --errhmm {input.model} --length-mean {PBSIM_LEN}   --depth {PBSIM_DEPTH}  --prefix sd_{this_mut_len}  --genome {input_file_out_remakeFa}'
                  ' --accuracy-mean {PBSIM_ACCURACY} --length-sd {PBSIM_LEN_SD}  ' )

rule merge_pbsim_fq:
    input:
        pbsim_fq = expand(PBSIM_FQ_PART,mut_len=MUTlEN,num=NUM_PBSIM),
    output:
        pbsim_fq_all=expand(PBSIM_FQ,mut_len=MUTlEN),
    run:
        for i in range(len(MUTlEN)):
            this_mut_len = MUTlEN[i]
            input_fq = input.pbsim_fq[i]
            output_fq = output.pbsim_fq_all[i]
            shell('cat {input_fq} > {output_fq}')


rule merge_pbsim_ref:
    input:
        pbsim_ref = expand(PBSIM_REF_PART,mut_len=MUTlEN,num=NUM_PBSIM)
    output:
        pbsim_ref_all=expand(PBSIM_REF,mut_len=MUTlEN)
    run:
        for i in range(len(MUTlEN)):
            this_mut_len = MUTlEN[i]
            input_fq = input.pbsim_ref[i]
            output_fq = output.pbsim_ref_all[i]
            shell('cat {input_fq} > {output_fq}')

rule minimap2_index:
    input:
        fa = FA_INITIAL
    output:
        index = MINIMAP2_INDEX
    shell:
        "minimap2 -H -d {output.index}  {input.fa}"

rule align:
    input:
        fa=FA_INITIAL,
        pbsim_fq = PBSIM_FQ
    output:
        sam = MINIMAP2_SAM
    run:
        if COMPAR_SOFTWARE == 'minimap2':

            shell("minimap2  -ax map-hifi  --eqx " 
            "{input.fa} {input.pbsim_fq} > {output.sam} -t 8")

        elif COMPAR_SOFTWARE == 'winnowmap2':
            print('ERROR:winnowmap2不行！')
        elif COMPAR_SOFTWARE == 'ngmlr':
            shell("ngmlr -t 4 -r {input.fa}  -q {input.pbsim_fq} -o {output.sam} -x pacbio")
        elif COMPAR_SOFTWARE == 'lra':
            shell('lra index -CCS {input.fa}')
            shell('lra align -CCS {input.fa} {input.pbsim_fq} -t 8 -p s > {output.sam}')

        else:
            print('No such compare software!')

rule samToBam:
    input:
        sam = MINIMAP2_SAM
    output:
        bam=MINIMAP2_BAM
    run:
        if COMPAR_SOFTWARE == 'minimap2' or COMPAR_SOFTWARE == 'winnowmap2' or COMPAR_SOFTWARE == 'lra':
            shell('samtools view -q 20 -bS  {input.sam}  -o {output.bam}')
        elif COMPAR_SOFTWARE == 'ngmlr':
            shell("""
            grep "^@" {input.sam} >{output}
            grep -v "^@"  {input.sam} |awk '{{if ($5>=20) print $0 }}' >>{output}
            """)
        else:
            print('No such compare software!')

rule bamAddSort:
    input:
        bam=MINIMAP2_BAM
    output:
        bam_sorted=MINIMAP2_BAM_SORTED
    run:
        shell('samtools sort {input.bam} -o {output.bam_sorted}')

rule bamIndex:
    input:
        bam_sorted=MINIMAP2_BAM_SORTED
    output:
        bam_sorted_bai = MINIMAP2_BAM_SORTED_BAI
    run:
        shell('samtools index {input}')

rule refFaidx:
    input:
        pbsim_ref=PBSIM_REF
    output:
        pbsim_ref_fai=PBSIM_REF_FAI
    run:
        shell('samtools faidx {input.pbsim_ref}')

rule runMshunter:
    input:
        list=LIST_4,
        fa_ori = FA_INITIAL,

        bam=expand(MINIMAP2_BAM_SORTED,mut_len=MUTlEN),
        bai=expand(MINIMAP2_BAM_SORTED_BAI,mut_len=MUTlEN)
    output:
        vcfgz = expand(VCF_GZ,mut_len=MUTlEN),
        vcf = expand(VCF,mut_len=MUTlEN)
    run:
        for i in range(len(MUTlEN)):
            this_mut_len = MUTlEN[i]
            input_bam = input.bam[i]
            output_vcfgz = output.vcfgz[i]
            output_vcf = output.vcf[i]

            output_preifx = "/".join(output.vcfgz[i].split("/")[:-1])
            shell('rm -rf {output_preifx}')
            shell('{MSHunter_PYTHON} {MSHunter_MAIN} genotype -m {input.list} -i {input_bam} -r {input.fa_ori} -o {output_preifx} -tech ccs  '
                  '-cs {COMPAR_SOFTWARE} -t 8 -s mini2_{this_mut_len}')
            shell('gunzip -k -c {output_vcfgz} > {output_vcf}')

rule draw:
    input:
        file_predict = expand(VCF,mut_len=MUTlEN),
        file_true = expand(REMAKED_TXT,mut_len=MUTlEN)
    output:
        image_path = expand(IMAGE_PATH,mut_len=MUTlEN),
        image_error_path = expand(IMAGE_ERROR_PATH,mut_len=MUTlEN),
        image_error_path_all = IMAGE_ERROR_PATH_ALL,
        image_errorMean_path = IMAGE_ERROR_MEAN_PATH,
        df_excel_error= expand(DF_EXCEL_ERROR,motifLen=motifLen)
    run:
        error_list_all,R_list_all = [],[]
        mutlen_list,error_mean_list = [],[]

        for i in range(len(MUTlEN)):
            this_mut_len = MUTlEN[i]
            print('this variaLen is ',this_mut_len)

            file_true = input.file_true[i]
            file_predict = input.file_predict[i]
            output_image_path = output.image_path[i]

            draw_heatmap(file_true,file_predict,output_image_path,this_mut_len,motifLen,mode='homo')

            output_image_error_path = output.image_error_path[i]

            error_list,R_list = get_error(file_true,file_predict,output_image_error_path,this_mut_len,mode='homo',error_way='1/2RMSE')

            error_list_all.append(error_list)
            R_list_all.append(R_list)

            mutlen_list.append(this_mut_len)
            error_mean_list.append(sum(error_list)/len(error_list))

        plt.figure()
        plt.title('Error Curve',fontsize=15)
        plt.xlabel('Number of repeat units (real)')
        plt.ylabel('RMSE')
        for i in range(len(error_list_all)):
            plt.scatter(R_list_all[i],error_list_all[i],c=colors[i],s=10)
            plt.plot(R_list_all[i],error_list_all[i],label='variaLen is %s' % MUTlEN[i],c=colors[i])
        plt.legend(fontsize=8)
        plt.savefig(output.image_error_path_all)
        plt.close()

        plt.figure()
        plt.title('Mean Error Curve by Mutation Length',fontsize=15)
        plt.xlabel('Mutation Length(bp)')
        plt.ylabel('Mean Error')
        plt.scatter(mutlen_list,error_mean_list)
        plt.plot(mutlen_list,error_mean_list)
        plt.savefig(output.image_errorMean_path)
        plt.close()

        df_excel = pd.DataFrame(columns=['error','mean error'])
        df_excel['error'] = sum(error_list_all,[])
        error_mean_list_df = []
        for i in range(len(sum(error_list_all,[]))):
            error_mean_list_df.append(error_mean_list[0])
        df_excel['mean error'] = error_mean_list_df
        df_excel['repeat num'] = R_list_all[0]
        try:
            pd.concat(map(df_excel.to_excel,output.df_excel_error))
        except Exception as a:
            print(f'错误{a}')

rule draw_Histogram:
    input:
        file_predict = expand(VCF,mut_len=MUTlEN),
        file_true = expand(REMAKED_TXT,mut_len=MUTlEN)
    output:
        image_path_real = expand(IMAGE_Histogram_PATH_REAL,mut_len=MUTlEN),
        image_path_predict= expand(IMAGE_Histogram_PATH_PREDICT,mut_len=MUTlEN),
        image_path_real_all = IMAGE_PATH_REAL_ALL,
        image_path_predict_all= IMAGE_PATH_PREDICT_ALL,
        df_excel_histogram = expand(DF_EXCEL_HISTOGRAM,motifLen=motifLen)
    run:
        mutlen_list,str_num_real,str_num_predict = [],[],[]

        for i in range(len(MUTlEN)):
            this_mut_len = MUTlEN[i]
            print('(Histogram)this variaLen is ',this_mut_len)

            file_true = input.file_true[i]
            STR_num_all_real,STR_num_all_predict = draw_histogram(file_true,input.file_predict[i],output.image_path_real[i],output.image_path_predict[i],this_mut_len,motifLen,mode='homo') #同时画base='real'和'predict'
            mutlen_list.append(this_mut_len)
            str_num_real.append(STR_num_all_real)
            str_num_predict.append(STR_num_all_predict)

        plt.figure()
        plt.title('STR Num by Mutation Length(real)',fontsize=15)
        plt.xlabel('Mutation Length(bp)')
        plt.ylabel('STR Num')
        plt.bar(mutlen_list,str_num_real)
        plt.savefig(output.image_path_real_all)
        plt.close()

        plt.figure()
        plt.title('STR Num by Mutation Length(predict)',fontsize=15)
        plt.xlabel('Mutation Length(bp)')
        plt.ylabel('STR Num')
        plt.bar(mutlen_list,str_num_predict)
        plt.savefig(output.image_path_predict_all)
        plt.close()

        df_excel = pd.DataFrame(columns=['mutation length','str num(real)','str num(predict)'])
        df_excel['mutation length'] = mutlen_list
        df_excel['str num(real)'] = str_num_real
        df_excel['str num(predict)'] = str_num_predict
        try:
            pd.concat(map(df_excel.to_excel,output.df_excel_histogram))
        except Exception as a:
            print(f'错误{a}')

