#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import numpy as np
from sklearn import mixture
from collections import Counter
from sklearn.cluster import KMeans

def get_more_times(repeat_list):

    counts = Counter(repeat_list)
    f = sorted(zip(counts.values(), counts.keys()))
    if len(f) == 1:
        return f[-1][-1]
    else:
        if abs(f[-1][0] - f[-2][0]) / len(repeat_list) < 0.3:
            return int(round(np.mean(repeat_list)))
        else:
            return f[-1][-1]

def get_repeat_gmm(dis, target=1,k_para=3):

    support = sum(dis.values())
    cluster = len(dis)
    if cluster == 1:
        genotype = list(dis.keys()) if target == 1 else list(dis.keys()) * 2
        qual = support
        return {"genotype": genotype, "qual": qual}

    repeat_list = []
    for k, v in dis.items():
        repeat_list.extend([k] * v)
    repeat_list = np.array(repeat_list).reshape(-1, 1)

    dpgmm = mixture.GaussianMixture(n_components=cluster,
                                            covariance_type='full',
                                            tol=0.0001,
                                            max_iter=400)
    dpgmm.fit(repeat_list)
    pre_dis = {}
    pre_num = {}
    for k, v in dis.items():
        k_pre = dpgmm.predict(np.array([[k]]))
        if k_pre[0] not in pre_dis:
            pre_dis[k_pre[0]] = []
            pre_num[k_pre[0]] = 0
        pre_dis[k_pre[0]].extend([k] * v)
        pre_num[k_pre[0]] += v
    m = sorted(pre_num.keys(), key=(lambda x: pre_num[x]))

    if target == 1:
        if len(m) > 1:
            genotype = [get_more_times(pre_dis[m[-1]] + pre_dis[m[-2]])]
            qual = (pre_num[m[-1]] + pre_num[m[-2]]) / (1 + np.std(pre_dis[m[-1]]))
        else:
            genotype = [get_more_times(pre_dis[m[-1]])]
            qual = (pre_num[m[-1]]) / (1 + np.std(pre_dis[m[-1]]))

    elif target == 2:
        support_para = 0.8
        if pre_num[m[-1]] > support * support_para:

            genotype = [get_more_times(pre_dis[m[-1]])] * 2
            qual = support / (1 + np.std(pre_dis[m[-1]]))

        else:
            v_list = list(dis.values())
            v_list.sort()

            hap1 = get_more_times(pre_dis[m[-1]])
            qual1 = pre_num[m[-1]] / (1 + np.std(pre_dis[m[-1]]))
            hap2 = get_more_times(pre_dis[m[-2]])
            qual2 = pre_num[m[-2]] / (1 + np.std(pre_dis[m[-2]]))
            genotype = [hap1, hap2]
            qual = (qual1 + qual2) / 2

    return {"genotype": genotype, "qual": qual}


