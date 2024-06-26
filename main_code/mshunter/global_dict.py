#!/usr/bin/env python
# -*- coding: UTF-8 -*-
import logging

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
consoleHandler = logging.StreamHandler()
consoleHandler.setLevel(logging.DEBUG)
formatter = logging.Formatter('[%(levelname)s]\t %(message)s')
consoleHandler.setFormatter(formatter)
logger.addHandler(consoleHandler)
version_id = "0.1.5"

def global_init():
    global _global_dict
    _global_dict = {}
    _global_dict["ms_number"] = 0
    _global_dict["tools_version"] = version_id
    _global_dict["tools_name"] = "MSHunter"
    _global_dict["author"] = "Peng Jia et al."
    _global_dict["email"] = "pengjia@stu.xjtu.edu.cn"
    _global_dict["description"] = "Microsatellite "
    _global_dict["chrom_list"] = [str(i) for i in range(1, 23)] + \
                                 ["chr" + str(i) for i in range(1, 23)] + \
                                 ["X", "Y", "chrX", "chrY", "chrM", "MT"]
    _global_dict["default"] = {
        "genotype": {
            "reference": ".",
            "threads": 4,
            "minimum_mapping_quality": 1,
            "minimum_support_reads": 2,
            "batch": 2000,
            "debug": "False",
            "microsatellite_region_format": "msisensor_scan",
            "only_homopolymers": False,
            "only_simple": "False",
            "using_phasing_info": True,
            "allow_mismatch": True,
            "minimum_repeat_times": "1:8;2-5:5",
            "maximum_repeat_times": "1-5:100",
            "prefix_len": 5,
            "suffix_len": 5,
            "sequencing_error": 0.001,
            "maximum_distance_of_two_complex_events": 5,

            "minimum_phasing_reads": 3,
            "hap": False,
            "min_allele_fraction": 0.2,

            "minimum_support_reads_ratio": 0.2,

            'mut_information_visual':False,
        },

        "qc": {
            "reference": ".",
            "threads": 4,
            "minimum_mapping_quality": 1,
            "minimum_support_reads": 2,
            "batch": 2000,
            "debug": "False",
            "microsatellite_region_format": "msisensor_scan",
            "only_homopolymers": "False",
            "allow_mismatch": "True",
            "minimum_repeat_times": "1:8;2-5:5",
            "maximum_repeat_times": "1-5:100",
            "prefix_len": 20,
            "suffix_len": 20,
            "minimum_phasing_reads": 3,
            "sequencing_error": 0.001,
            "hap": False,
            "min_allele_fraction": 0.2,
        },

    }


def set_value(name, value):
    _global_dict[name] = value


def get_value(name, defValue=None):
    try:
        return _global_dict[name]
    except KeyError:
        print("[ERROR] No variable", name, "in global_dict")
        return defValue
