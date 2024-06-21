#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from units import *
from pre_stat import pre_stat
from estimate_error import estimate_error
from call_variants import call_variants
from global_dict import *

import yaml


def load_model(paras):
    model_path = paras["output_model"]
    stream = open(model_path, "r")
    model = yaml.load(stream.read(), Loader=yaml.FullLoader)
    set_value("model", model)


def genotype_ccs(paras):
    df_microsatellites = load_microsatellites(paras)
    call_variants(df_microsatellites)
