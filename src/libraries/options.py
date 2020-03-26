# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 22:34:29 2020

Author: Rounak
"""

from recordtype import recordtype as rt
out = rt("Output", field_names=["All","sys_sum","area_sum","bus","branch","gen","lim",
                                "force", "suppress_detail"])
limit = rt("Limit", field_names=["All","v","line","pg","qg"])
options = rt("Options", field_names=["verbose","output"])


default_lim = limit(All=-1,v=1,line=1,pg=1,qg=1)
default_out = out(All=-1,sys_sum=1,area_sum=0,bus=1,branch=1,gen=0,
                         lim=default_lim, force=0, suppress_detail=-1)
default_opt = options(verbose=False,output=default_out)