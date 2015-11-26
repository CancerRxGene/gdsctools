# -*- python -*-
# -*- coding utf-8 -*-
#
#  This file is part of GDSCTools software
#
#  Copyright (c) 2015 - Wellcome Trust Sanger Institute
#  All rights reserved
#
#  File author(s): Thomas Cokelaer <cokelaer@gmail.comWE HERE>
#
#  Distributed under the BSD 3-Clause License.
#  See accompanying file LICENSE.txt distributed with this software
#
#  website: http://github.com/CancerRxGene/gdsctools
#
##############################################################################
"""Multicore implementation of the ANOVA code"""
import time
from gdsctools.anova import ANOVA


__all__ = ['multicore_anova']



def multicore_anova(ic50, genomic_features, maxcpu=2):
    """Using 4 cores, the entire analysis took 15 minutes using
    4 CPUs (16 Oct 2015).

    :param ic50: a filename or :class:`IC50` instance.
    :return: the anova instance itself (not the results); see example below.

    ::

        from gdsctools.anova import multicore
        master = multicore(dataset, maxcpu=2)
        results = master.anova_all()

        from gdsctools import ANOVAReport()
        report = ANOVAReport(master, results)
        report.create_html_pages(0

    .. warning:: experimental. Seems to work but sometimes hangs forever.
    """
    print("experimental code to run the analysis with several cores")
    print("May takes lots or resources and slow down your system")
    t1 = time.time()
    master = ANOVA(ic50, genomic_features=genomic_features, 
            low_memory=True)

    drugs = master.ic50.drugIds

    from easydev import MultiProcessing
    t = MultiProcessing(maxcpu=maxcpu)
    # add all jobs (one per drug)
    for i, drug in enumerate(drugs):
        t.add_job(analyse_one_drug, master, drug)
    t.run()

    # populate the ANOVA instance with the results
    for this in t.results:
        drug = this[0]
        result = this[1]
        master.individual_anova[drug] = result

    print("\nTook " + str(time.time() - t1) + "seconds.")
    return master


def analyse_one_drug(master, drug):
    res = master.anova_one_drug(drug_id=drug, animate=False)
    return (drug, res)


