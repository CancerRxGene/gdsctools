import os

from gdsctools.report import ReportMAIN, HTMLTable
from gdsctools.tissues import TCGA
from gdsctools.readers import DrugDecode

import pandas as pd


class DataPackagesSummary(object):
    """Scan a set of data packages and create a summary report.


    In a given directory, identifies data packages and create a summary page

    """

    def __init__(self, collaborator=None, directories=None,
            directory='.'):
        pass
        # TODO  move GDSCDirectorySummary in here.
