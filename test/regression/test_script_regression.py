from gdsctools.scripts import regression
from gdsctools import gdsctools_data
import os
import pytest


prog = "gdsctools_regression"

#@pytest.fixture(scope='session')
def test_analysis():
    GF = gdsctools_data("genomic_features_v5.csv.gz")
    IC = gdsctools_data("IC50_v5.csv.gz")

    # Test that database must be provided
    import tempfile
    pname = tempfile.mkdtemp()
    df = regression.main([prog, '-F', GF, '-I', IC, "-O", pname, "--force"])


def test_help():
    regression.main([prog, '1>/tmp/out', '2>/tmp/err'])

def test_help2():
    regression.main([prog, '--help', '1>/tmp/out', '2>/tmp/err'])

def test_version():
    regression.main([prog, '--version', '1>/tmp/out', '2>/tmp/err'])


def test_license():
    regression.main([prog, '--license', '1>/tmp/out', '2>/tmp/err'])
