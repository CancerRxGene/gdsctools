from gdsctools.cosmictools import COSMICFetcher, COSMICInfo
from gdsctools.datasets import cosmic_builder_test
import pytest

import os
skiptravis = pytest.mark.skipif( "TRAVIS_PYTHON_VERSION" in os.environ,
     reason="On travis")
 


def test_cosmic_info():
    c = COSMICInfo()
    assert c.get(930301, 'SAMPLE_NAME') == 'VMRC-MELG'


def test_cosmic_fetcher():
    c = COSMICFetcher(cosmic_builder_test.filename)


def test_cosmic():
    r = COSMICInfo()
    r.get(924100)
    r.get('924100')
    r.get(924100, colname='SAMPLE_NAME')
    r.get(1000000000)


@skiptravis
def test_onweb():
    r = COSMICInfo()
    r.on_web(924100)
