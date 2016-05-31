from gdsctools import *
from nose.plugins.attrib import attr

@attr("onweb")
def test_help():
    gdsctools_help()
    gdsctools_help(IC50)
    
