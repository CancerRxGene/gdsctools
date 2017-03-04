from gdsctools import *

try:
    from unittest.mock import patch
except:
    from mock import patch

def test_help(mocker):
    def func(*args, **kwargs):
        pass
    with patch('webbrowser.open', func):
        with patch('webbrowser.open_new', func):
            gdsctools_help()
            gdsctools_help(IC50)

