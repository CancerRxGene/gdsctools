import pkg_resources
try:
    version = pkg_resources.require("dreamtools")[0].version
    __version__ = version
except:
    # update this manually is possible when the version in the
    # setup changes
    version = "0.1"

