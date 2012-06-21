import os

csmm_config = {
        "CSMM_INSTALL_DIR":"./",
        }

def get_config(key):
    value = os.getenv(key)
    if not value:
        value = csmm_config[key]
    return value

def get_datafile(filename):
    datadir = os.path.join(get_config("CSMM_INSTALL_DIR"), "data")
    return os.path.join(datadir, os.path.basename(filename))
