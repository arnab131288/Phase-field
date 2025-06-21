"""
Manage output directory and write 2D fields to CSV.
"""
import os
import shutil
import numpy as np
from config import params

OUT = params['OUT_DIR']
# Remove existing output directory and recreate it
if os.path.exists(OUT):
    shutil.rmtree(OUT)
os.makedirs(OUT)

def write_csv(name, field):
    """Save a 2D numpy array as CSV in the output directory."""
    path = os.path.join(OUT, name)
    np.savetxt(path, field, delimiter=',')
