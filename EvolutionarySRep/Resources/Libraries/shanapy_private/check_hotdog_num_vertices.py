import vtk
import joblib
import sys
sys.modules['sklearn.externals.joblib'] = joblib
from vtk.util import numpy_support
from shanapy.utils import viz
# from shanapy.models import srep_fitter as sf, mean_curvature_flow as mcf
from shanapy.models import mean_curvature_flow as mcf
from shanapy.models import srep_fitter as sf
import numpy as np
from matplotlib import pyplot as plt
from glob import glob
import pyvista as pv
import vtk
import os
import pyacvd
from glob import glob
from centers_of_curvature import get_spine_endpoints_from_mesh
from shanapy.models.sreps.refiner import Refiner
from shanapy.visualization import SrepViewer
from shanapy.models.ellipsoid_srep_coords import *
import pandas as pd

# hotdog_files = sorted(glob("/home/nick/dev/libigl/data/hotdog/deformed_hotdog??.vtk"))[::-1]
hotdog_files = sorted(glob("/home/nick/dev/libigl/data/hotdog/deformed_hotdog??.vtk"))

for f in hotdog_files:
    print(f, pv.PolyData(f).GetNumberOfPoints())
