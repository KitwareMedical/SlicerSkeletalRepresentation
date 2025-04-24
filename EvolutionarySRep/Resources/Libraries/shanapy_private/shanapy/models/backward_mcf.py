### Usage:
### Slicer --no-main-window --python-script backward_mcf.py
import vtk, qt, ctk, slicer
import subprocess
import os
import re
import sys
import argparse

output_root = '/playpen/workspace/my_paper/linking/data/nonaligned_hipp_sreps/'
iter_num = 500

slicer.modules.skeletalrepresentationinitializer.logic().CLIBackwardFlow(sys.argv[1], iter_num)
