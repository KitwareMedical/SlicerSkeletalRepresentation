import os
import numpy as np

def save_points_vtk(filename: str, points: np.ndarray):
    """
    Save a N x 3 numpy array as a VTK points file.
    """

    if len(points.shape) == 1:
        points.reshape((-1, points.shape[0]))

    assert len(points.shape) == 2

    if os.path.isfile(filename):
        os.system(f"rm {filename}")
    
    os.system(f"touch {filename}")

    writestr = ""
    writestr += "# vtk DataFile Version 3.0\n"
    writestr += "generated point file\n"
    writestr += "ASCII\n"
    writestr += "DATASET POLYDATA\n"

    numpts = points.shape[0]

    writestr += f"POINTS {numpts} float\n"

    for p in points:
        foo = ""
        for e in p:
            foo += f"{e} "

        foo = foo[:-1]
        foo += "\n"
        writestr += foo
    
    f = open(filename, "a")
    f.write(writestr)





