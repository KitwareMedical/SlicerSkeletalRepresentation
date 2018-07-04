# -*- coding: utf-8 -*-
import numpy as np
import vtk


class spoke:
    def __init__(self, ux, uy, uz, r):
        self.U = np.array([ux, uy, uz])
        self.r = r


class hub:
    def __init__(self, x, y, z):
        self.P = np.array([x, y, z])


class atom:
    def __init__(self, hub):
        self.hub = hub

    def addSpoke(self, spoke, side):
        if side == 0:
            self.topSpoke = spoke
        elif side == 1:
            self.botSpoke = spoke
        else:
            self.crestSpoke = spoke

    def setLocation(self, row, col):
        self.row = row
        self.col = col

    def isCrest(self):
        return hasattr(self, 'crestSpoke')


class figure:
    def __init__(self, numRows, numCols):
        self.numRows = numRows
        self.numCols = numCols
        self.atoms = np.ndarray([numRows, numCols], dtype=object)

    def addAtom(self, row, col, atom):
        self.atoms[row, col] = atom

    def addAtomFromDict(self, row, col, primdict):
        px = float(primdict['x'])
        py = float(primdict['y'])
        pz = float(primdict['z'])

        newhub = hub(px, py, pz)
        newatom = atom(newhub)

        ux0 = float(primdict['ux[0]'])
        uy0 = float(primdict['uy[0]'])
        uz0 = float(primdict['uz[0]'])
        r0 = float(primdict['r[0]'])
        spoke0 = spoke(ux0, uy0, uz0, r0)
        newatom.addSpoke(spoke0, 0)

        ux1 = float(primdict['ux[1]'])
        uy1 = float(primdict['uy[1]'])
        uz1 = float(primdict['uz[1]'])
        r1 = float(primdict['r[1]'])
        spoke1 = spoke(ux1, uy1, uz1, r1)
        newatom.addSpoke(spoke1, 1)

        primtype = primdict['type']

        if primtype == 'EndPrimitive':
            ux2 = float(primdict['ux[2]'])
            uy2 = float(primdict['uy[2]'])
            uz2 = float(primdict['uz[2]'])
            r2 = float(primdict['r[2]'])
            spoke2 = spoke(ux2, uy2, uz2, r2)
            newatom.addSpoke(spoke2, 2)

        self.addAtom(row, col, newatom)


class srep:
    def __init__(self):
        self.fig = None

    def getSectionDict(self, lines, start, stop):
        section = [line.strip(';') for line in lines[start:stop]]
        sectiondict = dict(line.split(' = ') for line in section)
        return sectiondict

    def readSrepFromM3D(self, filename):
        f = open(filename, 'r')

        lines = [line.strip() for line in f.readlines()]

        figidx = lines.index('figure[0] {')
        coloridx = lines.index('color {')
        endidx = lines.index('}', coloridx)

        figparams = self.getSectionDict(lines, figidx + 1, coloridx)
        numrows = int(figparams['numRows'])
        numcols = int(figparams['numColumns'])

        fig = figure(numrows, numcols)

        primidx = endidx + 1

        for row in range(numrows):
            for col in range(numcols):
                endprimidx = lines.index('}', primidx)
                primsection = self.getSectionDict(lines, primidx + 1, endprimidx)
                fig.addAtomFromDict(row, col, primsection)
                primidx = endprimidx + 1

        self.fig = fig
