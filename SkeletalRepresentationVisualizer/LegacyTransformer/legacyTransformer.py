# author: Zhiyuan Liu
# date: 2018.6.1 Happy Children's Day
# this file transform legacy srep to new srep not only about the file format, but divide crest spoke and inner spokes
from __future__ import print_function
from xml.etree.ElementTree import Element, SubElement
from ElementTree_pretty import prettify
# import xml.etree.cElementTree as etree
import vtk
import srep
import os
import sys

class legacyTransformer:
    def __init__(self):
        self.outputFile = None

    def transformLegacySrep(self, m3d_filename, outPrefix, epsilon):
        print('tranform legacy srep')
        outputUp = outPrefix + '/up.vtp'
        outputDown = outPrefix + '/down.vtp'
        outputCrest = outPrefix + '/crest.vtp'
        """
        The purpose of this part of code is to play around m3d to fit into the new s-rep format framework
        """

        up_spokes_polydata = vtk.vtkPolyData()
        down_spokes_polydata = vtk.vtkPolyData()
        crest_spokes_polydata = vtk.vtkPolyData()

        # read in m3d file
        s = srep.srep()
        s.readSrepFromM3D(m3d_filename)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetDataModeToAscii()

        """
        Given an s-rep,
        1. read in medial points
        2. set spoke vectors as point data
        3. set spoke vector length as point data as well
        4. create a header that contains meta data (e.g., number of rows, number of columns, mesh type)
        5. create xml header
        """
        nRows = s.fig.numRows
        nCols = s.fig.numCols
        # create
        root = Element('s-rep')
        nRowsXMLElement = SubElement(root, 'nRows')
        nRowsXMLElement.text = str(nRows)

        nColsXMLElement = SubElement(root, 'nCols')
        nColsXMLElement.text = str(nCols)

        meshTypeXMLElement = SubElement(root, 'meshType')
        meshTypeXMLElement.text = 'Quad'

        colorXMLElement = SubElement(root, 'color')

        redXMLElement = SubElement(colorXMLElement, 'red')
        redXMLElement.text = str(0)

        greenXMLElement = SubElement(colorXMLElement, 'green')
        greenXMLElement.text = str(0.5)

        blueXMLElement = SubElement(colorXMLElement, 'blue')
        blueXMLElement.text = str(0)

        isMeanFlagXMLElement = SubElement(root, 'isMean')
        isMeanFlagXMLElement.text = 'False'

        meanStatPathXMLElement = SubElement(root, 'meanStatPath')
        meanStatPathXMLElement.text = ''
        # later this can be done as a parameter
        upSpokeXMLElement = SubElement(root, 'upSpoke')
        upSpokeXMLElement.text = os.path.join(outputUp)

        downSpokeXMLElement = SubElement(root, 'downSpoke')
        downSpokeXMLElement.text = os.path.join(outputDown)

        crestSpokeXMLElement = SubElement(root, 'crestSpoke')
        crestSpokeXMLElement.text = os.path.join(outputCrest)

        file_handle = open(outPrefix + '/header.xml','w')
        file_handle.write(prettify(root))
        file_handle.close()

        # medial points and polygons first
        medial_points = vtk.vtkPoints()
        medial_points.SetDataTypeToDouble() # important, this fix the bug that new srep has different precision with legacy one, you can find it by runningapplyTps2NewSrep program
        medial_polys = vtk.vtkCellArray()

        # This will be curves
        crest_points = vtk.vtkPoints()
        crest_polys = vtk.vtkCellArray()
        #
        up_spoke_directions = vtk.vtkDoubleArray()
        up_spoke_directions.SetNumberOfComponents(3)
        up_spoke_directions.SetName("spokeDirection")

        up_spoke_lengths = vtk.vtkDoubleArray()
        up_spoke_lengths.SetNumberOfComponents(1)
        up_spoke_lengths.SetName("spokeLength")

        down_spoke_directions = vtk.vtkDoubleArray()
        down_spoke_directions.SetNumberOfComponents(3)
        down_spoke_directions.SetName("spokeDirection")

        down_spoke_lengths = vtk.vtkDoubleArray()
        down_spoke_lengths.SetNumberOfComponents(1)
        down_spoke_lengths.SetName("spokeLength")

        crest_spoke_directions = vtk.vtkDoubleArray()
        crest_spoke_directions.SetNumberOfComponents(3)
        crest_spoke_directions.SetName("spokeDirection")

        crest_spoke_lengths = vtk.vtkDoubleArray()
        crest_spoke_lengths.SetNumberOfComponents(1)
        crest_spoke_lengths.SetName("spokeLength")

        """
        Read in
          1. medial points,
          2. up spokeLength, up spokeDirection
          3. down spokeLength, down spokeDirection
        """
        for r in range(nRows):
            for c in range(nCols):
                current_atom = s.fig.atoms[r, c]
                current_point = current_atom.hub.P
                current_id = medial_points.InsertNextPoint(current_point)

                if r < nRows - 1 and c < nCols - 1:
                    quad = vtk.vtkQuad()
                    quad.GetPointIds().SetId(0, current_id)
                    quad.GetPointIds().SetId(1, current_id + nCols)
                    quad.GetPointIds().SetId(2, current_id + nCols + 1)
                    quad.GetPointIds().SetId(3, current_id + 1)
                    medial_polys.InsertNextCell(quad)

                current_up_spoke = current_atom.topSpoke
                current_up_spoke_direction = current_up_spoke.U
                up_spoke_directions.InsertNextTuple(current_up_spoke_direction)
                current_up_spoke_length = current_up_spoke.r
                up_spoke_lengths.InsertNextTuple1(current_up_spoke_length)

                current_down_spoke = current_atom.botSpoke
                current_down_spoke_direction = current_down_spoke.U
                down_spoke_directions.InsertNextTuple(current_down_spoke_direction)
                current_down_spoke_length = current_down_spoke.r
                down_spoke_lengths.InsertNextTuple1(current_down_spoke_length)

        """
        Construct crest vtp
          1. read points, direction, and length in clockwise
          2. write it as a curve  
        """
        # \TODO: There is a lot repeating of the same code. Need to think if I can abstract the repetition into a function
        # first row
        r = 0
        for c in range(nCols - 1):
            current_crest_atom = s.fig.atoms[r, c]
            current_crest_spoke = current_crest_atom.crestSpoke
            current_crest_spoke_direction = current_crest_spoke.U

            current_crest_point = current_crest_atom.hub.P + epsilon * current_crest_spoke_direction
            current_crest_point_id = crest_points.InsertNextPoint(current_crest_point)
            line = vtk.vtkLine()
            if current_crest_point_id > 0:
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, current_crest_point_id - 1)
                line.GetPointIds().SetId(1, current_crest_point_id)
                crest_polys.InsertNextCell(line)
            crest_spoke_directions.InsertNextTuple(current_crest_spoke_direction)
            current_crest_spoke_length = current_crest_spoke.r
            crest_spoke_lengths.InsertNextTuple1(current_crest_spoke_length)

        # last column
        c = nCols - 1
        for r in range(nRows - 1):
            current_crest_atom = s.fig.atoms[r, c]
            current_crest_spoke = current_crest_atom.crestSpoke
            current_crest_spoke_direction = current_crest_spoke.U

            current_crest_point = current_crest_atom.hub.P + epsilon * current_crest_spoke_direction
            current_crest_point_id = crest_points.InsertNextPoint(current_crest_point)
            line = vtk.vtkLine()
            if current_crest_point_id > 0:
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, current_crest_point_id - 1)
                line.GetPointIds().SetId(1, current_crest_point_id)
                crest_polys.InsertNextCell(line)
            crest_spoke_directions.InsertNextTuple(current_crest_spoke_direction)
            current_crest_spoke_length = current_crest_spoke.r
            crest_spoke_lengths.InsertNextTuple1(current_crest_spoke_length)

        # bottom row
        r = nRows - 1
        for c in range(nCols-1, 0, -1):
            current_crest_atom = s.fig.atoms[r, c]
            current_crest_spoke = current_crest_atom.crestSpoke
            current_crest_spoke_direction = current_crest_spoke.U

            current_crest_point = current_crest_atom.hub.P + epsilon * current_crest_spoke_direction
            current_crest_point_id = crest_points.InsertNextPoint(current_crest_point)
            line = vtk.vtkLine()
            if current_crest_point_id > 0:
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, current_crest_point_id - 1)
                line.GetPointIds().SetId(1, current_crest_point_id)
                crest_polys.InsertNextCell(line)
            crest_spoke_directions.InsertNextTuple(current_crest_spoke_direction)
            current_crest_spoke_length = current_crest_spoke.r
            crest_spoke_lengths.InsertNextTuple1(current_crest_spoke_length)

        # first column
        c = 0
        for r in range(nRows-1, 0, -1):
            current_crest_atom = s.fig.atoms[r, c]
            current_crest_spoke = current_crest_atom.crestSpoke
            current_crest_spoke_direction = current_crest_spoke.U

            current_crest_point = current_crest_atom.hub.P + epsilon * current_crest_spoke_direction
            current_crest_point_id = crest_points.InsertNextPoint(current_crest_point)
            line = vtk.vtkLine()
            if current_crest_point_id > 0:
                line = vtk.vtkLine()
                line.GetPointIds().SetId(0, current_crest_point_id - 1)
                line.GetPointIds().SetId(1, current_crest_point_id)
                crest_polys.InsertNextCell(line)
            crest_spoke_directions.InsertNextTuple(current_crest_spoke_direction)
            current_crest_spoke_length = current_crest_spoke.r
            crest_spoke_lengths.InsertNextTuple1(current_crest_spoke_length)
        # to complete loop
        line = vtk.vtkLine()
        line.GetPointIds().SetId(0, current_crest_point_id)
        line.GetPointIds().SetId(1, 0)
        crest_polys.InsertNextCell(line)

        polyLine = vtk.vtkPolyLine()
        polyLine.GetPointIds().SetNumberOfIds(crest_points.GetNumberOfPoints()+1)
        for i in range(crest_points.GetNumberOfPoints()):
            # polyLine->GetPointIds()->SetId(i, i);
            polyLine.GetPointIds().SetId(i, i)

        polyLine.GetPointIds().SetId(crest_points.GetNumberOfPoints(),0)
        cells = vtk.vtkCellArray()
        cells.InsertNextCell(polyLine)

        crest_spokes_polydata.SetPoints(crest_points)
        crest_spokes_polydata.SetLines(cells)
        crest_spokes_polydata.GetPointData().AddArray(crest_spoke_directions)
        crest_spokes_polydata.GetPointData().SetActiveVectors("spokeDirection")
        crest_spokes_polydata.GetPointData().AddArray(crest_spoke_lengths)
        crest_spokes_polydata.GetPointData().SetActiveScalars("spokeLength")
        writer.SetFileName(outputCrest)
        writer.SetInputData(crest_spokes_polydata)
        writer.Update()

        up_spokes_polydata.SetPoints(medial_points)
        up_spokes_polydata.SetPolys(medial_polys)
        up_spokes_polydata.GetPointData().AddArray(up_spoke_directions)
        up_spokes_polydata.GetPointData().SetActiveVectors("spokeDirection")
        up_spokes_polydata.GetPointData().AddArray(up_spoke_lengths)
        up_spokes_polydata.GetPointData().SetActiveScalars("spokeLength")
        writer.SetFileName(outputUp)
        writer.SetInputData(up_spokes_polydata)
        writer.Update()

        down_spokes_polydata.SetPoints(medial_points)
        down_spokes_polydata.SetPolys(medial_polys)
        down_spokes_polydata.GetPointData().AddArray(down_spoke_directions)
        down_spokes_polydata.GetPointData().SetActiveVectors("spokeDirection")
        down_spokes_polydata.GetPointData().AddArray(down_spoke_lengths)
        down_spokes_polydata.GetPointData().SetActiveScalars("spokeLength")
        writer.SetFileName(outputDown)
        writer.SetInputData(down_spokes_polydata)
        writer.Update()


        print('Finished transformation and save to files:')
        print(outputUp)
        print(outputDown)
        print(outputCrest)

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Usage: ' + sys.argv[0] + '[input legacy m3d file name] [output prefix] [optional: epsilon]')
        exit(-1)

    epsilon = 0.02 # movement of crest spoke away from its nearby interior point
    if len(sys.argv) == 4:
        epsilon = float(sys.argv[3])

    m3d_filename = sys.argv[1]
    outPrefix = sys.argv[2]
    
    transformLegacySrep(m3d_filename, outPrefix, epsilon)
