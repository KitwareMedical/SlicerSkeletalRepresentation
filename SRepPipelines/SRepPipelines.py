from typing import Annotated

import vtk
import slicer
from slicer.parameterNodeWrapper import (
    Decimals,
    Default,
    Maximum,
    Minimum,
    SingleStep,
    WithinRange,
)
from slicer.ScriptedLoadableModule import *

from vtkSlicerSRepCreatorModuleLogicPython import vtkSlicerSRepCreatorLogic
from vtkSlicerSRepRefinementModuleLogicPython import vtkSlicerSRepRefinementLogic

try:
    from PipelineCreator import slicerPipeline, PipelineProgressCallback

    # Note: the values in Default() are GUI defaults only
    @slicerPipeline(name="SRep.CreateElliptical", dependencies=["SRepCreator"], categories=["SRep"])
    def srepCreator(
        model: slicer.vtkMRMLModelNode,
        numFoldPoints: Annotated[int, Minimum(0), Default(24)],
        numStepsToCrest: Annotated[int, Minimum(0), Default(2)],
        dt: Annotated[float, Decimals(5), SingleStep(1e-5), Default(0.001)],
        smoothAmount: Annotated[float, SingleStep(1e-2), Default(0.01)],
        maxIterations: Annotated[int, Default(500)],
        *,
        progressCallback: PipelineProgressCallback = PipelineProgressCallback(), 
    ) -> slicer.vtkMRMLEllipticalSRepNode:
        logic = vtkSlicerSRepCreatorLogic()
        logic.SetMRMLScene(slicer.mrmlScene)

        @vtk.calldata_type(vtk.VTK_DOUBLE)
        def progressCb(caller, event, progress):
            progressCallback.reportProgress("SRep.CreateElliptical", progress, 0, 1)
        tag = logic.AddObserver(vtk.vtkCommand.ProgressEvent, progressCb)

        val = logic.Run(model, numFoldPoints, numStepsToCrest, dt, smoothAmount, maxIterations)

        logic.RemoveObserver(tag)

        return val

    @slicerPipeline(name="SRep.Refinement", dependencies=["SRepRefinement"], categories=["SRep"])
    def srepRefinement(
        model: slicer.vtkMRMLModelNode,
        srep: slicer.vtkMRMLEllipticalSRepNode,
        interpolationLevel: Annotated[int, Minimum(0), Maximum(15), Default(3)],
        initialRegionSize: Annotated[float, SingleStep(0.1), Decimals(3), WithinRange(0.001, 15), Default(0.01)],
        finalRegionSize: Annotated[float, SingleStep(0.1), Decimals(3), WithinRange(0.0001, 15), Default(0.001)],
        maxIterations: Annotated[int, WithinRange(200, 10000), Default(2000)],
        imageMatchWeight: Annotated[float, Decimals(4), WithinRange(0, 1000), Default(0.004)],
        normalMatchWeight: Annotated[float, Decimals(4), WithinRange(0, 1000), Default(20)],
        geometricIllegalityWeight: Annotated[float, Decimals(4), WithinRange(0, 1000), Default(50)],
        *,
        progressCallback: PipelineProgressCallback = PipelineProgressCallback(),
    ) -> slicer.vtkMRMLEllipticalSRepNode:
        logic = vtkSlicerSRepRefinementLogic()
        logic.SetMRMLScene(slicer.mrmlScene)

        @vtk.calldata_type(vtk.VTK_DOUBLE)
        def progressCb(caller, event, progress):
            progressCallback.reportProgress("SRep.Refinement", progress, 0, 1)
        tag = logic.AddObserver(vtk.vtkCommand.ProgressEvent, progressCb)

        val = logic.Run(model, srep, initialRegionSize, finalRegionSize, maxIterations, interpolationLevel,
                  imageMatchWeight, normalMatchWeight, geometricIllegalityWeight)

        logic.RemoveObserver(tag)

        return val

except ImportError:
    pass

#
# SRepPipelines
#

class SRepPipelines(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "SRep Pipelines"
        self.parent.categories = ["SRep.Advanced"]
        self.parent.dependencies = ["SRep", "SRepCreator", "SRepRefinement"]
        self.parent.contributors = ["Connor Bowley (Kitware, Inc)"]
        self.parent.helpText = "This module exists to create pipelines for SRep related actions"
        self.parent.acknowledgementText = "This file was originally developed by Connor Bowley for SlicerSALT."
        self.parent.hidden = True


#
# SRepPipelinesLogic
#

class SRepPipelinesLogic(ScriptedLoadableModuleLogic):
    """This class should implement all the actual
    computation done by your module.  The interface
    should be such that other python code can import
    this class and make use of the functionality without
    requiring an instance of the Widget.
    Uses ScriptedLoadableModuleLogic base class, available at:
    https://github.com/Slicer/Slicer/blob/main/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self) -> None:
        """
        Called when the logic class is instantiated. Can be used for initializing member variables.
        """
        ScriptedLoadableModuleLogic.__init__(self)
