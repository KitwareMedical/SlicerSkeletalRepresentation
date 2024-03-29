# This file was generated by the PipelineCreator

from typing import Annotated
from MRMLCorePython import vtkMRMLModelNode
from PipelineCreator import PipelineProgressCallback
from PipelineCreator import slicerPipeline
from Widgets.PipelineProgressBar import PipelineProgressBar
from slicer import vtkMRMLNode
from slicer.ScriptedLoadableModule import ScriptedLoadableModule
from slicer.ScriptedLoadableModule import ScriptedLoadableModuleLogic
from slicer.ScriptedLoadableModule import ScriptedLoadableModuleWidget
from slicer.parameterNodeWrapper import createGui
from slicer.parameterNodeWrapper import isParameterPack
from slicer.parameterNodeWrapper import parameterNodeWrapper
from slicer.parameterNodeWrapper import parameterPack
from slicer.util import VTKObservationMixin
from typing import Optional
from vtkSlicerSRepModuleMRMLPython import vtkMRMLEllipticalSRepNode
import pickle
import qt
import slicer
from slicer.parameterNodeWrapper import (
    Decimals,
    Default,
    Maximum,
    Minimum,
    SingleStep,
    WithinRange,
)


#
# SRepCreateAndRefine
#

class SRepCreateAndRefine(ScriptedLoadableModule):
    """Uses ScriptedLoadableModule base class, available at:
    https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
    """

    def __init__(self, parent):
        ScriptedLoadableModule.__init__(self, parent)
        self.parent.title = "SRepCreateAndRefine"
        self.parent.categories = ["PipelineModules"]
        self.parent.dependencies = ["SRepCreator", "SRepRefinement"]
        self.parent.contributors = ["Connor Bowley (Kitware, Inc)", "Jared Vicory (Kitware, Inc)", "Harald Scheirich (Kitware, Inc)", "PipelineCreator"]
        self.parent.helpText = "Pipeline module combining s-reps creation and refinement. This module was created by the PipelineCreator."
        self.parent.acknowledgementText = "This module was created by the PipelineCreator."


#
# SRepCreateAndRefineInputs
#

@parameterPack
class SRepCreateAndRefineInputs:
    createEllipticalInputModel: vtkMRMLModelNode
    createEllipticalNumFoldPoints: Annotated[int, Minimum(0), Default(24)]
    createEllipticalNumStepsToCrest: Annotated[int, Minimum(0), Default(2)]
    createEllipticalDt: Annotated[float, Decimals(5), SingleStep(1e-5), Default(0.001)]
    createEllipticalSmoothAmount: Annotated[float, SingleStep(1e-2), Default(0.01)]
    createEllipticalMaxIterations: Annotated[int, Default(500)]
    refinementInterpolationLevel: Annotated[int, Minimum(0), Maximum(15), Default(3)]
    refinementInitialRegionSize: Annotated[float, SingleStep(0.1), Decimals(3), WithinRange(0.001, 15), Default(0.01)]
    refinementFinalRegionSize: Annotated[float, SingleStep(0.1), Decimals(3), WithinRange(0.0001, 15), Default(0.001)]
    refinementMaxIterations: Annotated[int, WithinRange(200, 10000), Default(2000)]
    refinementImageMatchWeight: Annotated[float, Decimals(4), WithinRange(0, 1000), Default(0.004)]
    refinementNormalMatchWeight: Annotated[float, Decimals(4), WithinRange(0, 1000), Default(20)]
    refinementGeometricIllegalityWeight:  Annotated[float, Decimals(4), WithinRange(0, 1000), Default(50)]


#
# SRepCreateAndRefineOutputs
#

@parameterPack
class SRepCreateAndRefineOutputs:
    RefinementResult: vtkMRMLEllipticalSRepNode


#
# SRepCreateAndRefineParameterNode
#

@parameterNodeWrapper
class SRepCreateAndRefineParameterNode:
    inputs: SRepCreateAndRefineInputs
    outputs: SRepCreateAndRefineOutputs

#
# SRepCreateAndRefineWidget
#

class SRepCreateAndRefineWidget(ScriptedLoadableModuleWidget, VTKObservationMixin):
    def __init__(self, parent=None):
        self.logic = None
        self._parameterNode = None
        self._parameterNodeGuiTag = None
        ScriptedLoadableModuleWidget.__init__(self, parent)

    def setup(self):
        ScriptedLoadableModuleWidget.setup(self)
        self.logic = SRepCreateAndRefineLogic()
        self.paramWidget = createGui(SRepCreateAndRefineParameterNode)
        self.paramWidget.setMRMLScene(slicer.mrmlScene)
        self.runButton = qt.QPushButton("Run")
        self.progressBar = PipelineProgressBar()

        self.layout.addWidget(self.paramWidget)
        self.layout.addWidget(self.runButton)
        self.layout.addWidget(self.progressBar)
        self.layout.addStretch()

        # These connections ensure that we update parameter node when scene is closed
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.StartCloseEvent, self.onSceneStartClose)
        self.addObserver(slicer.mrmlScene, slicer.mrmlScene.EndCloseEvent, self.onSceneEndClose)

        # Connect the run button
        self.runButton.clicked.connect(self._onRun)

        # Make sure parameter node is initialized (needed for module reload)
        self.initializeParameterNode()

    def cleanup(self) -> None:
        """
        Called when the application closes and the module widget is destroyed.
        """
        self.removeObservers()

    def enter(self) -> None:
        """
        Called each time the user opens this module.
        """
        # Make sure parameter node exists and observed
        self.initializeParameterNode()

    def exit(self) -> None:
        """
        Called each time the user opens a different module.
        """
        # Do not react to parameter node changes (GUI will be updated when the user enters into the module)
        if self._parameterNode:
            self._parameterNode.disconnectGui(self._parameterNodeGuiTag)
            self._parameterNodeGuiTag = None

    def onSceneStartClose(self, caller, event) -> None:
        """
        Called just before the scene is closed.
        """
        # Parameter node will be reset, do not use it anymore
        self.setParameterNode(None)

    def onSceneEndClose(self, caller, event) -> None:
        """
        Called just after the scene is closed.
        """
        # If this module is shown while the scene is closed then recreate a new parameter node immediately
        if self.parent.isEntered:
            self.initializeParameterNode()

    def initializeParameterNode(self) -> None:
        """
        Ensure parameter node exists and observed.
        """
        # Parameter node stores all user choices in parameter values, node selections, etc.
        # so that when the scene is saved and reloaded, these settings are restored.

        self.setParameterNode(SRepCreateAndRefineParameterNode(self.logic.getParameterNode()))

    def setParameterNode(self, inputParameterNode: Optional[SRepCreateAndRefineParameterNode]) -> None:
        """
        Set and observe parameter node.
        Observation is needed because when the parameter node is changed then the GUI must be updated immediately.
        """

        if self._parameterNode:
            self._parameterNode.disconnectGui(self._parameterNodeGuiTag)
        self._parameterNode = inputParameterNode
        if self._parameterNode:
            self._parameterNodeGuiTag = self._parameterNode.connectGui(self.paramWidget)

    def _copyNode(self, src, dest):
        # Clones src into dest, but keeps dest's display and storage nodes, if any
        # If neither src nor dest has display nodes, the default are created
        if src is not None and dest is not None:
            name = dest.GetName()
            if dest.IsA('vtkMRMLDisplayableNode'):
                displayNodesIDs = [dest.GetNthDisplayNodeID(n) for n in range(dest.GetNumberOfDisplayNodes())]
                storageNodesIDs = [dest.GetNthStorageNodeID(n) for n in range(dest.GetNumberOfStorageNodes())]

            dest.Copy(src)
            dest.SetName(name)

            if dest.IsA('vtkMRMLDisplayableNode'):
                dest.RemoveAllDisplayNodeIDs()
                for n, displayNodeID in enumerate(displayNodesIDs):
                    dest.SetAndObserveNthDisplayNodeID(n, displayNodeID)
                for n, storageNodeID in enumerate(storageNodesIDs):
                    dest.SetAndObserveNthStorageNodeID(n, storageNodeID)

    def _copyParameterPack(self, from_, to):
        for paramName in from_.allParameters:
            if isinstance(from_.getValue(paramName), vtkMRMLNode):
                self._copyNode(from_.getValue(paramName), to.getValue(paramName))
            elif isParameterPack(from_.getValue(paramName)):
                self._copyParameterPack(from_.getValue(paramName), to.getValue(paramName))
            else:
                to.setValue(paramName, from_.getValue(paramName))

    def _removeNodes(self, item):
        if isinstance(item, vtkMRMLNode):
            slicer.mrmlScene.RemoveNode(item)
        elif isParameterPack(item):
            for paramName in item.allParameters.keys():
                self._removeNodes(item.getValue(paramName))

    def _onRun(self):
        outputValue = self.logic.run(
            createEllipticalInputModel=self._parameterNode.inputs.createEllipticalInputModel,
            createEllipticalNumFoldPoints=self._parameterNode.inputs.createEllipticalNumFoldPoints,
            createEllipticalNumStepsToCrest=self._parameterNode.inputs.createEllipticalNumStepsToCrest,
            createEllipticalDt=self._parameterNode.inputs.createEllipticalDt,
            createEllipticalSmoothAmount=self._parameterNode.inputs.createEllipticalSmoothAmount,
            createEllipticalMaxIterations=self._parameterNode.inputs.createEllipticalMaxIterations,
            refinementInterpolationLevel=self._parameterNode.inputs.refinementInterpolationLevel,
            refinementInitialRegionSize=self._parameterNode.inputs.refinementInitialRegionSize,
            refinementFinalRegionSize=self._parameterNode.inputs.refinementFinalRegionSize,
            refinementMaxIterations=self._parameterNode.inputs.refinementMaxIterations,
            refinementImageMatchWeight=self._parameterNode.inputs.refinementImageMatchWeight,
            refinementNormalMatchWeight=self._parameterNode.inputs.refinementNormalMatchWeight,
            refinementGeometricIllegalityWeight=self._parameterNode.inputs.refinementGeometricIllegalityWeight,
            progress_callback=self.progressBar.getProgressCallback())

        # Copy the output. Use CopyContent for nodes and do a normal copy for non-nodes.
        # For parameterPacks, need to recurse into them though so CopyContent can be used for
        # node members.
        if isinstance(outputValue, SRepCreateAndRefineOutputs):
            self._copyParameterPack(outputValue, self._parameterNode.outputs)
        elif isParameterPack(outputValue):
            # A parameter pack, but not the output one
            paramName = next(iter(self._parameterNode.outputs.allParameters.keys()))
            self._copyParameterPack(outputValue, self._parameterNode.outputs.getValue(paramName))
        elif isinstance(outputValue, vtkMRMLNode):
            # if the output is not a parameter pack, there is only one output
            paramName = next(iter(self._parameterNode.outputs.allParameters.keys()))
            self._copyNode(outputValue, self._parameterNode.outputs.getValue(paramName))
        else:
            # single value that is not a node
            paramName = next(iter(self._parameterNode.outputs.allParameters.keys()))
            self._parameterNode.outputs.setValue(paramName, outputValue)

        self._removeNodes(outputValue)


#
# SRepCreateAndRefineLogic
#

def _nodeReferencedBy(node, listOfNodes):
    for option in listOfNodes:
        if option is not None:
            roles = []
            option.GetNodeReferenceRoles(roles)
            for role in roles:
                if option.HasNodeReferenceID(role, node.GetID()):
                    return True
    return False

class SRepCreateAndRefineLogic(ScriptedLoadableModuleLogic):
    def __init__(self):
        ScriptedLoadableModuleLogic.__init__(self)

    @staticmethod
    @slicerPipeline(name="SRep.CreateAndRefine", dependencies=['SRepCreator', 'SRepRefinement'], categories=['SRep'])
    def run(createEllipticalInputModel: vtkMRMLModelNode, 
            createEllipticalNumFoldPoints: Annotated[int, Minimum(0), Default(24)], 
            createEllipticalNumStepsToCrest: Annotated[int, Minimum(0), Default(2)],
            createEllipticalDt: Annotated[float, Decimals(5), SingleStep(1e-5), Default(0.001)], 
            createEllipticalSmoothAmount: Annotated[float, SingleStep(1e-2), Default(0.01)],
            createEllipticalMaxIterations: Annotated[int, Default(500)],
            refinementInterpolationLevel: Annotated[int, Minimum(0), Maximum(15), Default(3)],
            refinementInitialRegionSize: Annotated[float, SingleStep(0.1), Decimals(3), WithinRange(0.001, 15), Default(0.01)],
            refinementFinalRegionSize: Annotated[float, SingleStep(0.1), Decimals(3), WithinRange(0.0001, 15), Default(0.001)],
            refinementMaxIterations: Annotated[int, WithinRange(200, 10000), Default(2000)],
            refinementImageMatchWeight: Annotated[float, Decimals(4), WithinRange(0, 1000), Default(0.004)],
            refinementNormalMatchWeight: Annotated[float, Decimals(4), WithinRange(0, 1000), Default(20)],
            refinementGeometricIllegalityWeight: Annotated[float, Decimals(4), WithinRange(0, 1000), Default(50)],
            *,
            progress_callback: PipelineProgressCallback = PipelineProgressCallback(),
            delete_intermediate_nodes: bool=True) -> vtkMRMLEllipticalSRepNode:
        progress_callback.reportProgress("", 0, 0, 2)
        # declare needed variables so they exist in the finally clause
        step_1_SRep_CreateElliptical_return = None
        step_2_SRep_Refinement_return = None

        try:
            # step 1 - SRep.CreateElliptical
            progress_callback.reportProgress("SRep.CreateElliptical", 0, 0, 2)
            function_1_SRep_CreateElliptical = pickle.loads(b'\x80\x04\x95!\x00\x00\x00\x00\x00\x00\x00\x8c\rSRepPipelines\x94\x8c\x0bsrepCreator\x94\x93\x94.')
            step_1_SRep_CreateElliptical_return = function_1_SRep_CreateElliptical(
                dt=createEllipticalDt,
                maxIterations=createEllipticalMaxIterations,
                model=createEllipticalInputModel,
                numFoldPoints=createEllipticalNumFoldPoints,
                numStepsToCrest=createEllipticalNumStepsToCrest,
                smoothAmount=createEllipticalSmoothAmount,
            progressCallback=progress_callback.getSubCallback(0, 2))

            # step 2 - SRep.Refinement
            progress_callback.reportProgress("SRep.Refinement", 0, 1, 2)
            function_2_SRep_Refinement = pickle.loads(b'\x80\x04\x95$\x00\x00\x00\x00\x00\x00\x00\x8c\rSRepPipelines\x94\x8c\x0esrepRefinement\x94\x93\x94.')
            step_2_SRep_Refinement_return = function_2_SRep_Refinement(
                finalRegionSize=refinementFinalRegionSize,
                geometricIllegalityWeight=refinementGeometricIllegalityWeight,
                imageMatchWeight=refinementImageMatchWeight,
                initialRegionSize=refinementInitialRegionSize,
                interpolationLevel=refinementInterpolationLevel,
                maxIterations=refinementMaxIterations,
                model=createEllipticalInputModel,
                normalMatchWeight=refinementNormalMatchWeight,
                srep=step_1_SRep_CreateElliptical_return,
            progressCallback=progress_callback.getSubCallback(1, 2))

            pass
        finally:
            if delete_intermediate_nodes:
                trueReturns = [step_2_SRep_Refinement_return]
                slicer.mrmlScene.RemoveNode(step_1_SRep_CreateElliptical_return)

        progress_callback.reportProgress("", 0, 2, 2)

        return step_2_SRep_Refinement_return
