from __future__ import annotations

import logging
import json
import math
import numpy as np
import os
from pathlib import Path
from scipy.sparse.linalg import spsolve
import subprocess
from typing import Optional

import vtk
import qt
import slicer
from slicer.ScriptedLoadableModule import (
    ScriptedLoadableModule, ScriptedLoadableModuleLogic,
    ScriptedLoadableModuleWidget, ScriptedLoadableModuleTest,
)

from EvolutionarySRepUtil import ensure_requirements
ensure_requirements()

import lapy
import pyvista as pv

import sys
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),os.path.join('Resources','Libraries')))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)),os.path.join('Resources','Libraries','shanapy_private')))

from nick_srep import Srep
from smooth_cyclic_curves import curve_pd_from_points
from shanapy_private.fitting_srep_to_quasi import get_perfect_ellipsoid_and_srep_from_quasi

#
# EvolutionarySRep
#

class EvolutionarySRep(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "Evolutionary S-rep fitting"
    self.parent.categories = ["Shape Creation"]
    self.parent.dependencies = []
    self.parent.contributors = ["Jared Vicory (Kitware)"]
    self.parent.helpText = (
      "Fit s-rep to object via evolutionary approach"
    )
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = (
      "This file was originally developed by Jean-Christophe Fillion-Robin, "
      "Kitware Inc., Andras Lasso, PerkLab, and Steve Pieper, Isomics, Inc. "
      "and was partially funded by NIH grant 3P41RR013218-12S1."
    )


#
# EvolutionarySRepWidget
#

class EvolutionarySRepWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent=None):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.__init__(self, parent)
    self.logic: Optional[EvolutionarySRepLogic] = None

  def setup(self):
    """
    Called when the user opens the module the first time and the widget is initialized.
    """
    ScriptedLoadableModuleWidget.setup(self)

    # Load widget from .ui file (created by Qt Designer)
    uiWidget = slicer.util.loadUI(self.resourcePath('UI/EvolutionarySRep.ui'))
    self.layout.addWidget(uiWidget)
    self.ui = slicer.util.childWidgetVariables(uiWidget)

    self.ui.inputModel.setMRMLScene(slicer.mrmlScene)

    self.logic = EvolutionarySRepLogic(self.ui.fittingProgressBar)

    self.ui.ApplyButton.connect('clicked(bool)', self.onApplyButton)

  def onApplyButton(self):
    """
    Run processing when user clicks "Apply" button.
    """
    try:
      import deformetrica
    except ImportError:
      message = ("This module requires the use of Deformetrica which is distributed under the INRIA Non-Commercial "
      "License Agreement. This license allows for the use of Deformetrica for educational, research, or evaluation "
      "purposes only.\n\n"
      "The full license is available at https://gitlab.com/icm-institute/aramislab/deformetrica \n\n"   
      "Deformetrica is not packaged with SlicerSALT. If you agree to the license described above, "
      "click Yes below to automatically download and install Deformetrica and continue with s-rep fitting. "
      "If you do not agree, click No.")

      agree = slicer.util.confirmYesNoDisplay(message,windowTitle="Deformetrica License Notification")
      if (agree):
        import tempfile

        with tempfile.TemporaryDirectory() as tempdir:
          tempdir = Path(tempdir)

          wheel_url = "https://github.com/slicersalt/deformetrica/releases/download/v4.3/deformetrica-4.3.0-py3-none-any.whl"
          target_path = tempdir / "deformetrica-4.3.0-py3-none-any.whl"
          checksum = "SHA256:a1a32a47fb2050cc7c568e55fcb7c3fb6f54eb967b90ea9a89338b3a087bb6aa"

          slicer.util.downloadFile(wheel_url,target_path,checksum)
          slicer.util.pip_install(str(target_path))
      else:
         slicer.util.errorDisplay("Unable to continue without Deformetrica")
         return
       

    try:
      self.logic.run( self.ui.inputModel.currentNode(), Path(self.ui.outputPath.directory) )
    except Exception as e:
      slicer.util.errorDisplay("Failed to compute results: {}".format(e))
      import traceback
      traceback.print_exc()


#
# EvolutionarySRepLogic
#

class EvolutionarySRepLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, progressBar):
        """
        Called when the logic class is instantiated. Can be used for initializing member variables.
        """
        ScriptedLoadableModuleLogic.__init__(self)
        self.progressBar = progressBar

  def lapy_to_vtk(self, mesh, translation = 0.0, scale = 1.0, normalize = False):
    npd = vtk.vtkPolyData()
    temp = lapy.TriaMesh(mesh.v,mesh.t)

    if (normalize):
      center = np.mean(temp.v,axis=0)
      temp.v = temp.v - center

      temp.v = temp.v / np.linalg.norm(temp.v)

    npts = vtk.vtkPoints()
    npts.SetData(vtk.util.numpy_support.numpy_to_vtk( (scale*temp.v)+translation ))
    npd.SetPoints(npts)
    
    polys = 3*np.ones([temp.t.shape[0],4],dtype=np.int64)
    polys[:,1:] = temp.t

    npolys = vtk.vtkCellArray()
    npolys.SetCells(polys.shape[0],vtk.util.numpy_support.numpy_to_vtkIdTypeArray(polys,deep=True))
    npd.SetPolys(npolys)
    
    return npd

  def vtk_to_lapy(self, pd):
    points = vtk.util.numpy_support.vtk_to_numpy(pd.GetPoints().GetData())
    
    polys = vtk.util.numpy_support.vtk_to_numpy(pd.GetPolys().GetData())
    polys = polys.reshape([int(polys.size/4),4])[:,1:]
    
    m = lapy.TriaMesh(points,polys)
    return m
  
  def srep_to_json(self, srep, outfile):
        # Skeletal point indices are even numbered
    all_skel = []
    for i in range(0,26,2):
      skel = []
      # Spine
      # Up Spoke
      data = {}
      
      ind = i
      usp = {}
      usp["CoordinateSystem"] = "LPS"
      usp["Value"] = srep.GetPoints().GetPoint(ind)

      usd = {}
      usd["CoordinateSystem"] = "LPS"
      usd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+1)) - np.array(srep.GetPoints().GetPoint(ind))).tolist()

      us = {}
      us["SkeletalPoint"] = usp
      us["Direction"] = usd

      data = {}
      data["UpSpoke"] = us

      # Down Spoke
      dsp = {}
      dsp["CoordinateSystem"] = "LPS"
      dsp["Value"] = srep.GetPoints().GetPoint(ind+122)

      dsd = {}
      dsd["CoordinateSystem"] = "LPS"
      dsd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+122+1)) - np.array(srep.GetPoints().GetPoint(ind+122))).tolist()

      ds = {}
      ds["SkeletalPoint"] = dsp
      ds["Direction"] = dsd

      data["DownSpoke"] = ds

      skel.append(data)

      ########################################################################
      # 1st skin
      # Up Spoke
      data = {}

      ind = ind + 26
      usp = {}
      usp["CoordinateSystem"] = "LPS"
      usp["Value"] = srep.GetPoints().GetPoint(ind)

      usd = {}
      usd["CoordinateSystem"] = "LPS"
      usd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+1)) - np.array(srep.GetPoints().GetPoint(ind))).tolist()

      us = {}
      us["SkeletalPoint"] = usp
      us["Direction"] = usd

      data["UpSpoke"] = us

      # Down Spoke
      dsp = {}
      dsp["CoordinateSystem"] = "LPS"
      dsp["Value"] = srep.GetPoints().GetPoint(ind+122)

      dsd = {}
      dsd["CoordinateSystem"] = "LPS"
      dsd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+122+1)) - np.array(srep.GetPoints().GetPoint(ind+122))).tolist()

      ds = {}
      ds["SkeletalPoint"] = dsp
      ds["Direction"] = dsd

      data["DownSpoke"] = ds
      skel.append(data)

      ########################################################################
      # 2nd skin
      # Up Spoke
      data = {}

      ind = ind + 48
      usp = {}
      usp["CoordinateSystem"] = "LPS"
      usp["Value"] = srep.GetPoints().GetPoint(ind)

      usd = {}
      usd["CoordinateSystem"] = "LPS"
      usd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+1)) - np.array(srep.GetPoints().GetPoint(ind))).tolist()

      us = {}
      us["SkeletalPoint"] = usp
      us["Direction"] = usd

      data["UpSpoke"] = us

      # Down Spoke
      dsp = {}
      dsp["CoordinateSystem"] = "LPS"
      dsp["Value"] = srep.GetPoints().GetPoint(ind+122)

      dsd = {}
      dsd["CoordinateSystem"] = "LPS"
      dsd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+122+1)) - np.array(srep.GetPoints().GetPoint(ind+122))).tolist()

      ds = {}
      ds["SkeletalPoint"] = dsp
      ds["Direction"] = dsd

      data["DownSpoke"] = ds

      # Crest Spoke
      csp = {}
      csp["CoordinateSystem"] = "LPS"
      csp["Value"] = srep.GetPoints().GetPoint(ind+170)

      csd = {}
      csd["CoordinateSystem"] = "LPS"
      csd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+170+1)) - np.array(srep.GetPoints().GetPoint(ind+170))).tolist()

      cs = {}
      cs["SkeletalPoint"] = csp
      cs["Direction"] = csd

      data["CrestSpoke"] = cs
      
      skel.append(data)
      all_skel.append(skel)


    # Now we start over and go back around the other side
    for i in range(22,0,-2):
      # Spine
      # Up Spoke
      skel = []
      data = {}
      
      ind = i
      usp = {}
      usp["CoordinateSystem"] = "LPS"
      usp["Value"] = srep.GetPoints().GetPoint(ind)

      usd = {}
      usd["CoordinateSystem"] = "LPS"
      usd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+1)) - np.array(srep.GetPoints().GetPoint(ind))).tolist()

      us = {}
      us["SkeletalPoint"] = usp
      us["Direction"] = usd

      data = {}
      data["UpSpoke"] = us

      # Down Spoke
      dsp = {}
      dsp["CoordinateSystem"] = "LPS"
      dsp["Value"] = srep.GetPoints().GetPoint(ind+122)

      dsd = {}
      dsd["CoordinateSystem"] = "LPS"
      dsd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+122+1)) - np.array(srep.GetPoints().GetPoint(ind+122))).tolist()

      ds = {}
      ds["SkeletalPoint"] = dsp
      ds["Direction"] = dsd

      data["DownSpoke"] = ds

      skel.append(data)

      ########################################################################
      # 1st skin
      # Up Spoke
      data = {}

      ind = ind + 26 + 24
      usp = {}
      usp["CoordinateSystem"] = "LPS"
      usp["Value"] = srep.GetPoints().GetPoint(ind)

      usd = {}
      usd["CoordinateSystem"] = "LPS"
      usd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+1)) - np.array(srep.GetPoints().GetPoint(ind))).tolist()

      us = {}
      us["SkeletalPoint"] = usp
      us["Direction"] = usd

      data["UpSpoke"] = us

      # Down Spoke
      dsp = {}
      dsp["CoordinateSystem"] = "LPS"
      dsp["Value"] = srep.GetPoints().GetPoint(ind+122)

      dsd = {}
      dsd["CoordinateSystem"] = "LPS"
      dsd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+122+1)) - np.array(srep.GetPoints().GetPoint(ind+122))).tolist()

      ds = {}
      ds["SkeletalPoint"] = dsp
      ds["Direction"] = dsd

      data["DownSpoke"] = ds
      skel.append(data)

      ########################################################################
      # 2nd skin
      # Up Spoke
      data = {}

      ind = ind + 48
      usp = {}
      usp["CoordinateSystem"] = "LPS"
      usp["Value"] = srep.GetPoints().GetPoint(ind)

      usd = {}
      usd["CoordinateSystem"] = "LPS"
      usd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+1)) - np.array(srep.GetPoints().GetPoint(ind))).tolist()

      us = {}
      us["SkeletalPoint"] = usp
      us["Direction"] = usd

      data["UpSpoke"] = us

      # Down Spoke
      dsp = {}
      dsp["CoordinateSystem"] = "LPS"
      dsp["Value"] = srep.GetPoints().GetPoint(ind+122)

      dsd = {}
      dsd["CoordinateSystem"] = "LPS"
      dsd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+122+1)) - np.array(srep.GetPoints().GetPoint(ind+122))).tolist()

      ds = {}
      ds["SkeletalPoint"] = dsp
      ds["Direction"] = dsd

      data["DownSpoke"] = ds

      # Crest Spoke
      csp = {}
      csp["CoordinateSystem"] = "LPS"
      csp["Value"] = srep.GetPoints().GetPoint(ind+170)

      csd = {}
      csd["CoordinateSystem"] = "LPS"
      csd["Value"] = (np.array(srep.GetPoints().GetPoint(ind+170+1)) - np.array(srep.GetPoints().GetPoint(ind+170))).tolist()

      cs = {}
      cs["SkeletalPoint"] = csp
      cs["Direction"] = csd

      data["CrestSpoke"] = cs
      skel.append(data)
      all_skel.append(skel)

    j = {}
    es = {}
    es["CrestPoints"] = 24
    es["Steps"] = 2
    es["Skeleton"] = all_skel

    j["EllipticalSRep"] = es

    with(open(outfile,"w")) as outfile:
      outfile.write(json.dumps(j,indent=2))
  
  def tria_mean_curvature_flow_mod(self, tria, max_iter=30, stop_eps=1e-13, step=1.0):
    """
    Adapted from LaPy tria_mean_curvature_flow to save intermediate meshes     
    """
    meshes = []
    # pre-normalize
    trianorm = lapy.TriaMesh(tria.v, tria.t)
    trianorm.normalize_()
    meshes.append(lapy.TriaMesh(trianorm.v,trianorm.t))
    # compute fixed A
    lump = True  # for computation here and inside loop
    fem = lapy.Solver(trianorm, lump)
    a_mat = fem.stiffness
    for x in range(max_iter):
        # store last position (for delta computation below)
        vlast = trianorm.v
        # get current mass matrix and Mv
        mass = lapy.Solver.fem_tria_mass(trianorm, lump)
        mass_v = mass.dot(trianorm.v)
        
        trianorm.v = spsolve(mass + step * a_mat, mass_v)
        # normalize updated mesh
        trianorm.normalize_()
        meshes.append(lapy.TriaMesh(trianorm.v, trianorm.t))
        # compute difference
        dv = trianorm.v - vlast
        diff = np.trace(np.square(np.matmul(np.transpose(dv), mass.dot(dv))))

        if diff < stop_eps:
            print(f"Converged after {x + 1} iterations.")
            break
    return trianorm, meshes
  
  def run_backward_stage(self, start_file, stop_file, reg_dir):
    # Create deformetrica config files
    model_xml = f"""<model>

    <model-type>Registration</model-type>
    <dimension>3</dimension>

    <template>
        <object id="mesh">
            <deformable-object-type>SurfaceMesh</deformable-object-type>
            <attachment-type>Varifold</attachment-type>
            
            <noise-std>0.01</noise-std>
            
            <kernel-width>0.5</kernel-width>
            <kernel-type>torch</kernel-type>
            
            <filename>{os.path.abspath(start_file)}</filename>
        </object>

    </template>
    
    <deformation-parameters>
        <kernel-width>0.1</kernel-width>
        <kernel-type>torch</kernel-type>
        <number-of-timepoints>10</number-of-timepoints>
    </deformation-parameters>
</model>"""
    with open( reg_dir / 'model.xml', "w" ) as model_file:
        model_file.write(model_xml)

    dataset_xml = f"""<data-set>
    <subject id="subj1">
        <visit id="experiment">
            
            <filename object_id="mesh">{os.path.abspath(stop_file)}</filename>
        </visit>
    </subject>
</data-set>
"""
    with open( reg_dir / 'data_set.xml', "w" ) as dataset_file:
        dataset_file.write(dataset_xml)

    opt_xml = """<?xml version="1.0"?>
<optimization-parameters>
    <optimization-method-type>ScipyLBFGS</optimization-method-type>
</optimization-parameters>"""
    with open( reg_dir / 'optimization_parameters.xml', "w" ) as opt_file:
        opt_file.write(opt_xml)
    
    # Run deformetrica
    cwd = os.getcwd()
    os.chdir(os.path.abspath(reg_dir))
    subprocess.run(["deformetrica","estimate","model.xml","data_set.xml","-p","optimization_parameters.xml"])
    os.chdir(cwd)

  def run_shooting(self, reg_path, shooting_path, srep_file):
    reg_out_dir = reg_path / 'output'
    model_xml = f"""<model>

    <model-type>Shooting</model-type>
    <dimension>3</dimension>

    <template>
        <object id="srep">
            <deformable-object-type>PolyLine</deformable-object-type>
            <filename>{srep_file}</filename>
        </object>
    </template>

    <initial-control-points>{reg_out_dir / "DeterministicAtlas__EstimatedParameters__ControlPoints.txt"}</initial-control-points>
    <initial-momenta>{reg_out_dir / "DeterministicAtlas__EstimatedParameters__Momenta.txt"}</initial-momenta>
    
    <deformation-parameters>
        
        <kernel-width>0.1</kernel-width>
        <kernel-type>torch</kernel-type>
        <number-of-timepoints>10</number-of-timepoints>
    </deformation-parameters>
</model>"""
    with open( shooting_path /'model.xml', "w" ) as model_file:
      model_file.write(model_xml)

    # Run deformetrica
    cwd = os.getcwd()
    os.chdir(shooting_path)
    subprocess.run(["deformetrica","compute","model.xml"])

    # Run refinement
    srep = Srep()
    srep.load_polydata_file('./output/Shooting__GeodesicFlow__srep__tp_10__age_1.00.vtk')

    surf = pv.read( reg_out_dir / "DeterministicAtlas__Reconstruction__mesh__subject_subj1.vtk" )

    srep.regularize_crest(0.1)
    srep.regularize_skeleton()
    # srep.medialize_skeleton()
    srep.refine_spoke_lengths(surf)

    pv.wrap(srep.polydata).save("regularized_output_srep.vtk")

    os.chdir(cwd)

  def run(self, input_mesh, output_path:Path):
    """
    Run the processing algorithm.
    Can be used without GUI widget.
    :param inputCSV:
    :param template:
    :param output:
    """
    logging.info('Processing started')

    logging.info('Running curvature flow')
    self.progressBar.setFormat('0%: Running curvature flow')

    # Run cMCF via lapy

    sn = input_mesh.CreateDefaultStorageNode()
    sn.SetFileName(str(output_path / "input.vtp"))
    sn.WriteData(input_mesh)

    r = vtk.vtkXMLPolyDataReader()
    r.SetFileName(str(output_path / "input.vtp"))
    r.Update()

    pd = r.GetOutput()

    # Convert to lapy format 
    mesh = self.vtk_to_lapy(pd)

    # Normalize
    norm_mesh = lapy.TriaMesh(mesh.v,mesh.t)
    norm_mesh.normalize_()

    # Save transform to unnormalize later
    translation = np.mean(mesh.v,axis=0) - np.mean(norm_mesh.v,axis=0)
    scale = np.linalg.norm(mesh.v - translation) / np.linalg.norm(norm_mesh.v)

    # Run the flow
    flow_mesh, meshes = self.tria_mean_curvature_flow_mod(norm_mesh,max_iter=40,stop_eps=1e-12,step=0.001)

    mesh_dir = output_path / "meshes"
    if not os.path.isdir(mesh_dir):
        os.mkdir(mesh_dir)

    for i,m in enumerate(meshes):
      temp = self.lapy_to_vtk(m)
      pv.wrap(temp).save( mesh_dir / f"{i:02d}.vtk" )
      
    new_pd = self.lapy_to_vtk(flow_mesh)

    # Fit PE to smoothed mesh    
    self.progressBar.setFormat('10%: Fitting perfect ellipsoid')
    self.progressBar.setValue(10)
    logging.info('Fitting perfect ellipsoid')

    pe, srep, coc, vertices, crest_indices, p_crest_curve = get_perfect_ellipsoid_and_srep_from_quasi(pv.wrap(new_pd), num_spine_points=13, num_skeletal_onionskins=4)

    pe.save( mesh_dir / 'PE.vtk' )

    srep_dir = output_path / "sreps"
    if not os.path.isdir(srep_dir):
      os.mkdir(srep_dir)
    pv.wrap(srep).save( srep_dir / 'PE_srep_regular.vtk' )

    # Downsample PE mesh (these are typically huge)
    dec = vtk.vtkDecimatePro()
    dec.SetInputData(pe)
    dec.SetTargetReduction(0.9)
    dec.Update()

    PE_downsample_file = mesh_dir / 'PE_downsample.vtk'
    PE_w = vtk.vtkPolyDataWriter()
    PE_w.SetFileName(PE_downsample_file)
    PE_w.SetInputData(dec.GetOutput())
    PE_w.Update()

    # Generate crest and vertex info

    PE_crest = curve_pd_from_points(np.array([pe.GetPoint(i) for i in crest_indices]))
    PE_crest = pv.wrap(PE_crest)

    crest_dir = output_path / 'crest_curves'
    if not os.path.isdir(crest_dir):
      os.mkdir(crest_dir)

    PE_crest.save(crest_dir / "crestPE.vtk")

    # Backwards flow
    self.progressBar.setFormat('33%: Running backwards flow')
    self.progressBar.setValue(33)
    logging.info('Running backwards flow')

    # Flow from PE to first mesh
    start_file = mesh_dir / 'PE_downsample.vtk'

    start_num = len(meshes) - 1
    stop_file = mesh_dir / f'{start_num:02d}.vtk'

    reg_dir = output_path / 'registrations'
    if not os.path.isdir(reg_dir):
      os.mkdir(reg_dir)

    reg_out_dir = reg_dir / f'PE_{start_num}'
    if not os.path.isdir(reg_out_dir):
      os.mkdir(reg_out_dir)
    self.run_backward_stage(start_file, stop_file, reg_out_dir)

    used_nums = []
    used_nums.append('PE')
    used_nums.append(start_num)

    last_dir = reg_out_dir / "output"
    last_num = start_num
    next_num = math.floor(start_num / 10.0) * 10
    if start_num % 10 == 0:
      next_num = start_num - 5

    while next_num > -1:
      start_file = last_dir / "DeterministicAtlas__Reconstruction__mesh__subject_subj1.vtk"
      stop_file = mesh_dir / f'{next_num:02d}.vtk'
      
      reg_out_dir = reg_dir / f'{last_num}_{next_num}'
      if not os.path.isdir(reg_out_dir):
        os.mkdir(reg_out_dir)
      self.run_backward_stage(start_file, stop_file, reg_out_dir)

      last_dir = reg_out_dir / "output"
      last_num = next_num
      used_nums.append(last_num)

      # We skip a lot of stages > 10 because the deformations are small
      if last_num > 10:
        next_num = last_num - 5
      elif last_num == 10:
        next_num = 7
      elif last_num == 7:
        next_num = 4
      else:
        next_num = last_num - 1

    self.progressBar.setFormat('66%: Fitting s-rep via shooting')
    self.progressBar.setValue(66)
    logging.info('Fitting s-rep via shooting')

    # Run shooting
    shoot_dir = output_path / 'shootings'
    if not os.path.isdir(shoot_dir):
      os.mkdir(shoot_dir)

    dir_name = str(used_nums[0]) + "_" + str(used_nums[1])
    reg_dir_name = reg_dir / dir_name

    shoot_dir_name = shoot_dir / dir_name
    if not os.path.isdir(shoot_dir_name):
      os.mkdir(shoot_dir_name)

    self.run_shooting(reg_dir_name,shoot_dir_name, srep_dir / 'PE_srep_regular.vtk')

    last_shoot_dir = shoot_dir_name

    for i in range(1,len(used_nums)-1):
      dir_name = str(used_nums[i]) + "_" + str(used_nums[i+1])
      reg_dir_name = reg_dir / dir_name

      shoot_dir_name = shoot_dir / dir_name
      if not os.path.isdir(shoot_dir_name):
          os.mkdir(shoot_dir_name)
      self.run_shooting(reg_dir_name,shoot_dir_name,last_shoot_dir / 'regularized_output_srep.vtk')

      last_shoot_dir = shoot_dir_name

    self.progressBar.setFormat('100%: Processing completed')
    self.progressBar.setValue(100)
    logging.info('Processing completed')

    r = vtk.vtkPolyDataReader()
    r.SetFileName(last_shoot_dir / 'regularized_output_srep.vtk')
    r.Update()
    srep = r.GetOutput()

    # Undo normalization from the beginning
    srep_pts = vtk.util.numpy_support.vtk_to_numpy(srep.GetPoints().GetData())
    srep_pts = (srep_pts*scale) + translation
    srep.GetPoints().SetData(vtk.util.numpy_support.numpy_to_vtk(srep_pts))
    
    outfile = output_path / "final_srep.srep.json"
    self.srep_to_json(srep,outfile)

    slicer.util.loadNodeFromFile(outfile)

#
# EvolutionarySRepTest
#

class EvolutionarySRepTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_EvolutionarySRep1()

  def test_EvolutionarySRep1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")

    import tempfile
    import os

    logic = EvolutionarySRepLogic()

    with tempfile.TemporaryDirectory() as tempdir:
      tempdir = Path(tempdir)

      content1 = os.urandom(32)
      content2 = os.urandom(32)

      data = tempdir / 'data'
      data.mkdir()
      (data / 'file').write_bytes(content1)
      (data / 'sub').mkdir()
      (data / 'sub' / 'file').write_bytes(content2)

      output = tempdir / 'output'

      logic.run(data, output)

      self.assertTrue(output.exists())
      self.assertTrue((output / 'file').exists())
      self.assertEqual((output / 'file').read_bytes(), content1)

      self.assertTrue((output / 'sub').exists())
      self.assertTrue((output / 'file').exists())
      self.assertEqual((output / 'sub' / 'file').read_bytes(), content2)

    self.delayDisplay('Test passed')
