import os
import unittest
from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy

#
# T1_Mapping
#

class T1_Mapping(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "T1_Mapping" # TODO make this more human readable by adding spaces
    self.parent.categories = ["Quantification"]
    self.parent.dependencies = []
    self.parent.contributors = ["Xiao Da (MGH), Artem Mamonov (MGH), Jayashree Kalpathy-Cramer (MGH), Andriy Fedorov (BWH)"] # replace with "Firstname Lastname (Organization)"
    self.parent.helpText = """
    Estimate effective T1 from multi-spectral FLASH MRI scans with arbitrary number of flip angles.
    """
    self.parent.acknowledgementText = """
    <img src=':Logos/QIICR.png'><br><br>
    Supported by NIH U24 CA180918, NIH U01 CA154601
    """ # replace with organization, grant and thanks.
    self.parent = parent
#
# T1_MappingWidget
#

class T1_MappingWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    #ScriptedLoadableModuleWidget.setup(self)

    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    #self.inputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.inputSelector.nodeTypes = ['vtkMRMLMultiVolumeNode']
    self.inputSelector.addAttribute( "vtkMRMLMultiVolumeNode", "LabelMap", 0 )
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

    #
    # output volume selector
    #
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.outputSelector.addAttribute( "vtkMRMLScalarVolumeNode", "LabelMap", 0 )
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = True
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the output to the algorithm." )
    parametersFormLayout.addRow("Output Volume: ", self.outputSelector)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply T1 mapping")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = True
    parametersFormLayout.addRow(self.applyButton)

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    #self.onSelect()

  def cleanup(self):
    pass

#  def onSelect(self):
#    self.applyButton.enabled = self.inputSelector.currentNode() and #self.outputSelector.currentNode()

  def onApplyButton(self):
    inputVolume = self.inputSelector.currentNode()
    outputVolume = self.outputSelector.currentNode()
    if not (inputVolume and outputVolume):
      qt.QMessageBox.critical(
          slicer.util.mainWindow(),
          'T1 Mapping', 'Input and output volumes are required for T1 Mapping')
      return
# run T1Mapping
# Create numpy array from images
    s=slicer.util.array(inputVolume.GetID())

# Create zero numpy arrays wich has the same dimension as the input images
    d=s.shape
    ssin=numpy.zeros(d,dtype=numpy.float)
    stan=numpy.zeros(d,dtype=numpy.float)

# Read info from VTK node (FA and TR)
    v=slicer.util.getNode(inputVolume.GetID())
    FA=string.split(v.GetAttribute('MultiVolume.FrameLabels'),',')
    nFrames=len(FA)
    for l in range(nFrames):
      FA[l]=float(FA[l])

    TR=str(v.GetAttribute('MultiVolume.DICOM.RepetitionTime'))
    TR=float(TR)

# convert degrees to radians
    FArdg=numpy.radians(FA)

# create output array
    o=numpy.zeros((d[0],d[1],d[2]))

# calculate the slope of the linear regression
    for i in range(nFrames):
      ssin[:,:,:,i]=s[:,:,:,i]/numpy.sin(FArdg[i])
      stan[:,:,:,i]=s[:,:,:,i]/numpy.tan(FArdg[i])
    
    for i in range(d[0]):
      for j in range(d[1]):
        for k in range(d[2]):
          [coef,res]=numpy.polyfit(stan[i,j,k,:],ssin[i,j,k,:],1)
          o[i,j,k]=coef

# compute the real T1    
    T1=TR/numpy.log(o)*(-1)
    T1=numpy.uint16(T1)
    importer = vtk.vtkImageImport()
    data_string = T1.tostring() 
    importer.CopyImportVoidPointer(data_string, len(data_string)) 
    setDataType = 'importer.SetDataScalarTypeTo' + 'UnsignedShort' + '()' 
    eval(setDataType) 
    importer.SetNumberOfScalarComponents(1) 
    importer.SetWholeExtent(0,T1.shape[2]-1,0,T1.shape[1]-1,0,T1.shape[0]-1) 
    importer.SetDataExtentToWholeExtent() 
    print importer.GetDataExtent() 
    importer.Update() 
    ijkToRAS = vtk.vtkMatrix4x4() 
    inputVolume.GetIJKToRASMatrix(ijkToRAS) 
    outputVolume.SetIJKToRASMatrix(ijkToRAS) 
    outputVolume.SetAndObserveImageData(importer.GetOutput())
    slicer.mrmlScene.AddNode(Volume)
    volumeDisplayNode = slicer.vtkMRMLScalarVolumeDisplayNode() 
    slicer.mrmlScene.AddNode(volumeDisplayNode) 
    greyColorTable = slicer.util.getNode('Grey') 
    volumeDisplayNode.SetAndObserveColorNodeID(greyColorTable.GetID())
# make the output volume appear in all the slice views
    #selectionNode = slicer.app.applicationLogic().GetSelectionNode()
    #selectionNode.SetReferenceActiveVolumeID(outputVolume.GetID())
    #slicer.app.applicationLogic().PropagateVolumeSelection(0) 

