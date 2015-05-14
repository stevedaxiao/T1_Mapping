from __main__ import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import os, sys, getopt, tempfile, shutil, string
import unittest
import logging
import numpy

"""
Estimate effective T1 from multi-spectral FLASH MRI scans with arbitrary number of flip angles.
Contributors: Xiao Da (MGH), Artem Mamonov (MGH), Jayashree Kalpathy-Cramer (MGH), Andriy Fedorov (BWH)
"""

def usage():
    print """ ./Slicer --no-splash --no-main-window --python-script Batch_T1_Mapping.py -i dicom_dir -p T1mapping -o output_dir"""
    exit()

def main():
    rOpts = 0
    verbose = 0
    
    # the defaults
    inputdir = 'None'
    prefix = 'None'
    outputdir = 'None'

    opts, files = getopt.gnu_getopt(sys.argv[1:], "ui:p:o:",
      ["usage", "inputdir=","prefix=","outputdir="]) # parameters with : or = will expect an argument!

    for o, a in opts:
      if o in ("-u","--usage"):
        usage()
        sys.exit(0)
      elif o in ("-i", "--inputdir"):
        inputdir = a
        rOpts+=1 # fore required options
      elif o in ("-p", "--prefix"):
        prefix = a
        rOpts+=1 # fore required options  
      elif o in ("-o", "--outputdir"):
        outputdir = a
        rOpts+=1 # fore required options     
      else:
        assert False, "unhandled option"
    if rOpts < 2:
      usage()
      # expand the files into absolute paths
    inputdir = os.path.realpath(inputdir)
    if not os.path.exists(inputdir):
      print "inputdir not found"
      sys.exit()
  
    if not os.path.exists(outputdir): 
      outputdir = inputdir
      
    
    # make output name and full path
    #outputFile  = os.path.join(outputdir, prefix + ".nrrd")
    #outputFile  = os.path.realpath(outputFile)   

    # --- START WORK HERE ---
    db = slicer.dicomDatabase

    DicomFolders=os.listdir(inputdir)
    for index,DicomFolder in enumerate(DicomFolders):
      fullPath = os.path.join(inputdir,DicomFolder)
      fileNames = os.listdir(fullPath)
      for index,fileName in enumerate(fileNames):
        fileNames[index]=os.path.join(inputdir,DicomFolder,fileName)
      plugin = slicer.modules.dicomPlugins['MultiVolumeImporterPlugin']()
      print '-------------------------------------------'
      if plugin:
        loadables = plugin.examine([fileNames])
      if len(loadables) == 0:
        print('plugin failed to interpret this series')
      else:
        patientID = db.fileValue(loadables[0].files[0],'0010,0020')
        seriesDescription = db.fileValue(loadables[0].files[0],'0008,103e')
        outputPath = outputdir + '/' + patientID +  ".nrrd" 
        print outputPath
        volume = plugin.load(loadables[0])
        slicer.util.saveNode(volume,outputPath)

    #convert to numpy array"
    s=slicer.util.array(volume.GetID())
    d=s.shape
    ssin=numpy.zeros(d,dtype=numpy.float)
    stan=numpy.zeros(d,dtype=numpy.float)
    # Read info from VTK node (FA and TR)
    v=slicer.util.getNode(volume.GetID())
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

    # write output
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
    volume.GetIJKToRASMatrix(ijkToRAS) 
    outputVolume = slicer.vtkMRMLScalarVolumeNode() 
    outputVolume.SetIJKToRASMatrix(ijkToRAS) 
    outputVolume.SetAndObserveImageData(importer.GetOutput())
    slicer.mrmlScene.AddNode(outputVolume)
    outputPath = outputdir + '/' + prefix + ".nrrd" 
    print outputPath
    slicer.util.saveNode(outputVolume,outputPath)
    sys.exit()
    
if __name__ == "__main__": main()

      

