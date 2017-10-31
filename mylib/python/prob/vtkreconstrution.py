"""
We use this script for ploting 1d vtk files with matplotlib instead
of paraview. We first use the `vtk` library to read the data,
and then we convert the data to `numpy` vector, to be ploted with
`matplotlib`.

You have to install python-vtk:
  sudo apt-get install python-vtk

The package is documented in
  http://www.vtk.org/doc/nightly/html/index.html

The function that we are using is documented in
  http://www.vtk.org/doc/nightly/html/classvtkUnstructuredGridReader.html

And the data object returned is documented here
  http://www.vtk.org/doc/nightly/html/classvtkUnstructuredGrid.html
"""

from os import listdir
from os.path import isfile, join

import sys
import numpy as np
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt

def get_data(inputf="input.vtk"):
  filename = inputf
  reader = vtkUnstructuredGridReader()
  reader.SetFileName(filename)
  reader.ReadAllVectorsOn()
  reader.ReadAllScalarsOn()
  reader.Update()

  data = reader.GetOutput()

  x = np.zeros(data.GetNumberOfPoints())
  y = np.zeros(data.GetNumberOfPoints())
  z = np.zeros(data.GetNumberOfPoints())

  for i in range(data.GetNumberOfPoints()):
    x[i],y[i],z[i] = data.GetPoint(i)

  n_arrays = data.GetPointData().GetNumberOfArrays()
  n_points = data.GetNumberOfPoints()

  array = np.zeros((n_arrays,n_points))
  names = [0]*n_arrays

  for i in range(n_arrays):
    array[i] = VN.vtk_to_numpy(data.GetPointData().GetArray(i))
    names[i] = data.GetPointData().GetArray(i).GetName()

  #n_phi_power = sum("power" in names[i] or "Phi" in names[i] for i in range(n_arrays))
  #array_phi_power = np.zeros((n_phi_power,n_points))
  #names_phi_power = [0]*n_phi_power

  array_phi_power = np.concatenate([[array[i]] for i in range(n_arrays) if "power" in names[i] or "Phi" in names[i]])
  names_phi_power = [names[i] for i in range(n_arrays) if "power" in names[i] or "Phi" in names[i]]

  #plt.rc('text', usetex=True) # this is for latex rendering
  #plt.title("$"+names[i]+"$")
  return(x,array_phi_power,names_phi_power)

def reconstruct_data(rootf="./c5g7heter1d/input"):
  file_het = rootf+".vtk"
  x_het,array_het,names_het = get_data(file_het)
  file_hom = rootf+"HomAss.vtk"
  x_hom,array_hom,names_hom = get_data(file_hom)
  mypath = "."
  #print [f for f in listdir(mypath) if isfile(join(mypath, f)) and ]
  files = [f for f in listdir(mypath) if f.endswith((".vtk")) and rootf+"Ass" in f and not "_geom" in f]
  files.sort()
  n_ass = len(files)

  all_x = [[get_data(f)[0]] for f in files]
  all_array = [[get_data(f)[1]] for f in files]
  all_names = [get_data(f)[2] for f in files]

  print all_names[0]
  print [name for name in all_names[0] if "Phi" in name]
  n_Phi = len([name for name in all_names[0] if "Phi" in name])
  print n_Phi

  ass_x = np.concatenate(all_x[i] for i in range(1,n_ass+1))
  print ass_x
  #for i in range(1,n_ass+1):


  ass_x = np.concatenate([get_data(f)[0] for f in files])
  #print ass_x
  ass_array = np.concatenate([[get_data(f)[1]] for f in files])
  #print ass_array
  ass_names = [[get_data(f)[2]] for f in files]
  #print ass_names

  array_hom_interp = np.concatenate([[np.interp(x_het, x_hom, v)] for v in array_hom],axis=0)

  array_rec = np.concatenate([[np.interp(x_het, x_hom, v)] for v in array_hom],axis=0)

  print "array_het = \n", array_het, "\n"
  print "array_hom_interp = \n", array_hom_interp, "\n"

  #print type(array_hom_interp)
  #print array_hom_interp.size
  #print type(array_het)
  #print array_het.size

  array = np.concatenate([array_het,array_hom_interp])
  names = names_het+names_hom
  return(x_het,array,names)

def plot_power(x,array,names):
  n_arrays = len(array)

  fig0 = plt.figure()
  ax0 = fig0.add_subplot(111)
  for i in range(n_arrays):
    if "power" not in names[i]:
      continue
    ax0.plot(x, array[i], 'r-', linewidth=2.0, label=names[i])
  colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired
  colors = [colormap(i) for i in np.linspace(0, 0.8,len(ax0.lines))]
  for i,j in enumerate(ax0.lines):
    j.set_color(colors[i])
  ax0.legend(loc=2)
  plt.axis('tight'); vx0,vx1,vy0,vy1 = plt.axis();
  plt.axis((vx0,vx1,0,1.1*vy1))
  plt.show()

def plot_data(x,array,names):
  n_arrays = len(array)

  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
  for i in range(n_arrays):
    if "Phi" not in names[i]:
      continue
    ax1.plot(x, array[i], 'b-', linewidth=2.0, label=names[i])
    #plt.title(names[i])
  colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired
  colors = [colormap(i) for i in np.linspace(0, 0.8,len(ax1.lines))]
  for i,j in enumerate(ax1.lines):
    j.set_color(colors[i])
  ax1.legend(loc=2)
  plt.axis('tight'); vx0,vx1,vy0,vy1 = plt.axis();
  plt.axis((vx0,vx1,0,1.1*vy1))
  plt.show()

if __name__ == "__main__":
  if len(sys.argv)==1:
    print "No file is provided, so input.vtk is the default used."
    x,array,names = reconstruct_data()
    print x.size, len(x)
    print array.size, len(array)
    print names
    plot_data(x,array,names)
  elif len(sys.argv)==2:
    print "File ",sys.argv[1]," will be printed."
    x,array,names = reconstruct_data(sys.argv[1])
    plot_data(x,array,names)
  else:
    print "Too many inputs provided. "
    "Run the script again with right number of inputs."



