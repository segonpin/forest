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

import sys
import numpy as np
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt

def plot_data(inputf="input.vtk"):
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

  #plt.rc('text', usetex=True) # this is for latex rendering
  #plt.title("$"+names[i]+"$") 

  for i in range(n_arrays):
    if "power" not in names[i]:
      continue
    plt.plot(x, array[i], 'r-', linewidth=2.0)
    plt.title(names[i]) 
    plt.axis('tight'); vx0,vx1,vy0,vy1 = plt.axis(); 
    plt.axis((vx0,vx1,0,1.1*vy1))
    plt.show()
    
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
    plot_data()
  elif len(sys.argv)==2:
    print "File ",sys.argv[1]," will be printed."
    plot_data(sys.argv[1])
  else:
    print "Too many inputs provided. "
    "Run the script again with right number of inputs."

