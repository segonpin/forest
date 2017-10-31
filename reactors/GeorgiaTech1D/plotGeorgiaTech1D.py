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

import numpy as np
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt

filename = "GeorgiaTech1D.vtk"
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
  plt.plot(x, array[i], 'r-o', linewidth=2.0)
  plt.title(names[i]) 
  plt.show()

