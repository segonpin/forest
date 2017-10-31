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
#from os.path import isfile, join

import sys
import numpy as np
from vtk import vtkUnstructuredGridReader
from vtk.util import numpy_support as VN
import matplotlib.pyplot as plt

# Different nomalizations for the reconstruction
def norm1(x,f,ass,j):
  fint = np.trapz(f[ass][j],x[ass])
  dx = x[ass][-1]-x[ass][0]
  return f[ass][j]/(fint/dx)
def norm2(x,f,ass,j):
  return f[ass][j]/np.mean(f[ass][j])
def norm3(x_hom,f_hom,x,f):
  return f*np.trapz(f_hom,x_hom)/np.trapz(f,x)

# for the moment we only work with power or Phi
def this_fun(fun_name,i):
  if fun_name=="Phi":
    new_name = fun_name+"_"+str(i)
  elif fun_name=="power":
    new_name = fun_name
  else:
    print "The function name specified is wrong"
    raise
  return new_name

# Getting the data
def get_data(inputf="input.vtk"):
  # open file and load data into reader object
  filename = inputf
  reader = vtkUnstructuredGridReader()
  reader.SetFileName(filename)
  reader.ReadAllVectorsOn()
  reader.ReadAllScalarsOn()
  reader.Update()
  # move data to data and extract the numpy vectors
  data = reader.GetOutput()
  x = np.zeros(data.GetNumberOfPoints())
  y = np.zeros(data.GetNumberOfPoints())
  z = np.zeros(data.GetNumberOfPoints())
  for i in range(data.GetNumberOfPoints()):
    x[i],y[i],z[i] = data.GetPoint(i)
  # how many vectors and with how many points each
  n_arrays = data.GetPointData().GetNumberOfArrays()
  n_points = data.GetNumberOfPoints()
  array = np.zeros((n_arrays,n_points))
  names = [0]*n_arrays
  # load the data in the numpy structures we have created
  for i in range(n_arrays):
    array[i] = VN.vtk_to_numpy(data.GetPointData().GetArray(i))
    names[i] = data.GetPointData().GetArray(i).GetName()
  return(x,array,names)

def get_fun(fun_name,x,array,names):
  n_arrays = len(names)
  array = np.concatenate([[array[i]] for i in range(n_arrays) if fun_name in names[i]])
  names = [names[i] for i in range(n_arrays) if fun_name in names[i]]
  return(x,array,names)

def reconstruct_data(fun_name="Phi",subgeom="Ass",rootf="input"):
  # get heterogeneous flux
  file_het = rootf+".vtk"
  x_het,array_het,names_het = get_fun(fun_name,*get_data(file_het))
  # get homogeneous flux
  file_hom = rootf+"Hom"+subgeom+".vtk"
  x_hom,array_hom,names_hom = get_fun(fun_name,*get_data(file_hom))
  # define condition for assembly files and get assembly-wise files
  def is_ass_file(f):
      return(f.endswith((".vtk")) and rootf+subgeom in f and not "_geom" in f)
  mypath = "."
  files = [f for f in listdir(mypath) if is_ass_file(f)]
  files.sort()
  print files
  files = ["inputAss2.vtk","inputAss1.vtk","inputAss2.vtk"]
  # read x, fields and names for each assembly
  n_ass = len(files)
  all_x = [get_fun(fun_name,*get_data(f))[0] for f in files]
  all_array = [get_fun(fun_name,*get_data(f))[1] for f in files]
  all_names = [get_fun(fun_name,*get_data(f))[2] for f in files]
  # concatenate the abcisas assembly-wise
  for i in range(1,n_ass):
    all_x[i] = all_x[i-1][-1]+all_x[i]
  ass_x = np.concatenate([all_x[i] for i in range(n_ass)])
  # concatenate the functions assembly-wise
  n_fun = sum([fun_name in name for name in all_names[0]])
  fun = [0]*n_fun
  for i in range(n_fun):
    fun[i] = np.concatenate([norm1(all_x,all_array,ass,j) \
      for ass in range(n_ass) for j in range(len(all_names[ass])) \
      if this_fun(fun_name,i) in all_names[ass][j]])
  # interpolate homogeneous and assembly-wise functions at heterogeneous mesh
  array_hom_interp = np.concatenate([[np.interp(x_het, x_hom, v)] for v in array_hom],axis=0)
  array_rec = np.concatenate([[np.interp(x_het, ass_x, v)] for v in fun],axis=0)
  names_rec = [this_fun(fun_name,i) for i in range(n_fun)]
  # reconstructing from the homogeneous and assembly-wise heterogeneous
  array_rec = [array_hom_interp[i]*array_rec[i] for i in range(n_fun)]
  array_rec = [norm3(x_hom,array_hom[i],ass_x,array_rec[i]) for i in range(n_fun)]
  names_rec = [name+"_rec" for name in names_rec]
  names_het = [name+"_het" for name in names_het]
  names_hom = [name+"_hom" for name in names_hom]
  # Put all things together and returning
  #array = [array_het,array_hom_interp]
  #names = [names_het+['-'],names_hom+['.-']]
  array = [array_het,array_rec]
  names = [names_het+['-'],names_rec+['--']]
  #array = [array_het,array_hom_interp,array_rec]
  #names = [names_het+['-'],names_hom+['.-'],names_rec+['--']]
  return(x_het,array,names)

def plot_data(x,array,names):
  n_figs = len(array)
  n_array = len(array[0])
  fig1 = plt.figure()
  ax1 = fig1.add_subplot(111)
  for f in range(n_figs):
    for i in range(n_array):
      ax1.plot(x, array[f][i], 'b'+names[f][-1], linewidth=2.0, label=names[f][i])
      #plt.title(names[i])
  colormap = plt.cm.gist_ncar #nipy_spectral, Set1,Paired
  #colors = [colormap(i) for i in np.linspace(0, 0.8,len(ax1.lines))]
  colors = [colormap(i) for i in np.linspace(0, 0.8,n_array)]
  colors = colors+colors+colors

  for i,j in enumerate(ax1.lines):
    j.set_color(colors[i])
  ax1.legend(loc=2)
  plt.axis('tight'); vx0,vx1,vy0,vy1 = plt.axis();
  plt.axis((vx0,vx1,0,1.1*vy1))
  plt.show()

if __name__ == "__main__":
  if len(sys.argv)==1:
    print "No file is provided, so input.vtk is the default used."
    plot_data(*reconstruct_data(fun_name="Phi"))
    plot_data(*reconstruct_data(fun_name="power"))
  elif len(sys.argv)==2:
    print "File ",sys.argv[1]," will be printed."
    x,array,names = reconstruct_data(sys.argv[1])
    plot_data(x,array,names)
  else:
    print "Too many inputs provided. "
    "Run the script again with right number of inputs."



