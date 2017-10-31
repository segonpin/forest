from vtk import *
from vtk import vtkUnstructuredGridReader, vtkDataSetMapper, \
    vtkActor, vtkRenderer, vtkRenderWindow, vtkRenderWindowInteractor, \
    vtkLookupTable, VTK_MAJOR_VERSION
def plot_materials(problem_name):
    # The source file
    #file_name = "uGridEx.vtk"
    #file_name = "2d1g_heter.mesh.vtk"
    file_name = problem_name+".mesh.vtk"

    # Read the source file.
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(file_name)

    reader.ReadAllScalarsOn()
    #reader.SetScalarsName("UserIndex")
    reader.SetScalarsName("material")
    #reader.SetScalarsName("Phi_0")

    reader.Update() # Needed because of GetScalarRange
    output = reader.GetOutput()
    scalar_range = output.GetScalarRange()

    # Create the mapper that corresponds the objects of the vtk file
    # into graphics elements
    mapper = vtkDataSetMapper()
    if VTK_MAJOR_VERSION <= 5:
        mapper.SetInput(output)
    else:
        mapper.SetInputData(output)
    mapper.SetScalarRange(scalar_range)
    # Create the Actor
    actor = vtkActor()
    actor.SetMapper(mapper)

    # Create the Renderer
    renderer = vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(1, 1, 1) # Set background to white

    #--------------------------------------------------------------------------
    # Adding colors for categories --------------------------------------------
    cl = []
    cl.append([float(cc)/255.0 for cc in [27, 158, 119]])	# Colorbrewer Dark2
    cl.append([float(cc)/255.0 for cc in [217, 95, 2]])
    cl.append([float(cc)/255.0 for cc in [117, 112, 179]])
    cl.append([float(cc)/255.0 for cc in [231, 41, 138]])
    cl.append([float(cc)/255.0 for cc in [102, 166, 30]])
    cl.append([float(cc)/255.0 for cc in [230, 171, 2]])
    cl.append([float(cc)/255.0 for cc in [166, 118, 29]])
    cl.append([float(cc)/255.0 for cc in [102, 102, 102]])
    
    lut = vtkLookupTable()
    lutNum = len(cl)
    lut.SetNumberOfTableValues(lutNum)
    lut.Build()
    for ii,cc in enumerate(cl):
    	lut.SetTableValue(ii,cc[0],cc[1],cc[2],1.0)
    numCats = int(scalar_range[1])
    lut.SetRange(0,numCats-1)
    
    mapper.SetLookupTable(lut)
    mapper.SelectColorArray("Category_ids")
    # Until here --------------------------------------------------------------
    #--------------------------------------------------------------------------

    # Create the RendererWindow
    renderer_window = vtkRenderWindow()
    renderer_window.AddRenderer(renderer)

    # Create the RendererWindowInteractor
    interactor = vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderer_window)

    # Display the vtk_file
    interactor.Initialize()
    interactor.Start()

def plot_fluxes(problem_name):
    # The source file
    #file_name = "uGridEx.vtk"
    #file_name = "2d1g_heter.mesh.vtk"
    file_name = problem_name+".vtk"

    # Read the source file.
    reader = vtkUnstructuredGridReader()
    reader.SetFileName(file_name)

    reader.ReadAllScalarsOn()
    #reader.SetScalarsName("UserIndex")
    #reader.SetScalarsName("material")
    reader.SetScalarsName("Phi_0")

    reader.Update() # Needed because of GetScalarRange
    output = reader.GetOutput()
    scalar_range = output.GetScalarRange()

    # Create the mapper that corresponds the objects of the vtk file
    # into graphics elements
    mapper = vtkDataSetMapper()
    if VTK_MAJOR_VERSION <= 5:
        mapper.SetInput(output)
    else:
        mapper.SetInputData(output)
    mapper.SetScalarRange(scalar_range)
    # Create the Actor
    actor = vtkActor()
    actor.SetMapper(mapper)

    # Create the Renderer
    renderer = vtkRenderer()
    renderer.AddActor(actor)
    renderer.SetBackground(1, 1, 1) # Set background to white

    #--------------------------------------------------------------------------
    # Adding colors for scalars -----------------------------------------------
    
    lut = vtkLookupTable()
    refLut = vtkLookupTable()
    lut.Build()
    refLut.Build()
    for i in range(256):
        lut.SetTableValue(i, refLut.GetTableValue(255-i))
    mapper.SetLookupTable(lut)
    #mapper.SelectColorArray("Category_ids")
    # Until here --------------------------------------------------------------
    #--------------------------------------------------------------------------

    # Create the RendererWindow
    renderer_window = vtkRenderWindow()
    renderer_window.AddRenderer(renderer)

    # Create the RendererWindowInteractor
    interactor = vtkRenderWindowInteractor()
    interactor.SetRenderWindow(renderer_window)

    # Display the vtk_file
    interactor.Initialize()
    interactor.Start()
    
    
    
    
    
    
    