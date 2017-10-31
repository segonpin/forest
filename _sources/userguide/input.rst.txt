.. _input:

======================================
Input Data
======================================

.. contents:: Table of Contents
   :depth: 3


.. graphviz::

    digraph G
    {
      graph[bgcolor="transparent"];
      rankdir=TB;
      compound=true;
      edge [fontname="FreeSans",fontsize=15,labelfontname="FreeSans",labelfontsize=10];
      node [fontname="FreeSans",fontsize=15,
      shape=record,height=0.2,width=0.4,
      color="black", fillcolor="white", style="filled"];
      subgraph cluster0
      {
        style="rounded,filled"; color=lightgrey;
        {
          rank=same;
          INsolver [label="Solver Settings"];
          INprob [label="Problem Data"];
        }
        URL="\ref input";
        label="Input Data"
      }
    }

We have different input files summarized as:

  - Solver Settings (\*.settings.xml) to be decided by the user.
  
  - Problem Data (\*.geom.xml and \*.mat.xml) defining the benchmark.

The main reason for splitting the information in several files
is the encapsulation of information. A detailed explanation of
each one of the files follows below.

----------------
Problem geometry
----------------

.. figure:: ../../img/geom_cascade.png
   :width: 80 %
   
   Cascade of scales for geometry definition


The description is inside the block "geometry", which also defines the 
dimension of the problem, 

  .. code-block:: xml
  
    <geometry dim="2">
      ....
    </geometry>

Inside this module we have to specify the geometrical configuration of 
the problem through several modules defining a hierarchy of scales
by setting:

  - The **core** module, consisting on
  
    .. code-block:: xml

      <core composed="true" name="mini core">
        <nnodes> 1 1 </nnodes> <!-- xnodes ynodes znodes -->
        <length> 
          42.84; <!-- X- X+ AXIS -->
          21.42; <!-- Y- Y+ AXIS -->
        </length>
        <components> <!-- LAT or BOX or PCNW PCNE PCSW PCSE -->
          0;
        </components>
        <boundary> <!-- ALB or VOID or REFL -->
          2 2 2 2; <!-- X- X+ Y- Y+ -->
        </boundary>
      </core>
    
    - name: name of the reactor. Optional.
    - dimension: spatial dimension of the problem.
    - nnodes: number of nodes in each axial direction,
      where nodes refer to assemblies.
    - length: length of the nodes in the specific axis.
    - components: The id of the assembly that defines each node at the
      next level.
    - boundary: label setting the type of boundary condition for each 
      boundary, with 0=VOID, 2=REFLECTIVE.
      
  - The **lattices** module (optional), defining the geometry of each
    assembly, which depending on each assembly type consists of:
    
    - **grid** type: id,[name,]nnodes,components,type,length

      .. code-block:: xml
      
        <lattice id="0" type="grid" name="assembly">
          <nnodes> 8 8 </nnodes> <!-- xnodes ynodes znodes -->
          <components> <!-- PIN or BOX or PCNW PCNE PCSW PCSE -->
            0 0 0 0 0 0 0 0;
            0 1 0 1 1 0 1 0;
            0 0 0 0 0 0 0 0;
            0 1 0 1 1 0 1 0;
            0 1 0 1 1 0 1 0;
            0 0 0 0 0 0 0 0;
            0 1 0 1 1 0 1 0;
            0 0 0 0 0 0 0 0;
          </components>
          <length> 
            2.7 2.4 2.7 1.2 1.2 2.7 2.4 2.7; <!-- X- X+ AXIS -->
            2.7 2.4 2.7 1.2 1.2 2.7 2.4 2.7; <!-- Y- Y+ AXIS -->
          </length>
        </lattice>

      - id: where should we use this lattice inside the core composition.
      - name: name of the assembly. Optional.
      - nnodes: number of nodes in each axial direction.
      - components: The id identifying each node at the next level.
      - type: the way we define the assembly
      - length: length of the nodes in the specific axis.
      
    - **pin_map** type: id,[name,]nnodes,components,type,water_gap,pitch

      .. code-block:: xml
      
        <lattice id="0" type="pin_map" name="moderator assembly">
          <nnodes> 2 1 </nnodes> <!-- xnodes ynodes znodes -->
          <components> <!-- PIN or BOX or PCNW PCNE PCSW PCSE -->
            1 0;
          </components>
          <pitch>21.42</pitch> 
        </lattice>

      - id: where should we use this lattice inside the core composition.
      - name: name of the assembly. Optional.
      - nnodes: number of nodes in each axial direction.
      - components: The id identifying each node at the next level.
      - type: the way we define the assembly
      - pitch: length of node along each axis (one number if square node).
      
  - The **pins** module (optional), defining the geometry of each pin,
    which depending on each type of pin consists on:
    
    - **box** type: [name, ] materials (material id for the pin)
    
      .. code-block:: xml
      
        <pin id="0" type="box" name="moderator">
          <materials>0;</materials>
        </pin>

    - **pin** type: [name, ] materials (out/inner materials id),
      fuel_radius
      
      .. code-block:: xml
      
        <pin id="1" type="pin" name="fuel-cladding mixture">
          <fuel_radius>9.18</fuel_radius>
          <materials>0 1;</materials>
        </pin>

As an example we can show the following file:

.. literalinclude:: ../../../mylib/python/examples/2d1g_heter/2d1g_heter.geom.xml
   :language: xml
   :linenos:

-----------------
Problem materials
-----------------

As an example we can show the following file:

.. literalinclude:: ../../../mylib/python/examples/2d1g_heter/2d1g_heter.mat.xml
   :language: xml
   :linenos:

----------------
Solver Settings
----------------

As an example we can show the following file:

.. literalinclude:: ../../../mylib/python/examples/2d1g_heter/2d1g_heter.settings.xml
   :language: xml
   :linenos:


.. note::

    Different input files should be used for the different physics,
    maybe sharing the geometric information, so the \*.geom file can
    be split in two in the future, one part with the geometric
    information and other one with the material data, so the geometric
    information can be used also for the thermal-hydraulic module.
    Another option is to include the thermal-hydraulic properties
    in the \*.geom file. Still have to decide about that.







