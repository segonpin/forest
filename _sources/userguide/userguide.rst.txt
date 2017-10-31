.. _userguide:

======================================
User Guide
======================================

..  toctree::
    :maxdepth: 2
   
    input
    output
    examples

The following graph represent the work-flow of the library, with more
detailed information about each module by clicking the graph.

.. graphviz::

    digraph G {
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
        URL="input.html";
        label="Input Data"
      }
      subgraph cluster
      {
        style=rounded; color=gold;
        subgraph cluster1
        {
          style="rounded,filled"; color=skyblue;
          TPF [label="Two Phase Flow"];
          {node[shape=none, width=0, height=0, label=""] pi };
          HTE [label="Heat Transfer"];
          URL="../methods/thermalhydraulics.html";
          label = "Thermal-Hydraulics";
        }
        subgraph cluster2
        {
          style="rounded,filled"; color=orange;
          LATT [label="Lattice"];
          {
            rank=same;
            HOM [label="Homogenization"];
            REC [label="Reconstruction"];
          }
          CORE [label="Core"];
          URL="../methods/neutronics.html";
          label = "Neutronics";
        }
        URL="../methods/multiphysics.html";
        label = "Multi-Physics";
      }
      subgraph cluster3
      {
        style="rounded,filled"; color=yellowgreen;
        {
          rank=same;
          OUTgraph [label="Graphical Output"];
          OUTdata [label="Data Output"];
        }
        URL="output.html";
        label="Output Results";
      }
      {edge[style=invis] {TPF -> pi -> HTE} };
      HTE -> TPF [constraint = false];
      {node[shape=none, width=0, height=0, label=""] IN};
      {edge[style=invis] {INsolver -> IN -> TPF} };
      {edge[style=invis] {INprob -> IN -> LATT} };
      INsolver -> TPF [ltail=cluster0, lhead=cluster];
      INprob -> LATT [ltail=cluster0, lhead=cluster];
      LATT -> HOM -> CORE;
      {edge[style=invis] LATT -> REC -> CORE };
      CORE -> REC -> LATT [constraint = false];
      CORE -> HTE [constraint = false];
      TPF -> LATT [constraint = false];
      {node[shape=none, width=0, height=0, label=""] OUT };
      {edge[style=invis] {CORE -> OUT -> OUTdata} };
      {edge[style=invis] {HTE -> OUT -> OUTgraph} };
      CORE -> OUTdata [ltail=cluster, lhead=cluster3];
      HTE -> OUTgraph [ltail=cluster, lhead=cluster3];
    }

Brief explanation of each module:

  - Input
  
    - **Problem data**: The data defining the problem to be solved.
    - **Solver settings**: The parameters used by the numerical methods to solve the
      problem (tolerances,  maximum number of iterations,  etc).

  - Neutronics
  
    - **Core**: Approximation of the neutron flux inside  the reactor core (modeled 
      for example with the neutron diffusion equation).
    - **Lattice**: Approximation of the neutron flux inside an assembly (modeled for
      example with the neutron transport equation).
    - **Homogenization**: This procedure is used to homogenize and condense the cross
      sections from the Lattice problem to be solved with a lower order operator.
    - **Reconstruction**: Reconstructing a detailed solution of the neutron flux from
      the coarse mesh solution.
       
  - Thermal-Hydraulics
  
    - **Heat transfer**: Here we calculate the heat that will reach the moderator taking
      into account the power generated inside the fuel.
    - **Two Phase Flow**: This module solves the thermal-hydraulics for the given problem
      using the heat provided by the heat transfer equation.
       
  - Output
  
    - **Data output**: Data related to the solution.
    - **Graphical output**: vtk file used to plot the solution (for example with Paraview).




