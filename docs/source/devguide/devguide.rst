.. _devguide:

======================================
Developer Guide
======================================

.. toctree::
   :maxdepth: 1

   tutorial

   
We provide an overview for the developers about the hierarchy followed
by the different abstract modules. The blocks are grouped by different 
levels or interfaces:

  - "User's interface": requiring no knowledge of programming
  
  - "Programmer's interface": requiring some knowledge of programming and
    basic knowledge of the methods used, this is the level that should be 
    modified or extended in order to add or change particular capabilities.
    
  - "Low level interface": Here we find the libraries used to facilitate
    the development of the previous block. These are advances libraries which
    are not intended to be modified, but that can be used in order to add
    state of the art algorithms and solvers implemented in the future in
    these libraries.

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
        style=rounded; color=gold;
        subgraph cluster1
        {
          style="rounded,filled"; color=lightgrey;
          {
            INsolver [label="Solver Settings"];
            INprob [label="Problem Data"];
          }
          URL="";
          label="Input Data"
        }
        subgraph cluster2
        {
          style="rounded,filled"; color=lightgrey;
          {
            OUTgraph [label="Graphical Output"];
            OUTdata [label="Data Output"];
          }
          URL="";
          label="Output Results";
        }
        #URL="";
        label = "User's Interface";  
      }
      subgraph cluster3
      {
        style=rounded; color=gold;
        subgraph cluster4
        {
          style="rounded,filled"; color=yellowgreen;
          {
            rank=same;
            PROBsource [label="Source Problem"];
            PROBeigen [label="Eigenvalue Problem"];
          }
          URL="";
          label="Problems";
        }
        subgraph cluster5
        {
          style="rounded,filled"; color=skyblue;
          EQsolver [label="Linear Solver"];
          {node[shape=none, width=0, height=0, label=""] solver0 };
          {node[shape=none, width=0, height=0, label=""] solver1 };
          EIGsolver [label="Eigenvalue Solver"];
          URL="";
          label = "Algebra";
        }
        subgraph cluster6
        {
          style="rounded,filled"; color=skyblue;
          LATT [label="Lattice"];
          HOM [label="Homogenization"];
          REC [label="Reconstruction"];
          CORE [label="Core"];
          URL="";
          label = "Physical Operators";
        }
        label = "Programmer's Interface";  
      }
      subgraph cluster7
      {
        style=rounded; color=gold;
        subgraph cluster8
        {
          style="rounded,filled"; color=orange;
          {
            DEALII [label="Deal.II"];
          }
          URL="";
          label="Modelling";
        }
        subgraph cluster9
        {
          style="rounded,filled"; color=orange;
          {
          rank=same;
          PETSc [label="PETSc"];
          SLEPc [label="SLEPc"];
          }
          URL="";
          label = "Algebra";
        }
        label = "Low level Interface";  
      }
      {edge[style=invis] {INsolver -> INprob} };
      {edge[style=invis] {OUTgraph -> OUTdata} };
      {edge[style=invis] {INprob -> PROBsource} };
      {edge[style=invis] {OUTdata -> PROBeigen} };
      {edge[style=invis] {PROBsource -> EQsolver} };
      {edge[style=invis] {PROBeigen -> LATT} };
      {edge[style=invis] {LATT -> HOM -> REC -> CORE} };
      {edge[style=invis] {EQsolver -> solver0 -> solver1 -> EIGsolver} };
      {edge[style=invis] {CORE -> DEALII} };
      {edge[style=invis] {EIGsolver -> SLEPc -> PETSc} };
    }




