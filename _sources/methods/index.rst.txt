.. _methods:

======================================
Methods
======================================

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

    subgraph cluster
    {
      style=rounded; color=gold;

      subgraph cluster1
      {
        style="rounded,filled"; color=skyblue;
        TPF [label="Two Phase Flow"];
        {node[shape=none, width=0, height=0, label=""] pi };
        HTE [label="Heat Transfer"];
        URL="thermalhydraulics.html";
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
        URL="neutronics.html";
        label = "Neutronics";
      }
      label = "Multi-Physics";
    }

    {edge[style=invis] {TPF -> pi -> HTE} };
    HTE -> TPF [constraint = false];
    LATT -> HOM -> CORE;
    {edge[style=invis] LATT -> REC -> CORE };
    CORE -> REC -> LATT [constraint = false];
    CORE -> HTE [constraint = false];
    TPF -> LATT [constraint = false];
    }

.. toctree::
   :maxdepth: 2

   neutronics
   thermalhydraulics
   notation



