
.. _thermalhydraulics:

======================================
Thermal-Hydraulics
======================================

.. contents:: Table of Contents
   :depth: 3

Not ready yet.

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
      subgraph cluster1
      {
        style="rounded,filled"; color=skyblue;
        TPF [label="Two Phase Flow"];
        {node[shape=none, width=0, height=0, label=""] pi };
        HTE [label="Heat Transfer"];
        label = "Thermal-Hydraulics";
      }
      {edge[style=invis] {TPF -> pi -> HTE} };
      HTE -> TPF [constraint = false];
    }





