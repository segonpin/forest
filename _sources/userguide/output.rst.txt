.. _output:

======================================
Output Results
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
      subgraph cluster3
      {
        style="rounded,filled"; color=yellowgreen;
        {
          rank=same;
          OUTgraph [label="Graphical Output"];
          OUTdata [label="Data Output"];
        }
        label="Output Results";
      }
    }


We have different of input files to be used. A first approach can be
summarized as:

  - Data Output (\*.out)
  - Auxiliary Output (cross sections or power)??  )
  - Graphical Output (\*.vtk)

A detailed explanation of each one of the files follows below.

----------------
Standard output
----------------

The same as the program is printing to the screen

.. literalinclude:: ../../../reactors/heter2d/heter2d.log
   :language: text
   :linenos:

----------------
Data output
----------------

As an example we can show the following file:

.. literalinclude:: ../../../reactors/heter2d/heter2d.out.xml
   :language: xml
   :linenos:
   :tab-width: 2
   
----------------
Graphical output
----------------

A \*.vtk file (virtual tool-kit) to be used with a post-processing tool
(for instance \b Paraview or \b Visit).

In the following figure we can see the shape of the neutron flux for
the problem **heter2d**

.. figure:: ../../../reactors/heter2d/heter2d.png
   :width: 100 %
   
   Modeled fuel bundle
   





