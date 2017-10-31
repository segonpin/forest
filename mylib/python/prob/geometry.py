# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:22:51 2016

@author: segonpin
"""

# import ElementTree to deal with xml files
import xml.etree.ElementTree as ET
import utils
import numpy as np

# define the class
class Geometry(object):
    """Class Settings."""
    tabs = '  '
    def __init__(self, **kwargs):
        """Return a new Geometry object.
        """
        # defaults for description
        self.description = (
            "Settings for the solution of the problem")
        # defaults for dimension
        self.dim = 0
        # not inside the core
        self.core = dict()
        # defaults for lattices
        self.lattices = dict()
        # defaults for pins
        self.pins = dict()
        # parse the arguments
        for prop, value in kwargs.iteritems():
            setattr(self, prop, value)

    def set_core(self, core):
        """Set the core characteristics."""
        # composed, name
        core.set('composed', str(self.core['composed']).lower())
        if "name" in self.core:
            core.set('name', self.core['name'])
        # nnodes
        nnodes = ET.SubElement(core, "nnodes")
        nnodes.text = '  '.join([str(n) for n in self.core['nnodes']])
        #length
        length = ET.SubElement(core, "length")
        val = utils.np2str(self.core['length'])
        length.text  = val.replace("\n", "\n      ")[:-2]
        #components
        components = ET.SubElement(core, "components")
        val = utils.np2str(np.asarray(self.core['components']))
        components.text  = val.replace("\n", "\n      ")[:-2]
        #components
        boundary = ET.SubElement(core, "boundary")
        boundary.text = '  '.join([str(n) for n in self.core['boundary']])
        #val = utils.np2str(np.asarray(self.core['boundary']))
        #boundary.text  = val.replace("\n", "\n      ")[:-2]

    def set_lattices(self, lattices):
        """Fill the lattices in the tree."""
        lattices.append(ET.Comment(text="lattices defining the core"))
        for ilatt in self.lattices:
            newlatt = ET.SubElement(lattices, "lattice")
            newlatt.set('id', str(ilatt))
            newlatt.set('type', self.lattices[ilatt]['type'])
            if "name" in self.lattices[ilatt]:
                newlatt.set('name', self.lattices[ilatt]['name'])
            # nnodes
            nnodes = ET.SubElement(newlatt, "nnodes")
            nnodes.text = '  '.join([str(n) for n in self.lattices[ilatt]['nnodes']])
            #components
            components = ET.SubElement(newlatt, "components")
            val = utils.np2str(np.asarray(self.lattices[ilatt]['components']))
            components.text  = val.replace("\n", "\n        ")[:-2]
            # water gap
            if "water_gap" in self.lattices[ilatt]:
                water_gap = ET.SubElement(newlatt, "water_gap")
                water_gap.text  =  '  '.join([str(n) for n in self.lattices[ilatt]['water_gap']])
            length = ET.SubElement(newlatt, "length")
            val = utils.np2str(self.lattices[ilatt]['length'])
            length.text  = val.replace("\n", "\n        ")[:-2]

    def set_pins(self, pins):
        """Fill the pins in the tree."""
        pins.append(ET.Comment(text="pins defining the lattices"))
        for ipin in self.pins:
            newpin = ET.SubElement(pins, "pin")
            newpin.set('id', str(ipin))
            newpin.set('type', self.pins[ipin]['type'])
            if "name" in self.pins[ipin]:
                newpin.set('name', self.pins[ipin]['name'])
            # name
            materials = ET.SubElement(newpin, "materials")
            materials.text  =  '  '.join([str(n) for n in self.pins[ipin]['materials']])
            # type
            # for grid
            if self.pins[ipin]['type'].lower() in ("pin", "pinnew", "pinbox"):
                # pitch
                fuel_radius = ET.SubElement(newpin, "fuel_radius")
                fuel_radius.text = '{:.6e}'.format(float(self.pins[ipin]['fuel_radius']))

    def build_tree(self):
        """build the tree."""
        # we open the settings section
        geometry = ET.Element("geometry")
        # set the dimension of the problem
        geometry.set('dim', str(self.dim)) # dimension
        # Let write a description at the beginning
        xmldescription = self.description.replace("\n", "\n      ")
        geometry.append(ET.Comment(text=xmldescription))
        # Now we fill the information about the core
        core = ET.SubElement(geometry, "core")
        self.set_core(core)
        # Now we fill the information about the lattices
        lattices = ET.SubElement(geometry, "lattices")
        self.set_lattices(lattices)
        # Now we fill the information about the pins
        pins = ET.SubElement(geometry, "pins")
        self.set_pins(pins)
        # indenting the file for the pretty printing
        utils.indent(geometry)
        # we put this into the tree
        self.geom_xml = ET.ElementTree(geometry)

    def print_xml(self, finp):
        self.geom_xml.write(finp, xml_declaration=True, encoding='utf-8')

if __name__ == "__main__":

    PROB_NAME = "xmltest.geom.xml"
    DESCRIPTION = (
        "This is just a test for the class Geometry \n"
        "that can be used to create the *.geom.xml file \n"
        "needed to run FOREST.")
    CORE1D = dict([('composed', True),
                   ('name', 'mini core'),
                   ('nnodes', [1]),
                   ('length', [[18.0]]),
                   ('components', [0]),
                   ('boundary', [0, 0])])
    print 'First problem\n'
    GEOMETRY1D = Geometry(dim=1, description=DESCRIPTION, core=CORE1D)
    GEOMETRY1D.build_tree()
    GEOMETRY1D.print_xml(PROB_NAME)

    PROB_NAME = "xmltest.geom.xml"
    DESCRIPTION = (
        "This is just a test for the class Geometry \n"
        "that can be used to create the *.geom.xml file \n"
        "needed to run FOREST.")
    CORE1D = dict([('composed', True),
                   ('name', 'mini core'),
                   ('nnodes', [1]),
                   ('length', [[18.0]]),
                   ('components', [0]),
                   ('boundary', [0, 0])])
    LATT1D = dict()
    LATT1D[0] = dict([('type', 'grid'),
                      ('name', 'assembly'),
                      ('nnodes', [7]),
                      ('components', [0, 1, 0, 1, 0, 1, 0]),
                      ('length', [2.7, 2.4, 2.7, 2.4, 2.7, 2.4, 2.7])])
    print 'Second problem\n'
    GEOMETRY1D = Geometry(dim=1, description=DESCRIPTION,
                          core=CORE1D, lattices=LATT1D)
    GEOMETRY1D.build_tree()
    GEOMETRY1D.print_xml(PROB_NAME)

    CORE2D = dict([('composed', True),
                   ('name', 'name of the core'),
                   ('nnodes', [1, 1]),
                   ('length', [[42.84], [21.42]]),
                   ('components', [1]),
                   ('boundary', [2, 2, 2, 2])])
    LATT2D = dict()
    LATT2D[0] = dict([('type', 'pin_map'),
                      ('name', 'assembly'),
                      ('nnodes', [2, 1]),
                      ('components', [1, 0]),
                      ('length', [2*[21.42]])])
    LATT2D[1] = dict([('type', 'grid'),
                      ('name', 'assembly'),
                      ('nnodes', [1]),
                      ('length', [1]),
                      ('components', [1]),
                      ('length', [[21.42]])])
    PIN2D = dict()
    PIN2D[0] = dict([('type', 'box'),
                     ('name', 'moderator'),
                     ('materials', [0])])
    PIN2D[1] = dict([('type', 'pin'),
                     ('name', 'fuel-cladding mixture'),
                     ('fuel_radius', 9.18),
                     ('materials', [0, 1])])
    print 'Third problem\n'
    GEOMETRY = Geometry(dim=2, description=DESCRIPTION,
                        core=CORE2D, lattices=LATT2D, pins=PIN2D)
    GEOMETRY.build_tree()
    GEOMETRY.print_xml(PROB_NAME)

    print 'Done.'

