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
class Materials(object):
    """Class Materials."""

    def __init__(self, **kwargs):
        """Return a new Materials object."""
        # defaults for description
        self.description = ("Materials definition")
        # defaults for dimension
        self.ngroups = 1
        # defaults for mixtures
        self.mix = dict()
        # parse the arguments
        for prop, value in kwargs.iteritems():
            setattr(self, prop, value)

    def set_mix(self, materials):
        """fill the materials in the tree."""
        for idx in self.mix:
            newmix = ET.SubElement(materials, "mix")
            newmix.set('id', str(idx))
            if 'name' in self.mix[idx]:
                newmix.set('name', self.mix[idx]['name'])
            else:
                newmix.set('name', 'Mix')
            for fname in self.mix[idx]:
                if isinstance(self.mix[idx][fname], (int, long, float)):
                    f = ET.SubElement(newmix, fname)
                    f.text = '{:.6e};'.format(float(self.mix[idx][fname]))
                elif isinstance(self.mix[idx][fname], (np.ndarray)):
                    f = ET.SubElement(newmix, fname)
                    val = utils.np2str(self.mix[idx][fname])
                    f.text = val.replace("\n", "\n       ")[:-3]
                elif isinstance(self.mix[idx][fname], (list)):
                    f = ET.SubElement(newmix, fname)
                    val = utils.np2str(np.asarray(self.mix[idx][fname]))
                    f.text = val.replace("\n", "\n       ")[:-3]

    def build_tree(self):
        """build the tree."""
        # we open the settings section
        materials = ET.Element("materials")
        materials.set('ngroups', str(self.ngroups))
        # Let write a description at the beginning
        xmldescription = self.description.replace("\n", "\n      ")
        materials.append(ET.Comment(text=xmldescription))
        # Now we fill the information about the materials
        self.set_mix(materials)
        # indenting the file for the pretty printing
        utils.indent(materials)
        # we put this into the tree
        self.mat_xml = ET.ElementTree(materials)

    def attach_tree(self, tree):
        if isinstance(tree, (ET.ElementTree)):
            self.mat_xml = tree
        else:
            print "Nothing has been attached. Tree is not of the rigth class."
            #raise NameError('Tree is not of the rigth class.')

    def print_xml(self, finp):
        self.mat_xml.write(finp, xml_declaration=True, encoding='utf-8')

if __name__ == "__main__":

    PROB_NAME1 = "xmltest1"
    DESCRIPTION = (
        "This is just a test for the class Materials \n"
        "that can be used to create the *.mat.xml file \n"
        "needed to run FOREST.")
    MIX1G = dict()
    MIX1G[0] = dict([('name', 'moderator'),
                   ('SigmaT', np.array([0.37037000])),
                   ('Chi',    np.array([1.00000000])),
                   ('SigmaS', np.array([0.33400000])),
                   ('NuSigF', np.array([0.00000000]))])

    MIX1G[1] = dict([('name', 'fuel'),
                   ('SigmaT', 0.41666700),
                   ('Chi',    1.00000000),
                   ('SigmaS', 0.33400000),
                   ('NuSigF', 0.17800000)])
    MATERIALS = Materials(ngroups=1, description=DESCRIPTION, mix=MIX1G)
    MATERIALS.build_tree()
    MATERIALS.print_xml(PROB_NAME1+".mat.xml")

    PROB_NAME2 = "xmltest2"
    MATERIALS2 = Materials(ngroups=1, description=DESCRIPTION, mix=MIX1G)
    MATERIALS2.build_tree()
    MATERIALS2.print_xml(PROB_NAME2+".mat.xml")

    PROB_NAME3 = "xmltest3"
    MATERIALS3 = Materials(ngroups=1, description=DESCRIPTION, mix=MIX1G)
    MATERIALS3.build_tree()
    MATERIALS3.print_xml(PROB_NAME3+".mat.xml")

    prob_name_list = [PROB_NAME1, PROB_NAME2, PROB_NAME3]
    prob_name_list = [name+".mat.xml" for name in prob_name_list]


    #print utils.combine_names(fnames=prob_name_list,common_ext='.mat.xml')

    #print type(utils.combine_xml(prob_name_list))
    #print ET.tostring(utils.combine_xml(prob_name_list).getroot())

    #r = utils.XMLCombiner([PROB_NAME1, PROB_NAME2]).combine()
    #print '-'*20
    #print r


    #MIX2G = dict()
    #MIX2G[0] = dict([('name', 'moderator2g'),
    #               ('SigmaT', np.array([0.2203944178, 1.1188498224])),
    #               ('Chi',    np.array([1.0000000000, 0.0000000000])),
    #               ('SigmaS', np.array([[0.1928676778, 0.0000000000],
    #                                    [0.0190582000, 1.0550274224]])),
    #               ('NuSigF', np.array([0.0062056900, 0.1020900000]))])
    #MATERIALS = Materials(ngroups=2, description=DESCRIPTION, mix=MIX2G)
    #MATERIALS.build_tree(PROB_NAME)

