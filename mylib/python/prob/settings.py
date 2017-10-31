# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 10:22:51 2016

@author: segonpin
"""

# import ElementTree to deal with xml files
import xml.etree.ElementTree as ET
import utils

def set_root_name(rname='input'):
    if rname.endswith('.settings.xml'):
        rootname = rname[:-13]
    else:
        rootname = rname
    fnames = dict([
        ('settings', rootname + ".settings.xml"),
        ('geom', rootname + ".geom.xml"),
        ('mat', rootname + ".mat.xml"),
        ('out', rootname + ".out.xml"),
        ('mesheps', rootname + ".mesh.eps"),
        ('meshvtk', rootname + ".mesh.vtk"),
        ('vtk', rootname + ".vtk"),
    ])
    return fnames
# define the class
class Settings(object):
    """Class Settings."""

    def __init__(self, **kwargs):
        """Return a new Settings object."""
        # defaults for description
        self.description = ("Settings for the solution of the problem")
        # defaults for input files
        self.prob_name = "input"
        # defaults for output options
        self.output = dict([
            ('print_mesh', True),
            ('vtk_power', True),
            ('vtk_flux', True),
            ('vtk_mat', True),
            ('homogenize', False),
        ])
        # defaults for problem
        self.method = dict([('method', 'diffusion')])
        #self.method = dict([('method', 'transport'),
        #                    ('use_dsa', True),
        #                    ('quad', 'LevelSymType2'),
        #                    ('order', 2)])
        #self.method = dict([('method', 'transport'),
        #                    ('use_dsa', True),
        #                    ('quad', 'ChebyshevLegendre'),
        #                    ('order', [2, 2])])
        # defaults for algebra
        self.algebra = dict([('matrix_free', True),
                             ('form', 'FD'),
                             ('eig_solver', dict([('type', 'PI'),
                                                  ('tol', 1.e-7),
                                                  ('max_it', 1000)])),
                             ('mg_solver', dict([('type', 'GS'),
                                                  ('tol', 1.e-7),
                                                  ('max_it', 1)])),
                             ('inner_solver', dict([('type', 'Krylov'),
                                                    ('tol', 1.e-9),
                                                    ('max_it', 1000)]))])
        # defaults for fem settings
        self.fe_settings = dict([('degree', 0), ('n_ref', 0)])
        #self.test_dir = os.getcwd()
        # parse the arguments
        for prop, value in kwargs.iteritems():
            setattr(self, prop, value)
        # use the prob_name provided by the user to generate the input_files
        self.input_files = set_root_name(self.prob_name)

    def set_files(self, files):
        """set the files in the tree."""
        # subsection for the name of the files
        geo = ET.SubElement(files, "geom") # geom field
        geo.text = self.input_files["geom"] # field value
        mat = ET.SubElement(files, "mat") # mat field
        mat.text = self.input_files["mat"] # field value
        out = ET.SubElement(files, "out") # out field
        out.text = self.input_files["out"] # field value
        mesheps = ET.SubElement(files, "mesheps") # out field
        mesheps.text = self.input_files["mesheps"] # field value
        meshvtk = ET.SubElement(files, "meshvtk") # out field
        meshvtk.text = self.input_files["meshvtk"] # field value
        vtk = ET.SubElement(files, "vtk") # out field
        vtk.text = self.input_files["vtk"] # field value

    def set_output(self, output):
        """set the output flags in the tree."""
        #for key in self.output:
        #    field = ET.SubElement(output, key) # geom field
        #    field.text = str(self.output[key]).lower() # field value
        # subsection for the name of the files

        elem = ET.SubElement(output, "print_mesh") # print_mesh
        if "print_mesh" in self.output:
            elem.text = str(self.output["print_mesh"]).lower() # field value
        else:
            elem.text = "true"

        elem = ET.SubElement(output, "vtk_power") # vtk_power
        if "vtk_power" in self.output:
            elem.text = str(self.output["vtk_power"]).lower() # field value
        else:
            elem.text = "true"

        elem = ET.SubElement(output, "vtk_flux") # print mesh
        if "vtk_flux" in self.output:
            elem.text = str(self.output["vtk_flux"]).lower() # field value
        else:
            elem.text = "true"

        elem = ET.SubElement(output, "vtk_mat") # print mesh
        if "vtk_mat" in self.output:
            elem.text = str(self.output["vtk_mat"]).lower() # field value
        else:
            elem.text = "true"

        elem = ET.SubElement(output, "homogenize") # print mesh
        if "homogenize" in self.output:
            elem.text = str(self.output["homogenize"]).lower() # field value
        else:
            elem.text = "false"

    def set_method(self, method):
        """set the problem options in the tree."""
        # subsection for the problem discretization
        method.set('type', self.method['type'])
        if self.method['type'] == "transport":
            method.set('use_dsa', str(self.method['use_dsa']).lower())
            quad = ET.SubElement(method, "quadrature")
            quad.set('name', self.method['quad'])
            quad.set('sn', str(self.method['order']))
        elif self.method['type'] != "diffusion":
            print "the method used is not valid"


    def set_algebra(self, algebra):
        """set the algebra in the tree."""
        # subsection for the solver settings
        mat_free = ET.SubElement(algebra, "matrix_free")
        mat_free.text = str(self.algebra['matrix_free']).lower()
        eig_solver_form = ET.SubElement(algebra, "form")
        eig_solver_form.text = str(self.algebra['form']).lower()
        # eigen solver settings
        eig_solver = ET.SubElement(algebra, "eig_solver")
        eig_solver.set('type', self.algebra['eig_solver']['type'])
        tol = ET.SubElement(eig_solver, "tol")
        tol.text = str(self.algebra['eig_solver']['tol'])
        max_it = ET.SubElement(eig_solver, "max_it")
        max_it.text = str(self.algebra['eig_solver']['max_it'])
        # mg solver settings
        mg_solver = ET.SubElement(algebra, "mg_solver")
        mg_solver.set('type', self.algebra['mg_solver']['type'])
        tol = ET.SubElement(mg_solver, "tol")
        tol.text = str(self.algebra['mg_solver']['tol'])
        max_it = ET.SubElement(mg_solver, "max_it")
        max_it.text = str(self.algebra['mg_solver']['max_it'])
        # inner solver settings
        inner_solver = ET.SubElement(algebra, "inner_solver")
        inner_solver.set('type', self.algebra['inner_solver']['type'])
        tol = ET.SubElement(inner_solver, "tol")
        tol.text = str(self.algebra['inner_solver']['tol'])
        max_it = ET.SubElement(inner_solver, "max_it")
        max_it.text = str(self.algebra['inner_solver']['max_it'])

    def set_fem(self, fe_settings):
        """set the fem settings in the tree."""
        # subsection for the finite element
        fe_settings.set('degree', "%s"%(self.fe_settings['degree']))
        fe_settings.set('n_ref', "%s"%(self.fe_settings['n_ref']))

    def build_tree(self):
        """build the tree."""
        # we open the settings section
        settings = ET.Element("settings")
        # Let write a description at the beginning
        xmldescription = self.description.replace("\n", "\n      ")
        settings.append(ET.Comment(text=xmldescription))
        # subsection for the name of the files
        input_files = ET.SubElement(settings, "input_files")
        #input_files.append(ET.Comment(text="input files"))
        self.set_files(input_files)
        # subsection for the name of the files
        output = ET.SubElement(settings, "output")
        self.set_output(output)
        # subsection for the problem discretization
        method = ET.SubElement(settings, "problem")
        self.set_method(method)
        # subsection for the solver settings
        algebra = ET.SubElement(settings, "algebra")
        self.set_algebra(algebra)
        # subsection for the finite element
        fe_settings = ET.SubElement(settings, "fe_settings")
        self.set_fem(fe_settings)
        # indenting the file for the pretty printing
        utils.indent(settings)
        # we put this into the tree
        self.setts_xml = ET.ElementTree(settings)

    def print_xml(self, finp):
        self.setts_xml.write(finp, xml_declaration=True, encoding='utf-8')

def setmethod(mtype='diffusion', dsa=False, quad=['LevelSymType2',2]):
    if mtype in ('diffusion'):
        return dict([('type', mtype)])
    elif mtype in ('transport'):
        return dict([
            ('type', mtype),
            ('use_dsa', dsa),
            ('quad', quad[0]),
            ('order', quad[1])
            ])
    else:
        raise Exception(
            "Wrong type for method. Try 'diffusion' or 'transport' instead.")

def setalgebra(matfree=True, form = "FD", eig=["PI",1.e-7,1000], mg=["GS",1.e-7,10], inner=["Krylov",1.e-9,1000]):
    # By default it returns the following algebra dictionary
    # algebra = dict([
    #     ('matrix_free', True),
    #     ('eig_solver', dict([
    #         ('type', 'PI'),
    #         ('form', 'FD'),
    #         ('tol', 1.e-7),
    #         ('max_it', 1000),
    #     ])),
    #     ('inner_solver', dict([
    #         ('type', 'Krylov'),
    #         ('tol', 1.e-9),
    #         ('max_it', 1000),
    #     ]))
    # ])
    if not eig[0][0:6].lower() in ("pi", "krylov"):
        raise Exception("Wrong method. Try 'pi' or 'krylov' instead.")
    if not mg[0].lower() in ("gs", "krylov"):
        raise Exception("Wrong method. Try 'gs' or 'krylov' instead.")
    if not inner[0].lower() in ("richardson", "krylov"):
        raise Exception("Wrong method. Try 'richardson' or 'krylov' instead.")
    algebra = dict([
        ('matrix_free', matfree),
        ('form', form),
        ('eig_solver', dict([
            ('type', eig[0]),
            ('tol', eig[1]),
            ('max_it', eig[2]),
        ])),
        ('mg_solver', dict([
            ('type', mg[0]),
            ('tol', mg[1]),
            ('max_it', mg[2]),
        ])),
        ('inner_solver', dict([
            ('type', inner[0]),
            ('tol', inner[1]),
            ('max_it', inner[2]),
        ]))
    ])
    return algebra


if __name__ == "__main__":

    PROB_NAME = "xmltest.settings.xml"
    DESCRIPTION = (
        "This is just a test for the class Settings \n"
        "that can be used to create the *.setings.xml file \n"
        "needed to run FOREST.")
    METHOD = setmethod(mtype='transport',quad=['LevelSymType2', 2])
    OUTPUT = dict([("homogenize", False)])
    SETTINGS = Settings(description=DESCRIPTION,
                        prob_name=PROB_NAME, method=METHOD, output=OUTPUT)
    SETTINGS.build_tree()
    SETTINGS.print_xml(PROB_NAME)

    print setalgebra(matfree=True, eig=["PI",1.e-7,1000,"FD"], inner=["Krylov",1.e-9,1000])
    print setalgebra(matfree=True, eig=["Krylov:1",1.e-7,1000,"Standard"], inner=["Krylov",1.e-9,1000])
    print setmethod(mtype='transport', quad=['GaussLegendre', 96])
    print setmethod(mtype='transport',quad=['LevelSymType2', 2])
    print setmethod(mtype='diffusion')

