# -*- coding: utf-8 -*-
"""
Created on Fri Aug 12 18:14:05 2016

@author: segonpin
"""
# importing the classes for generating the input files
from geometry import Geometry
from materials import Materials
from settings import Settings

import sys
#from os import getcwd as os.getcwd
import os.path
#from os.path import relpath
# We use this library to execute an external command
from subprocess import call

class RunForest(object):
    """Class Settings."""

    def __init__(self, **kwargs):
        self.description = ("Description of the problem")
        self.prob_name = "input"
        # parse the arguments
        for prop, value in kwargs.iteritems():
            if prop == "prob_name" or prop == "description":
                setattr(self, prop, value)
            else:
                print "property " + prop + " can not be stored."
        # now we prepare the paths
        self.set_paths()

    def set_paths(self):
        # paths to some important files
        #self.proj_dir = os.path.join(os.sep, "home", "segonpin", "Codes",
        #                             "dealii-neutron-transport")
        self.proj_dir = os.path.abspath(os.path.join("..","..","..",".."))
        self.executable = os.path.join(self.proj_dir,
                                       os.path.join("testing","main.exe"))
        self.test_dir = os.path.join("mylib","python","examples","heter1d")
        self.full_test_dir = os.path.join(self.proj_dir, self.test_dir)
        # We now make everything relative to the current directory.
        self.here = os.getcwd()
        self.proj_dir = os.path.relpath(self.proj_dir, self.here)
        self.executable = os.path.relpath(self.executable, self.here)
        self.full_test_dir = os.path.relpath(self.full_test_dir, self.here)

    def set_settings(self, **kwargs):
        self.settings = Settings(prob_name=self.prob_name, **kwargs)
        self.settings.build_tree()

    def set_geometry(self, **kwargs):
        self.geometry = Geometry(description=self.description, **kwargs)
        self.geometry.build_tree()

    def set_materials(self, **kwargs):
        self.materials = Materials(description=self.description, **kwargs)
        self.materials.build_tree()

    def attach_materials(self, tree):
        self.materials = Materials(description=self.description)
        self.materials.attach_tree(tree)

    def print_inputs(self):
        # prepare the names for the files
        fsettings_xml = self.settings.input_files['settings']
        fgeom_xml = self.settings.input_files['geom']
        fmat_xml = self.settings.input_files['mat']
        # print the xml files
        self.settings.print_xml(fsettings_xml)
        self.geometry.print_xml(fgeom_xml)
        self.materials.print_xml(fmat_xml)

    def run(self):
        # if we are running a multi-thread program, we must not use multithread
        # also in openblas, as explained here
        # > https://github.com/xianyi/OpenBLAS/wiki/faq#multi-threaded
        # we can avoid it with the following call        
        os.environ['OPENBLAS_NUM_THREADS'] = '1'
        # first we print the input files. This will be removed if using SWIG
        self.print_inputs()
        # This is derived from the arguments
        #filename = PROB_NAME#os.path.join(full_test_dir, PROB_NAME)
        filename = self.settings.input_files['settings']
        if filename.endswith('.settings.xml'):
            filename = filename[:-13]
        ### THIS IS NOT USED YET ### test_log = "test_"+prob_name+".log"
        # This is derived from the arguments
        filename_err = filename + ".err"
        filename_bla = filename + ".bla"
        msg = "Running: "+self.executable+" -f "+filename+" ... "
        #sys.stdout.write('\r' + msg + ' ' * 20)
        sys.stdout.write(msg)
        sys.stdout.flush() # important
        fbla = open(filename_bla, 'w+')
        ferr = open(filename_err, 'w+')
        flag = call([self.executable, "-f", filename], stdout=fbla, stderr=ferr)
        ferr.close()
        fbla.close()
        if flag == 0: # If no problem is found, we remove auxiliary files
            os.remove(filename_bla)
            os.remove(filename_err)
            sys.stdout.write('done!\n')
            sys.stdout.flush() # important
        else:
            msg = 'Execution finished WITHOUT a satisfactory output.'
            sys.stdout.write('\n  WARNING: '+msg+'\n\n')
            sys.stdout.flush() # important
            print "  flag = ", flag

    def memory(self):
        # if we are running a multi-thread program, we must not use multithread
        # also in openblas, as explained here
        # > https://github.com/xianyi/OpenBLAS/wiki/faq#multi-threaded
        # we can avoid it with the following call        
        os.environ['OPENBLAS_NUM_THREADS'] = '1'
        # first we print the input files. This will be removed if using SWIG
        self.print_inputs()
        # This is derived from the arguments
        #filename = PROB_NAME#os.path.join(full_test_dir, PROB_NAME)
        filename = self.settings.input_files['settings']
        if filename.endswith('.settings.xml'):
            filename = filename[:-13]
        ### THIS IS NOT USED YET ### test_log = "test_"+prob_name+".log"
        # This is derived from the arguments
        filename_err = filename + ".mem"
        filename_bla = filename + ".bla"
        msg = "Checking memory for: "+self.executable+" -f "+filename+" ... "
        #sys.stdout.write('\r' + msg + ' ' * 20)
        sys.stdout.write(msg)
        sys.stdout.flush() # important
        fbla = open(filename_bla, 'w+')
        ferr = open(filename_err, 'w+')
        #flag = call(["valgrind", "--leak-check=yes", self.executable, "-f", filename], stdout=fbla, stderr=ferr)
        
        flag = call(["valgrind", "--tool=massif", self.executable, "-f", filename], stdout=fbla, stderr=ferr)
        # #If you want to use massif, you better look for a visualizer of the output, 
        # #for example the next one
        # $ sudo apt-get install massif-visualizer        
        # $ massif-visualizer
        
        #flag = call(["valgrind", "--tool=callgrind", self.executable, "-f", filename], stdout=fbla, stderr=ferr)
        # # For callgrind we have to plot the information with kcachegrind
        # $ sudo apt-get install kcachegrind        
        # $ kcachegrind callgrind.out.28047 
        
        ferr.close()
        fbla.close()
        if flag == 0: # If no problem is found, we remove auxiliary files
            #os.remove(filename_bla)
            #os.remove(filename_err)
            sys.stdout.write('done!\n')
            sys.stdout.flush() # important
        else:
            msg = 'Execution finished WITHOUT a satisfactory output.'
            sys.stdout.write('\n  WARNING: '+msg+'\n\n')
            sys.stdout.flush() # important
            print "  flag = ", flag

    def debug(self):
        # if we are running a multi-thread program, we must not use multithread
        # also in openblas, as explained here
        # > https://github.com/xianyi/OpenBLAS/wiki/faq#multi-threaded
        # we can avoid it with the following call        
        os.environ['OPENBLAS_NUM_THREADS'] = '1'
        # first we print the input files. This will be removed if using SWIG
        self.print_inputs()
        # This is derived from the arguments
        #filename = PROB_NAME#os.path.join(full_test_dir, PROB_NAME)
        filename = self.settings.input_files['settings']
        if filename.endswith('.settings.xml'):
            filename = filename[:-13]
        ### THIS IS NOT USED YET ### test_log = "test_"+prob_name+".log"
        # This is derived from the arguments
        msg = "\r  Debugging: "+self.executable+" -f "+filename+" ... "
        sys.stdout.write(msg)
        sys.stdout.flush() # important
        #import gdb
        #gdb.execute(self.executable + "-f" + filename)
        #flag = call(['gdb', "--args", self.executable + "-f" + filename])
        sys.stdout.write('done!\n')
        sys.stdout.flush() # important

    def get_keff(self):
        import xml.etree.ElementTree as ET
        tree = ET.parse(self.settings.input_files['out'])
        keff = tree.getroot()[0].text
        return float(keff)

    def clean(self):
        msg = "  Cleaning files ... "
        #sys.stdout.write('\r' + msg + ' ' * 20)
        sys.stdout.write(msg)
        sys.stdout.flush() # important

        os.remove(self.settings.input_files['settings'])
        os.remove(self.settings.input_files['geom'])
        os.remove(self.settings.input_files['mat'])
#        os.remove(self.settings.input_files['out'])
        os.remove(self.settings.input_files['mesheps'])
        os.remove(self.settings.input_files['meshvtk'])
        os.remove(self.settings.input_files['vtk'])
        #os.remove(self.settings.input_files['log'])

        #filename = PROB_NAME#os.path.join(full_test_dir, PROB_NAME)
        CLEANLOG=False
        if (CLEANLOG):
            filename = self.settings.input_files['settings']
            if filename.endswith('.settings.xml'):
                filename = filename[:-13]
            os.remove(filename+".log")
        sys.stdout.write('done!\n')###(results in "+flog+")\n"
        sys.stdout.flush() # important


def newpin(ptype='box',mat=[],fradius=0.0,name='pin'):
    if ptype in ('box'):
        return dict([
            ('type', 'box'),
            ('name', name),
            ('materials', [0, 1]),
        ])
    else:
        return dict([
            ('type', ptype),
            ('name', name),
            ('fuel_radius', fradius),
            ('materials', mat),
        ])

if __name__ == "__main__":
    problem = RunForest()
