# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 18:12:18 2016

@author: segonpin
"""

import xml.etree.ElementTree as ET
import numpy as np

def indent(elem, level=0, more_sibs=False):
    """Indent the tree in a recursive way."""
    i = '\n'
    if level:
        i += (level-1) * "  " # for level bigger than zero we add spaces
    num_kids = len(elem)
    if num_kids:
        if not elem.text or not elem.text.strip():
            elem.text = i + '  '
            if level:
                elem.text += '  '
        count = 0
        for kid in elem:
            indent(kid, level+1, count < num_kids - 1)
            count += 1
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
            if more_sibs:
                elem.tail += '  '
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
            if more_sibs:
                elem.tail += '  '

def np2str(matrix, level=0):
    """Convert a np.ndarray into a formatted string."""
    if isinstance(matrix, (list)):
        if isinstance(matrix[0], (list)):
            strmatrix = str()
            for row in matrix:
                strmatrix += np2str(row, level+1)
            strmatrix = "\n" + strmatrix
            return strmatrix
        else:
            strmatrix = np2str(np.asarray(matrix), level)
            return strmatrix
    else:
        if len(matrix.shape) == 1:
            if isinstance(matrix[0],(int)):
                strrows = ["%2s " % number for number in matrix]
            else:
                strrows = ["%.6e " % number for number in matrix]
            strrows = ''.join(strrows)
            strrows += ";\n"
            if level==0:
                strrows = "\n" + strrows
            return strrows
        elif len(matrix.shape) == 2:
            strmatrix = str()
            for row in matrix:
                strmatrix += np2str(row, level+1)
            strmatrix = "\n" + strmatrix
            return strmatrix
        else:
            raise("matrix with more than 2 indices not allowed yet.")

#def np2str(matrix, level=0):
#    """Convert a np.ndarray into a formatted string."""
#    if len(matrix.shape) == 1:
#        strrows = ["%.6e " % number for number in matrix]
#        strrows = ''.join(strrows)
#        strrows += ";\n"
#        if level==0:
#            strrows = "\n" + strrows
#        return strrows
#    elif len(matrix.shape) == 2:
#        strmatrix = str()
#        for row in matrix:
#            strmatrix += np2str(row, level+1)
#        strmatrix = "\n" + strmatrix
#        return strmatrix
#    else:
#        raise("matrix with more than 2 indices not allowed yet.")



class XMLCombiner(object):
    def __init__(self, filenames):
        assert len(filenames) > 0, 'No filenames!'
        # save all the roots, in order, to be processed later
        self.roots = [ET.parse(f).getroot() for f in filenames]

    def combine(self):
        for r in self.roots[1:]:
            # combine each element with the first one, and update that
            self.combine_element(self.roots[0], r)
        # return the string representation
        return ET.tostring(self.roots[0])

    def combine_element(self, one, other):
        """
        This function recursively updates either the text or the children
        of an element if another element is found in `one`, or adds it
        from `other` if not found.
        """
        # Create a mapping from tag name to element, as that's what we are fltering with
        mapping = {el.tag: el for el in one}
        for el in other:
            if len(el) == 0:
                # Not nested
                try:
                    # Update the text
                    mapping[el.tag].text = el.text
                except KeyError:
                    # An element with this name is not in the mapping
                    mapping[el.tag] = el
                    # Add it
                    one.append(el)
            else:
                try:
                    # Recursively process the element, and update it in the same way
                    self.combine_element(mapping[el.tag], el)
                except KeyError:
                    # Not in the mapping
                    mapping[el.tag] = el
                    # Just add it
                    one.append(el)

def combine_names(fnames=['input.xml'], common_ext='.xml'):
    # check that all the names share 'common_ext' and remove it
    if not all ([name.endswith(common_ext) for name in fnames]):
        raise NameError('All names do not share the extension you provided')
    nameslist = [name[:-len(common_ext)] for name in fnames]
    # store the common letters among all names
    nmax = min([len(name) for name in nameslist])
    new_name = []
    for i in range(nmax):
        this_letter = [n[i] for n in nameslist]
        if this_letter.count(this_letter[0]) == len(this_letter):
            new_name.append(nameslist[0][i])
        else:
            break
    # if common root has been found it is returned, otherwise 'input'
    if len(new_name)==0:
        return 'input'+common_ext
    else:
        return ''.join(new_name)+common_ext

def combine_xml(fnames=['input.mat.xml']):
    new_root = ET.parse(fnames[0]).getroot()
    for f in fnames[1:]:
        max_mix_id = 0
        for mix in new_root.iter('mix'):
            max_mix_id = max(int(mix.attrib['id']), max_mix_id)
        shift = max_mix_id + 1
        next_root = ET.parse(f).getroot()
        for mix in next_root.iter('mix'):
            mix.attrib['id'] = str(int(mix.attrib['id']) + shift)
        new_root.extend(next_root)
    indent(new_root) # fix the indentation
    new_xml = ET.ElementTree(new_root)
    #new_xml.write(combine_names(fnames)+".mat.xml",
    #              xml_declaration=True, encoding='utf-8')
    # r.write(sys.argv[-1], encoding="iso-8859-1", xml_declaration=True)
    #return ET.tostring(new_xml.getroot())
    return new_xml

