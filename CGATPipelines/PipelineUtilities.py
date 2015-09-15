"""
PipelineUtilities.py - helper functions for CGAT pipelines
=============================================================

Purpose
-------

To make excuting database commands and files IO as trivial and (TODO:
robust) as possible possible. Before changing the default behaviour of
any of these functions please discuss!. Better to add new functions
and then deprecate.

.. note::

   This module needs refactoring as it implements
   some functionality found elsewhere such as the
   Database commands and text IO commands.

Reference
---------


"""
import CGAT.IOTools as IOTools
import pickle


def write(outfile, lines, header=False):
    ''' expects [[[line1-field1],[line1-field2 ] ],... ]'''
    handle = IOTools.openFile(outfile, "w")

    if header:
        handle.write("\t".join([str(title) for title in header]) + "\n")

    for line in lines:
        handle.write("\t".join([str(field) for field in line]) + "\n")

    handle.close()

# ########################## Biomart Access


def biomart_iterator(attributes,
                     host,
                     biomart,
                     dataset,
                     filters=None,
                     values=None
                     ):
    ''' Modified from pipeline biomart... '''

    # Use PipelineBiomant.biomart_iterator
    raise NotImplementedError(
        "deprecated, use Biomart.biomart_iterator")

# ########################## File Utitilies


def txtToDict(filename, key=None, sep="\t"):
    '''make a dictionary from a text file keyed
    on the specified column '''

    # Please see function in IOTools: readDict()
    count = 0
    result = {}
    valueidx, keyidx = False, False
    field_names = []

    with open(filename, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            if count == 0:
                fieldn = 0
                for rawfield in line.split(sep):
                    field = rawfield.strip()
                    if field == key:
                        keyidx = fieldn
                    field_names.append(field)
                    fieldn += 1

                if not keyidx:
                    raise ValueError("key name not found in header")
                # if not valueidx:
                #   raise ValueError(
                #     "value name not found in header")
            else:
                fields = [x.strip() for x in line.split(sep)]
                fieldn = 0
                thiskey = fields[keyidx]
                result[thiskey] = {}
                for field in fields:
                    if fieldn == keyidx:
                        pass
                    else:
                        colkey = field_names[fieldn]
                        result[thiskey][colkey] = field
                    fieldn += 1
            count += 1

    return(result)

# ############################## Objects


def save(file_name, obj):
    '''dump a python object to a file using pickle'''
    with open(file_name, "wb") as pkl_file:
        pickle.dump(obj, pkl_file)
    return


def load(file_name):
    '''retrieve a pickled python object from a file'''
    with open(file_name, "r") as pkl_file:
        data = pickle.load(pkl_file)
    return data
