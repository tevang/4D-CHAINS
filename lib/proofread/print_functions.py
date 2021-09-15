import sys

from lib.peak import *

def unfold_peak(v):
    if type(v) == Peak:
        return v.__get_CHresons__()
    else:
        return v

def print_dict_contents(dictionary):
    """
    Print the contents of a dictionary. Can work with multi-dimensional dictionaries, too.
    :param dictionary: a dict or tree() object
    :return:
    """
    for k,v in list(dictionary.items()):
        sys.stdout.write(str(unfold_peak(k)) + " --> ")
        try:
            print_dict_contents(v)
        except AttributeError:
            if type(v) == list:
                value_string = ""
                for i in v:
                    if type(i) in [tuple, list, set]:
                        value_string += str([unfold_peak(j) for j in i]) + ", "
                    else:
                        value_string += str(unfold_peak(j)) + ", "
                print("[%s]" % value_string)
            elif type(v) == tuple:
                print(tuple([unfold_peak(j) for i in v for j in i]))
            elif type(v) == set:
                print(set([unfold_peak(j) for i in v for j in i]))
            else:
                print(unfold_peak(v))