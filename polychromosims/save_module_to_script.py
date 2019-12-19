import numpy as np

def str_extended(value):
    """
    A small helper function to convert a python object into executable code reproducing that object.

    Supported types
    ---------------
    None, bool, int, float, complex, str, list, dict, tuple, numpy.ndarray
    """
    np.set_printoptions(threshold=np.inf)
    
    def str_str(val):
        return "'"+val+"'"
    def str_list(val):
        if len(val) == 0:
            return "[]"
        ret = "["
        for v in val:
            ret += str_extended(v)+", "
        # Remove the last ", "
        return ret[:-2]+"]"
    def str_dict(val):
        if len(val) == 0:
            return "{}"
        ret = "{"
        for key in val.keys():
            ret += str_extended(key)+" : "+str_extended(val[key])+", "
        return ret[:-2]+"}"
    def str_tuple(val):
        if len(val) == 0:
            return "()"
        ret = "("
        for v in val:
            ret += str_extended(v)+", "
        # Remove the last ", "
        return ret[:-2]+")"
    def str_nparray(val):
        return "np."+repr(val)

    case_dict = { type(None) : str,
                  bool : str,
                  int : str,
                  float : str,
                  complex : str,
                  str : str_str,
                  list : str_list,
                  dict : str_dict,
                  tuple : str_tuple,
                  np.ndarray : str_nparray }

    try:
        return case_dict[type(value)](value)
    except KeyError as key:
        try:
            # Maybe it's some numpy type?
            return case_dict[type(value.item())](value)
        except:
            raise ValueError("Unsupported type: "+str(key)+" for attribute "+str(value))


def mod2py(mod, path, ignoreModules=True):
    """
    This function generates a python script containing all the values in the
    module. This is designed to print configuration modules in an
    easy-to-reload-and-inspect manner.

    Parameters
    ----------
    mod : a python module
        the module to save
    path : str
        the file to save to
    ignoreModules : bool
        skip anything that's itself a module.
        True by default.
    """

    to_write = [attr for attr in dir(mod) if not attr[:2] == "__"]

    with open(path, 'xt') as myfile:
        print("import numpy as np", file=myfile)
        print("", file=myfile)

        for attr in to_write:
            try:
                print(attr+" = "+str_extended(getattr(mod, attr)), file=myfile)
            except ValueError as VE:
                cur_type = type(getattr(mod, attr))
                if not (
                        (cur_type == type(np) and ignoreModules)
                     or (cur_type == type(mod2py) and getattr(mod, attr).__name__ == "gen_start_conf")
                    ):
                    raise VE
