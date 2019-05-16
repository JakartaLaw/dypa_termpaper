from modules.parameters import Parameters
import pickle
import datetime
import numpy

class Loader(object):
    """
    Loads and saves objects to .pkl files.

    NOTE: In our .gitignore we specify to not send .pkl files over git!
    This to not fuck up git. Has happened before.
    """

    @classmethod
    def save(cls, V_solution, C_solution, par, file_prefix='auto'):
        """
        Parameters
        ==========
        file_prefix = str. If auto -> creates a timestamp as prefix
        V_solution : dictionary with solution
        C_solution : dictionary with solution
        par : parameters (numba class)
        """
        if file_prefix == 'auto':
            #creates a timestamp
            file_prefix = datetime.datetime.now().isoformat()[:-7]

        V_name, C_name, par_name = cls._create_names(file_prefix)
        par_dict = cls._get_parameter_dict(par)

        pickle.dump(V_solution, open(V_name, 'wb'))
        pickle.dump(C_solution, open(C_name, 'wb'))
        pickle.dump(par_dict, open(par_name, 'wb'))

    @classmethod
    def load(cls, file_prefix):

        V_name, C_name, par_name = cls._create_names(file_prefix)

        V_solution = pickle.load(open(V_name, 'rb'))
        C_solution = pickle.load(open(C_name, 'rb'))
        par_dict = pickle.load(open(par_name, 'rb'))
        par = Parameters(**par_dict)

        return V_solution, C_solution, par

    @staticmethod
    def _create_names(file_prefix):
        V_name = file_prefix + 'V_sol.pkl'
        C_name = file_prefix + 'C_sol.pkl'
        par_name = file_prefix + 'par.pkl'

        return V_name, C_name, par_name

    @staticmethod
    def _get_parameter_dict(par):
        """turning parameters into dictionary"""

        filter_func = lambda x: x[0:2] != '__' and x != '_numba_type_'
        attrs = list(filter(filter_func, par.__dir__()))

        par_dict = {}
        for attr in attrs:
            par_dict[attr] = par.__getattribute__(attr)

        return par_dict
