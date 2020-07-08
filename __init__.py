"""
The SAMMpy package consists of a set of Python scripts to perform sensitivity
analysis with considering both parameteric uncertainty and process model 
uncertainty. 

This is the first version of SAMMpy and it is open source and any
assistance is welcomed. Please email us if you want to
contribute.

"""

__name__ = 'sammpy'
__author__ = 'Jing Yang, and Ming Ye'
__version__ = "0.1.0"
__maintainer__ = "Jing Yang"
__email__ = "mye@fsu.edu"


class model():
    """
    SAMMpy System Model Class.

    Parameters
    ----------
    name : string, optional
        Name of model. 
    frames : dict, required
        The system model framework contains the process names, the alternative 
        process model options and the model weights
    env : dict, optional
        The constant variables used in the model simulations
    pars : dict, required
        The random paramters used in the model simulations. Descriptions of 
        parameter names, dists and bounds should be given
    func : function, required
        A pre-defined system model function 
         """

    def __init__(self, name=None, frames={}, env={}, pars={}, func={}):
        self.name = name
        self.frames = frames
        self.env = env
        self.pars = pars
        self.func= func
        
   
    def sample(self, nobs, method='saltelli', seed=933090936):
        """
        Generate the random parameter values using SALib

        Parameters
        ----------
        nobs: int, required 
            Number of sample realizations 
        method: string, optional 
            The sampling method, default saltelli
        seed : int, optional
            The random seed

        Returns
        -------
        A NumPy matrix containing the model inputs using Saltelli's sampling
        scheme. The resulting matrix has size of N * D, where D is the 
        number of parameters.  
        """    
        
        sampling_methods_allowed = ['saltelli']
        if method.lower() not in sampling_methods_allowed:
            raise ValueError('Sampling method not supported: choose one of %s' 
                             %', '.join(sampling_methods_allowed))
            
        nvars = len(self.pars['names'])
        problem = self.pars.copy()
        problem['num_vars'] = nvars
        
        if method.lower()=='saltelli':
            from SALib.sample import saltelli
            values = saltelli.sample(problem, nobs, seed=seed, 
                                     calc_second_order=False)[::nvars + 2, :]
        
        return values