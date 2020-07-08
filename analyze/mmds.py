# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 20:47:34 2020

@author: Jing
"""
import numpy as np
import numba as nb
import itertools
from ..util import results

def analyze(model, Y, print_to_console=True):
    """
    Perform variance-based sensitivty analysis for each process.

    Parameters
    ----------
    model : object
        The model defined in the sammpy
    Y : numpy.array
        A NumPy array containing the model outputs
    print_to_console : bool
        Print results directly to console (default False)
   
    Returns
    ----------
    Returns a dictionary with keys 'PSK', 'PSTK', where
    each entry is a list of size of the number of process.
    """
    # Number of sample realizations
    obs = Y.shape[1]

    # Number of process and process models
    npros = len(model.frames['names'])
    
    # Creat a dict to store the results
    S = create_si_dict(npros)
    
    # Perfrom the difference-based process sensitivty anlaysis
    if print_to_console:
        print('Runing MMDS difference-based process sensitivy analysis...') 
    MMDS = mmds_mean_var(model, Y)
    
    # Save results to the dict
    for i in range(npros):
            S['mean'][i] = MMDS[0, i]
            S['variance'][i] = MMDS[1, i]

    # Print results to console
    if print_to_console:
        print_indices(model, S)
        
    return S

def mmds_mean_var(model, Y):
    '''
    Perform difference-based process sensitivty anlaysis 

    Returns
    ========
    numpy.array : A 2D numpy array storing the mean and variance process 
    sensitivy index for each process.
    '''
    # Number of sample realizations
    N = Y.shape[1]
    
    # Number of process and process models
    npros = len(model.frames['names'])
    npros_models_list = [len(model.frames['options'][i]) for i in range(npros)]
    
    # Initilize a array to store the results
    MMDS = np.zeros((2, npros))
    
    for ipros in range(npros):
        # locate the axis of the computed process
        axis_loc = ipros * 2
        
        # determine the number of process models
        n_K_models = npros_models_list[ipros]
        n_non_K_models = int(np.product(npros_models_list) / n_K_models)
        
        # swapaxes the computed process to head:
        tmp_Y = np.rollaxis(np.rollaxis(Y, axis_loc, 0), axis_loc + 1, 1)
        
        # reshape the ouput in to 4D-array
        for axis in range(2, len(tmp_Y.shape)):
            if tmp_Y.shape[axis] == N:
                tmp_Y = np.rollaxis(tmp_Y, axis, len(tmp_Y.shape))
        tmp_Y = tmp_Y.reshape((n_K_models, N, n_non_K_models, N**(npros - 1)))
        
        # compute the model weithts for process ~K and K
        PM_K = np.array(model.frames['weights'][ipros])
        weights = model.frames['weights'].copy()
        weights.pop(ipros)
        PM_non_K = np.array([np.product(np.array((list(i)))) for i in itertools.product(*weights)])
        
        @nb.njit
        def numba_mean_var(n_K_models, n_non_K_models, N, npros, PM_K, PM_non_K, tmp_Y):
            # compute the first-order process sensitivty index
            tmp_D = np.zeros((N, N))
            D_tc_Y = np.zeros((n_non_K_models, N**(npros - 1), n_K_models, n_K_models))
            D2_tc_Y = np.zeros((n_non_K_models, N**(npros - 1), n_K_models, n_K_models)) 
            D_c_Y = np.zeros((n_non_K_models, N**(npros - 1),  n_K_models))
            D2_c_Y = np.zeros((n_non_K_models, N**(npros - 1), n_K_models))
            D_tb_Y = np.zeros((n_non_K_models, N**(npros - 1)))
            D2_tb_Y = np.zeros((n_non_K_models, N**(npros - 1)))
            D_b_Y = np.zeros(n_non_K_models)
            D2_b_Y = np.zeros(n_non_K_models)
                       
            for i in range(int(n_non_K_models)):
                for j in range(N**(npros - 1)):
                    for k1 in range(n_K_models):
                        for k2 in range(n_K_models):
                            for l1 in range(N):
                                for l2 in range(N):
                                    tmp_D[l1, l2] = abs(tmp_Y[k1, l1, i, j] - tmp_Y[k2, l2, i, j])
                                    
                            D_tc_Y[i, j, k1, k2] = np.mean(tmp_D[:, :])
                            D2_tc_Y[i, j, k1, k2] = np.mean(tmp_D[:, :]**2)
                        
                        D_c_Y[i, j, k1] = np.sum(D_tc_Y[i, j, k1, :] * PM_K)
                        D2_c_Y[i, j, k1] = np.sum(D2_tc_Y[i, j, k1, :] * PM_K)
                    
                    D_tb_Y[i, j] = np.sum(D_c_Y[i, j, :] * PM_K)
                    D2_tb_Y[i, j] = np.sum(D2_c_Y[i, j, :] * PM_K)
                    
                D_b_Y[i] = np.mean(D_tb_Y[i, :])
                D2_b_Y[i] = np.mean(D2_tb_Y[i, :])
                
            D_a_Y = np.sum(D_b_Y[:] * PM_non_K)
            D2_a_Y = np.sum(D2_b_Y[:] * PM_non_K)
            
            mean = D_a_Y
            var = D2_a_Y - D_a_Y**2
            
            return mean, var
        
        mean, var = numba_mean_var(n_K_models, n_non_K_models, N, npros, PM_K, PM_non_K, tmp_Y)
        MMDS[0, ipros] = mean
        MMDS[1, ipros] = var
        
    return MMDS

def create_si_dict(npros):
    # initialize empty dict to store process sensitivity indices
    S = results.ResultDict((k, np.zeros(npros)) for k in ('mean', 'variance'))
    return S
            
def print_indices(model, S):
    # Output to console
    title = 'Process'
    names = model.frames['names']
    npros = npros = len(model.frames['names'])
    print('%s \t mean \t variance' % title)

    for i in range(npros):
        print('%s \t %.4f \t %.4f' % (names[i], S['mean'][i], S['variance'][i]))
        
    return

