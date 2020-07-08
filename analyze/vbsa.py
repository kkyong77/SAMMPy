# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 10:31:02 2020

@author: Jing
"""
import numpy as np
import numba as nb
import itertools
from ..util import results

def analyze(model, Y, print_to_console=True):
    """
    Perform variance-based sensitivty analysis (VBSA) for each process.

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
    
    # Perfrom the first-order process sensitivty anlaysis
    if print_to_console:
        print('Runing VBSA first-order process sensitivy analysis...') 
    PSI = first_order_PSK(model, Y)
    
    # Perfrom the total*effect process sensitivty anlaysis
    if print_to_console:
        print('Runing VBSA total-effect process sensitivy analysis...') 
    PST = total_effect_PSTK(model, Y)
    
    # Save results to the dict
    for i in range(npros):
            S['PSK'][i] = PSI[i]
            S['PSTK'][i] = PST[i]

    # Print results to console
    if print_to_console:
        print_indices(model, S)
        
    return S


def first_order_PSK(model, Y):
    '''
    Perform first-order process sensitivty anlaysis 

    Returns
    ========
    numpy.array : A 1D numpy array storing the first-order process sensitivy 
    index for each process.
    '''
    # Number of sample realizations
    N = Y.shape[1]
    
    # Number of process and process models
    npros = len(model.frames['names'])
    npros_models_list = [len(model.frames['options'][i]) for i in range(npros)]
    
    # Initilize a array to store the results
    PSI = np.zeros(npros)

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
        def numba_psi(n_K_models, n_non_K_models, N, PM_K, PM_non_K, tmp_Y):
            # compute the first-order process sensitivty index
            E_tb_Y = np.zeros((n_K_models, N, n_non_K_models))
            E2_tb_Y = np.zeros((n_K_models, N, n_non_K_models))
            E_b_Y  = np.zeros((n_K_models, N))   
            E_b_Y2  = np.zeros((n_K_models, N))   
            E2_b_Y  = np.zeros((n_K_models, N))
            E_ta_Y = np.zeros(n_K_models)
            E_ta_Y2 = np.zeros(n_K_models)
            E2_ta_Y = np.zeros(n_K_models)
                    
                               
            for i in range(n_K_models):
                for j in range(N):
                    for k in range(n_non_K_models):
                        E_tb_Y[i, j, k] = np.mean(tmp_Y[i, j, k, :])
                        E2_tb_Y[i, j, k] = np.mean(tmp_Y[i, j, k, :]**2)
                        
                    E_b_Y[i, j] = np.sum(E_tb_Y[i, j, :] * PM_non_K)
                    E_b_Y2[i, j] = E_b_Y[i, j]**2
                    E2_b_Y[i, j] = np.sum(E2_tb_Y[i, j, :] * PM_non_K)
                    
                E_ta_Y[i] = np.mean(E_b_Y[i, :])
                E_ta_Y2[i] = np.mean(E_b_Y2[i, :])
                E2_ta_Y[i] = np.mean(E2_b_Y[i, :])
                
            E_a_Y = np.sum(E_ta_Y[:] * PM_K)
            E_a_Y2 = np.sum(E_ta_Y2[:] * PM_K)
            E2_a_Y = np.sum(E2_ta_Y[:] * PM_K)
            
            Var_T = E2_a_Y - E_a_Y**2
            Var_Y = E_a_Y2 - E_a_Y**2
            
            Si = Var_Y / Var_T
         
            return Si
        
        Si = numba_psi(n_K_models, n_non_K_models, N, PM_K, PM_non_K, tmp_Y)
        PSI[ipros] = Si
        
    return PSI


def total_effect_PSTK(model, Y):
    '''
    Perform total-effect process sensitivty anlaysis 

    Returns
    ========
    numpy.array : A 1D numpy array storing the total-effect process sensitivy 
    index for each process.
    '''
    # Number of sample realizations
    N = Y.shape[1]
    
    # Number of process and process models
    npros = len(model.frames['names'])
    npros_models_list = [len(model.frames['options'][i]) for i in range(npros)]
    
    # Initilize a array to store the results
    PST = np.zeros(npros)
    
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
        def numba_pst(n_K_models, n_non_K_models, N, npros, PM_K, PM_non_K, tmp_Y):
            # compute the first-order process sensitivty index
            E_tb_Y = np.zeros((n_non_K_models, N**(npros - 1), n_K_models))
            E2_tb_Y = np.zeros((n_non_K_models, N**(npros - 1), n_K_models))
            E_b_Y  = np.zeros((n_non_K_models, N**(npros - 1)))   
            E_b_Y2  = np.zeros((n_non_K_models, N**(npros - 1)))   
            E2_b_Y  = np.zeros((n_non_K_models, N**(npros - 1)))
            E_ta_Y = np.zeros(n_non_K_models)
            E_ta_Y2 = np.zeros(n_non_K_models)
            E2_ta_Y = np.zeros(n_non_K_models)
                                    
            for i in range(n_non_K_models):
                for j in range(N**(npros - 1)):
                    for k in range(n_K_models):
                        E_tb_Y[i, j, k] = np.mean(tmp_Y[k, :, i, j])
                        E2_tb_Y[i, j, k] = np.mean(tmp_Y[k, :, i, j]**2)
                        
                    E_b_Y[i, j] = np.sum(E_tb_Y[i, j, :] * PM_K)
                    E_b_Y2[i, j] = E_b_Y[i, j]**2
                    E2_b_Y[i, j] = np.sum(E2_tb_Y[i, j, :] * PM_K)
                    
                E_ta_Y[i] = np.mean(E_b_Y[i, :])
                E_ta_Y2[i] = np.mean(E_b_Y2[i, :])
                E2_ta_Y[i] = np.mean(E2_b_Y[i, :])
                
            E_a_Y = np.sum(E_ta_Y[:] * PM_non_K)
            E_a_Y2 = np.sum(E_ta_Y2[:] * PM_non_K)
            E2_a_Y = np.sum(E2_ta_Y[:] * PM_non_K)
            
            Var_T = E2_a_Y - E_a_Y**2
            Var_Y = E_a_Y2 - E_a_Y**2
            
            St = 1 - Var_Y / Var_T
            
            return St
    
        St = numba_pst(n_K_models, n_non_K_models, N, npros, PM_K, PM_non_K, tmp_Y)
        PST[ipros] = St
            
    return PST

def create_si_dict(npros):
    # initialize empty dict to store process sensitivity indices
    S = results.ResultDict((k, np.zeros(npros)) for k in ('PSK', 'PSTK'))
    return S
            
def print_indices(model, S):
    # Output to console
    title = 'Process'
    names = model.frames['names']
    npros = npros = len(model.frames['names'])
    print('%s \t PSK \t PSTK' % title)

    for i in range(npros):
        print('%s \t %.4f \t %.4f' % (names[i], S['PSK'][i], S['PSTK'][i]))
        
    return
        
        

