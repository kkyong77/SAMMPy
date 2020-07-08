# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 14:54:40 2020

@author: Jing
"""
import math
import numpy as np
import matplotlib.pyplot as plt


def plot(model, values):
    '''
    Create histogram plots for sample realizations

    Parameters
    ----------
    model: object, defined in sammpy
    Ret: dict, of sensitivity results

    '''
    num_pars = values.shape[1]
    N = values.shape[0]
    
    # Figure configuration
    nrows = 2
    ncols = math.ceil(num_pars / 2)
    width = ncols * 5
    height = 7
    bins = 25
    
    # Plot the pdfs
    plt.figure(figsize=(width, height))

    for i in range(num_pars):
        plt.subplot(nrows, ncols, i + 1)
        plt.hist(values[:, i], density=True, bins=bins)
        plt.xlabel('$' + model.pars['names'][i] + '$', fontsize=14)
        plt.ylabel('$PDF$', fontsize=14)

    plt.show()
    return 