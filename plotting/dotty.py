"""
Created on Thu Jun 25 21:52:56 2020

@author: Jing
"""

import numpy as np
import matplotlib.pyplot as plt

def plot(model, Ret):
    '''
    Carte a scatter plot display the mean against variance of the output 
    difference

    Parameters
    ----------
    model: object, defined in sammpy
    Ret: dict, of sensitivity results

    '''
    plt.figure(figsize=(5, 5))

    plt.scatter(Ret['mean'], Ret['variance'], s=40)
    for i in range(len(model.frames['names'])):
        plt.text(Ret['mean'][i], Ret['variance'][i], model.frames['names'][i],
                 fontsize=14)
    plt.xlabel('Mean', fontsize=14)
    plt.ylabel('Variance', fontsize=14)
    
    plt.show()
    
    return 



