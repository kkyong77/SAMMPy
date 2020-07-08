
"""
Created on Thu Jun 25 20:26:55 2020

@author: Jing
"""
import numpy as np
import matplotlib.pyplot as plt

def plot(model, Ret):
    '''
    Create bar chart of variance-based sensitivity results

    Parameters
    ----------
    model: object, defined in sammpy
    Ret: dict, of sensitivity results

    '''
    plt.figure(figsize=(10, 4))
    x = np.arange(len(model.frames['names']))
    
    plt.subplot(1, 2, 1)
    plt.bar(x, Ret['PSK'])
    plt.xticks(x, model.frames['names'], fontsize=14)
    plt.ylabel('First-order process sensitivity')
        
        
    plt.subplot(1, 2, 2)
    plt.bar(x, Ret['PSTK'])
    plt.xticks(x, model.frames['names'], fontsize=14)
    plt.ylabel('Total-effect process sensitivity')
        
    plt.show()
    return 