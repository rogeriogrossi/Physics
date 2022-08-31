"""
Finds the best fit for the lattice parameter of a XPD experiment.

By means of Laue equation and Bragg's law this module fits the lattice parameter
of a X ray powder diffraction pattern.
For theorical information visit:
https://handwiki.org/wiki/Physics:Laue_equations

This module requires the following information:
    
    The radiation wavelength (wl) in angstroms.
    A first guess for the lattice parameter. One must use angstroms for unit.
    A list of the hkl reflections.
    A list of 2-theta positions which are observed information of the X ray 
    pattern.


@author: Rogerio Murilo Grossi
"""

import numpy as np

class Cubic():
    
    def __init__(self,reflections, peaks, l_parameter, wl):
        self.reflections = reflections
        self.peaks = peaks
        self.l_parameter = l_parameter
        self.wl = wl

     
    def f_theta(self):
        """
        Calculates the theorical 2theta position
        
        Returns: A list of 2theta angle

        """
        self.theta_2=[]
        for rfl in self.reflections:
             #plane_dist: Find the interplanar spacing between two planes.
             plane_spacing = np.sqrt(sum([i**2 for i in rfl]))
             #angle: find the 2theta position by using the bragg's law
             angle = np.round(2*np.arcsin(self.wl*plane_spacing/(self.l_parameter*2))*180/np.pi,5)
             self.theta_2.append(angle)
        return self.theta_2   
    
    def diff_per_obs(self):
        """
        Calculates the norm of the difference between theoretical and observed peaks
        
        Returns: The difference per observation

        """
        self.diff = []
        self.theta = self.f_theta()
        
        for i, pec in enumerate(self.peaks):
            self.diff.append(np.abs(pec - self.theta[i]))  
            #self.diff_norm = sum([np.abs(i) for i in self.diff])
            
        self.diff_final = sum(self.diff)/len(self.diff)
        return np.round(self.diff_final,5)

def refine_lp(obj,n_int=3, var = 1, re = 0.1,save=False):
    """
    Lattice parameter fit 
    
    
    Parameters
    ----------
    obj : Lattice object
    n_int : Number of interections
        type: int
    var : Variation of the lattice parameter 
        type: float
    re : Resolution
        type float
    save: Save the results in a .txt file
        type: bool
   
    Returns: The lattice parameter fit
    

    """
    a = obj.l_parameter
    
    for cicle in range(n_int):
        lp_values = [val for val in np.arange(-var,var,re)]
        dif_a = []
        for i in lp_values:
            obj.l_parameter = a+i
            dif_a.append(obj.diff_per_obs()) 
        # Parametro de rede se desloca para o mínimo
        obj.l_parameter = a + lp_values[dif_a.index(min(dif_a))]
        var /= 5
        re /= 5  
        a += lp_values[dif_a.index(min(dif_a))]
        
    if save==True:    
        theta_2 = [str(i) for i in obj.theta_2]
        reflections = [str(i) for i in obj.reflections]
        
        with open('refinament_Results.txt','a+') as file:
            file.write('-'*52)
            file.write(f'\n Lattice parameter: {a:.5f}\n Lambda: {obj.wl:.5f}\n ')
            file.write(f'Variation per peak: {min(dif_a)*100:.3f}%\n\n')
            file.write(f'Peaks: {"; ".join(theta_2)}\nReflections: {"; ".join(reflections)}\n\n')
            file.write('=-'*12 + 'End' + '=-'*12+'\n\n')
    return (obj.l_parameter,min(dif_a))
       

        
            