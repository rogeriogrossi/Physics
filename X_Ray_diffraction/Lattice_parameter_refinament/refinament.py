"""
@author: Rogerio Murilo Grossi
"""

import lattice as lc

wl = 0.621485
l_parameter = 19
reflections = [
    [4,3,1],
    [6,3,1],
    [8,5,3],
    [12,3,1],
    [14,1,1]
]
peaks = [9.3323,12.4103,18.1539,22.8179,25.926]


a1 = lc.Cubic(reflections, peaks, l_parameter, wl)
l = lc.refine_lp(a1, n_int=6,var=1,save=True)
print(f'Lattice parameter best fit {l[0]:.5f}, variation of  {l[1]*100:.2f}%')
print(f'Final peaks positions: {", ".join([str(i) for i in a1.theta_2])}')


