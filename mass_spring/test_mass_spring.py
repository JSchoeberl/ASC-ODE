import sys
sys.path.append('/Users/joachim/texjs/lva/ws2324/ScientificComputing/experiments/ASC-ODE/build/mass_spring')

from mass_spring import *


mss = MassSpringSystem2d()

mss.Add (Mass(1, (1,0)))


