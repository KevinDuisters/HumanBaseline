#===============================================================#
# Python training
height = 1.90
weight = 75
bmi = weight/height**2
print(bmi)

# printing strings (paste). Use the + for strings
print('look a string ' + str(bmi))

# Lists of lists
fam = [ ["Kev", 1.90], ["John", 1.81] , ["Emily", 1.65]  ]

# indexing!! Start at 0; Ends at -1. A range is from up to (not included!)
fam[1]
fam[0:2]
fam[:2]
fam[1:]

# adding and removing elements
fam+[["Harry",0.4]] 
print(fam) # does not change variable
del(fam[1]) # this changes the variable!
fam

#===============================================================#
# Pointers
# list points to memory containing values, so copying must be done explicit
y = fam
y[1] = ["New", 0.0]
fam[1]

famcopy = fam[:]  # or list(fam)
famcopy[1] = ["Emily",1.65]
fam[1]

#===============================================================#
# Functions and Methods
# function help(sorted)
# methods (associated with type) help(str)

#===============================================================#
import numpy as np 
#import math as math
from math import pi
x=np.array([1,2,3])
print(pi)
help(np.array)
# shape gives attribute info (same as str in R)
x.shape

#===============================================================#
# Distributions
x = np.random.normal(0,1,100)
y = np.random.normal(3,2,100)
z = np.column_stack([x,y])

# Visualisation
import matplotlib.pyplot as mpl 
plt = mpl.plot(z)
mpl.show(plt)

