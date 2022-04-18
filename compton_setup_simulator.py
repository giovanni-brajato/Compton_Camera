# import the module
from vpython import *


# the compton cone
comtpon_cone = cone(pos=vector(1, 1, 0),
     length=3,
     radius=1,
     color=vector(1, 0, 0))
# comtpon_cone.rotate(angle=math.pi/6,
#            axis=vec(0,0,1),
#            origin=vector(0,0,0),color=vector(1, 1, 0))

# the diffuser
diffuser = box(pos=vector(0,0,0), axis=vec(1,0,0), length=4,
    height=4, width=1, opacity=1)


# the absorber

absorber = box(pos=vector(0,0,-4), axis=vec(0,0,0), length=4,
    height=4, width=2, opacity=1)