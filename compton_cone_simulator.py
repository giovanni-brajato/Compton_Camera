import math
import numpy as np
import matplotlib.pyplot as plt

# define the functions to have a compact notation
from math import pi, cos, sin, tan, pow, sqrt, acos, atan

def sign(x):
    if x > 0: return 1
    elif x < 0: return -1
    else: return 0


# define the parameters

# CONSTANTS
mec2 = 0.511e6  # Energy of electron mass at rest

# angles
theta = pi / 4  # the co-latitude, in radiant
phi = pi / 5  # the longitude, in radiant

# energies
E1 = 1e6  # in EV
E2 = 1e6  # in EV

# locations
x1, y1, z1 = 0, 0, 0
x2, y2, z2 = 0, 0, -10

# camera plane location
z = 10

# Calculate angle compton alpha

alpha = acos(1 - E1 * mec2 / (E1 + E2) / E2)
print("Compton angle is : ", alpha, " rad, or", alpha / (pi * 2) * 360, "°")

# Calculate compton cone
a = pow(cos(theta), 2) - math.pow(sin(theta) * tan(alpha), 2)
b = z * (2 + cos(theta) * sin(theta) + 2 * cos(theta) * sin(theta) * pow(tan(alpha), 2))
c = pow(z, 2) * (pow(sin(theta), 2) - pow(cos(theta) * tan(alpha), 2))
d = cos(phi)
f = sin(phi)

A = pow(d, 2) + a * pow(f, 2)
B = d * f * (a - 1)
C = pow(f, 2) + a * pow(d, 2)
D = -b * f / 2
F = -b * d / 2
G = c

# Drawing the ellipse
x0 = (C * D - B * F) / (pow(B, 2) - A * C) + x1  # x-position of the center
y0 = (A * F - B * D) / (pow(B, 2) - A * C) + y1  # y-position of the center
a1 = sqrt((2*(A*pow(F,2) + C*pow(D,2) + G*pow(B,2) - 2*B*D*F - A*C*G))/((pow(B,2) - A*C)*(sqrt(pow(A-C,2) + 4*pow(B,2)) - A - C))) # radius on the x-axis
b1 = sqrt((2*(A*pow(F,2) + C*pow(D,2) + G*pow(B,2) - 2*B*D*F - A*C*G))/((pow(B,2) - A*C)*(-sqrt(pow(A-C,2) + 4*pow(B,2)) - A - C))) # radius on the y-axis
mu = sign(A-C)*(pi/4 + sign(A-C)*(pi/4 + 0.5*atan(2*B/(A-C))))  # rotation angle

eta = np.linspace(0, 2 * pi, 100)
Ell = np.array([a1 * np.cos(eta), b1 * np.sin(eta)]) # x0,y0 removed to keep the same center location
R_rot = np.array([[cos(mu), -sin(mu)], [sin(mu), cos(mu)]]) # 2-D rotation matrix


Ell_rot = np.zeros((2, Ell.shape[1]))
for i in range(Ell.shape[1]):
    Ell_rot[:, i] = np.dot(R_rot, Ell[:, i])

plt.plot(x0 + Ell[0, :], y0 + Ell[1, :] , label ='Initial Ellipse')  # initial ellipse
plt.plot(x0 + Ell_rot[0, :], y0 + Ell_rot[1, :], 'darkorange',  label ='Rotated Ellipse')  # rotated ellipse
plt.grid(color='lightgray', linestyle='--')

# naming the x axis
plt.xlabel('x - axis')
# naming the y axis
plt.ylabel('y - axis')

# giving a title to my graph
plt.title('Compton conic sections')
plt.legend()

# function to show the plot
plt.show()
