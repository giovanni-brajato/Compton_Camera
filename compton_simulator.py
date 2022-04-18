import math
import numpy as np
import matplotlib.pyplot as plt
import plotly

# define the functions to have a compact notation
from numpy import pi, cos, sin, tan, sqrt, arccos, arctan, sign, rad2deg, multiply, transpose, newaxis
rng = np.random.default_rng()
# define the parameters

# CONSTANTS
mec2 = 0.511e6  # Energy of electron mass at rest

N = int(1e2) # number of points

# energies with uncertainties
sigma_E1 = 5e2
sigma_E2 = 5e2
mean_E1 = 350e3
#mean_E2 = 200e3
mean_E2 = (- mean_E1 + sqrt(mean_E1**2 + 4/0.2*mean_E1*mec2))/2
E1 = mean_E1 + rng.normal(loc=0,scale=1,size=(N)) * sigma_E1  # in EV
E2 = mean_E2 + rng.normal(loc=0,scale=1,size=(N)) * sigma_E2 # in EV

# locations with uncertainties
sigma_xyz1 = np.identity(3)*1e-6
sigma_xyz2 = np.identity(3)*1e-7
mean_xyz1 = (1e-2, 3e-2, 13e-2 )
mean_xyz2 = (2e-2, 0, 10e-2)
XYZ1 = rng.multivariate_normal(mean=mean_xyz1, cov=(sigma_xyz1), size=N)
X1, Y1, Z1 = XYZ1[:,0], XYZ1[:,1], XYZ1[:,2]
XYZ2 = rng.multivariate_normal(mean=mean_xyz2, cov=(sigma_xyz2), size=N)
X2, Y2, Z2 = XYZ2[:,0], XYZ2[:,1], XYZ2[:,2]


# camera plane location
z = 30e-2

# calculate the angles
length = sqrt((X1 - X2) ** 2 + (Y1 - Y2) ** 2)
theta = arctan((Z1 - Z2) / length)  # the co-latitude, in radiant
phi = -pi/2*(sign(Y1-Y2 + sign(X1-X2)))+sign(Y1-Y2)*sign(X1-X2)*(pi/2*sign(Y1-Y2)-arctan((X1-X2)/(Y1-Y2)))  # the
# longitude, in radiant
theta_deg = np.rad2deg(theta)
phi_deg = np.rad2deg(phi)


# Calculate angle compton alpha

alpha = arccos(1 - mec2*E1 / ((E1 + E2) * E2))
alpha_deg = np.rad2deg(alpha)
# print("Compton angle is : ", alpha, " rad, or", alpha / (pi * 2) * 360, "°")

# Calculate compton cone/ellipse on the camera plane
a = cos(theta)**2 - (sin(theta) * tan(alpha))**2
b = z * (2 + cos(theta) * sin(theta) + 2 * cos(theta) * sin(theta) * tan(alpha)**2)
c = z**2 * (sin(theta)**2 - (cos(theta) * tan(alpha))**2)
d = cos(phi)
f = sin(phi)

A = d**2 + a * f**2
B = d * f * (a - 1)
C = f**2 + a * d**2
D = -b * f / 2
F = -b * d / 2
G = c

# Drawing the ellipse
x0 = (C * D - B * F) / (B**2 - A * C) + X1  # x-position of the center
y0 = (A * F - B * D) / (B**2 - A * C) + Y1  # y-position of the center

a1 = sqrt((2 * (A * F**2 + C * D**2 + G * B**2 - 2 * B * D * F - A * C * G)) / (
            (B**2 - A * C) * (sqrt((A - C)**2 + 4 * B**2) - A - C)))  # radius on the x-axis
b1 = sqrt((2 * (A * F**2 + C * D**2 + G * B**2 - 2 * B * D * F - A * C * G)) / (
            (B**2 - A * C) * (- sqrt((A - C)**2 + 4 * B**2) - A - C)))  # radius on the y-axis
mu = sign(A - C) * (pi / 4 + sign(A - C) * (pi / 4 + 0.5 * arctan(2 * B / (A - C))))  # rotation angle


N_points_ellipse = int(1e2)
eta = (np.linspace(0, 2 * pi,N_points_ellipse))[np.newaxis]
Ell = np.array([multiply(transpose(a1[newaxis]),cos(eta)) ,multiply(transpose(b1[newaxis]),sin(eta))])  # x0,y0 removed to keep the same center location
R_rot = np.array([[cos(mu), -sin(mu)], [sin(mu), cos(mu)]])  # 2-D rotation matrix

Ell_rot = np.zeros(Ell.shape)

for n in range(N):
    for i in range(N_points_ellipse):
        Ell_rot[:,n, i] = np.dot(R_rot[:,:,n], Ell[:,n, i])

plt.figure()
for n in range(N):
    plt.plot(x0[n] + Ell_rot[0, n,:], y0[n] + Ell_rot[1,n, :], color=np.random.rand(3,))  # rotated ellipse


# naming the x axis
plt.xlabel('x - axis')
# naming the y axis
plt.ylabel('y - axis')

# giving a title to my graph
plt.title('Compton conic sections')


plt.figure()
plt.hist(alpha_deg, bins='auto')  # arguments are passed to np.histogram
plt.title(chr(945))
# function to show the plot
plt.show()


