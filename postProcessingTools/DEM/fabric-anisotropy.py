from math import *
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from cmdLineArguments import *
from celluloid import Camera

EPS = 1.e-6

### --- Set equal axis for 3D plots. Code provided by Karlo on StackOverflow:
# https://stackoverflow.com/questions/13685386/matplotlib-equal-unit-length-with-equal-aspect-ratio-z-axis-is-not-equal-to
def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
###

class Point:
    """In this class, we simply define a 3D point as an array."""
    def __init__(self, my_x, my_y, my_z):
        self.coord =[my_x, my_y, my_z]
    def __getitem__(self, i):
        return self.coord[i]
    def __setitem__(self, i, x):
        self.coord[i] = x
    def __str__(self):
        return ("(" + str(self.coord[0]) + ", " + str(self.coord[1]) +
        ", " + str(self.coord[2]) + ")")

class Contact:
    """This class gathers information about a specific contact for post-processing purposes. For now, it only has information about the centers of the two colliding particles (or the center of the particle and the point of contact in case of a collision with a wall).
    """
    def __init__(self, pt0, pt1):
        self.pts = [pt0, pt1]

    def get_length(self):
        length = 0.
        for d in range(3):
            length += (self.pts[0][d] - self.pts[1][d])**2
        return sqrt(length)

    def get_xz_angle(self):
        """This function returns the angle of the contact line projected onto the xz plane (perpendicular to the ground and to the gate).
        """
        if (abs(self.pts[0][0] - self.pts[1][0]) < EPS):
            return pi/2
        else:
            angle = atan((self.pts[0][2] - self.pts[1][2])/
            (self.pts[0][0] - self.pts[1][0]))
            if angle < 0:
                return angle + pi
            return angle

    def get_azimuthal(self):
        """This function returns the azimuthal angle formed by the contact line in a spherical coordinates system
        """
        if (abs(self.pts[0][0] - self.pts[1][0]) < EPS):
            return 0.
        else:
            return atan((self.pts[0][1] - self.pts[1][1])/(self.pts[0][0] - self.pts[1][0]))

    def draw(self, ax, my_color = "black"):
        """Assuming a 3D figure has been initialized and a plt.show() will be called later, this function draws the line of this contact.
        """
        ax.plot([self.pts[i][0] for i in range(2)],
        [self.pts[i][1] for i in range(2)],[self.pts[i][2] for i in range(2)],
        color = my_color)

    def is_admissible(self, zmin):
        admissible = True
        for i in range(2):
            if (self.pts[i][2] <= zmin or self.pts[i][0]<=0.015+EPS or
                self.pts[i][0]>=0.145-EPS):
            # if (self.pts[i][2] <= zmin or self.pts[i][0]<=0.0015 or
            #     self.pts[i][0]>=0.0515):
                admissible = False
        return admissible

class ContactNetwork:
    """In this class we store all the contact lines and define some functions useful to the analysis of the contact network.

    Attributes:
        allContacts: a simple list of Contact objects
        nbContacts: an int storing the number of contacts *including the contacts with the walls*
    """
    def __init__(self):
        self.allContacts = []
        self.nbContacts = 0

    def add_contact(self, Contact):
        self.allContacts.append(Contact)

    def rm_contact(self, index):
        self.allContacts.pop(index)

    def read_vtp(self, filePath, nb_procs = 1, no_obs = False, zmin = -1.e10):
        """This function reads the VTP files storing contact informations written by Grains3D.

        For now, this function only reads line 7 of these files, which stores the coordinates of the centers of the colliding particles (or the contact point if the contact occurs with a wall).

        The data are simply a long line of floats separated by spaces, and they correspond to x, y, z coordinates of each contact points.
        """
        coord = np.array([])
        nb_obs = 0
        for i in range(nb_procs):
            current_file_path = filePath[:-5]+str(i)+filePath[-4:]
            try:
                myFile = open(current_file_path, "r")
            except IOError:
                print("Could not open file " + current_file_path)

            line_nb = 0
            for line in myFile:
                if line_nb == 3:
                    if no_obs:
                        nb_obs = int(line.split('=')[-1].split('\"')[1])
                    else:
                        nb_obs = 0
                if line_nb == 6:
                    if no_obs:
                        coord_i = np.array((line.strip(' \n')).split(" "))[:-(nb_obs)*6]
                    else:
                        coord_i = np.array((line.strip(' \n')).split(" "))
                    break
                line_nb += 1
            coord = np.hstack((coord,coord_i))
        coord = coord.astype(np.float)
        for i in range(0,int(len(coord)/6)):
            current_contact = Contact(Point(coord[6*i], coord[6*i+1],
            coord[6*i+2]), Point(coord[6*i+3], coord[6*i+4], coord[6*i+5]))
            if current_contact.is_admissible(zmin):
                self.add_contact(current_contact)
        self.nbContacts = len(self.allContacts)

    def draw_network(self):
        """This function builds on Contact.draw() to output the contact network in a 3D plot. We use the function set_axes_equal provided by Karlo (see at the very top of this file).
        """
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        for i in range(self.nbContacts):
            self.allContacts[i].draw(ax)
        set_axes_equal(ax)
        plt.show()

    def fabric_anisotropy(self, nbSamples = 51, time = 0):
        """This function computes and draw the probability density function P_n(theta) = Nc(theta)/Nc, where Nc(theta) is the number of contact which orientation is between theta and theta+dtheta, and Nc is the total number of contacts. See, e.g., D. Cantor et al., "Rheology and structure of polydisperse three-dimensional packings of spheres", Phys. Rev. E., 2018.

        This function allows to tune through the argument nbSamples, since dtheta = pi/nbSamples.
        """
        nbTheta = np.zeros((1,nbSamples))
        for i in range(self.nbContacts):
            theta = self.allContacts[i].get_xz_angle()
            if abs(theta - pi/2) < EPS:
                print("theta =", theta, "pt0 =",
                self.allContacts[i].pts[0], ", pt1 =",
                self.allContacts[i].pts[1])
            nbTheta[0,int(nbSamples*theta/pi)] += 1
        nbTheta /= self.nbContacts
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(np.array([(i*pi/nbSamples)%(2*nbSamples) for i in range(int(2*nbSamples+1))]), np.transpose(np.hstack((nbTheta, nbTheta, np.array([[nbTheta[0,0]]])))))
        plt.title("Fabric anisotropy, crosses, t="+str(time))
        ax.set_yticklabels([])
        plt.show()

# Example:
# Command to type in $PACIFIC_HOME/postProcessingTools:
# python DEM/fabric-anisotropy.py -file=./test-cases/crosses_10_procs/simul_ForceChain_T1_0.vtp -nb-procs=10
my_args = cmdLineArgs()
my_args.add_cmd_arg("file","file")
my_args.add_cmd_arg("nb-outputs","nb-outputs")
my_args.add_cmd_arg("nb-procs","nb-procs")
my_args.add_cmd_arg("nb-samples","nb-samples")
my_args.add_cmd_arg("draw-network","draw-network")
my_args.add_cmd_arg("no-obstacle","no-obstacle")
my_args.read_cmd_args()
file_path = my_args.get_attribute("file")
nb_out = my_args.get_attribute("nb-outputs")
nb_out = int(nb_out) if (nb_out != None) else 1
nb_procs = my_args.get_attribute("nb-procs")
nb_procs = int(nb_procs) if (nb_procs != None) else 1
nbSamples = my_args.get_attribute("nb-samples")
nbSamples = int(nbSamples) if (nbSamples != None) else 51
no_obs = my_args.get_attribute("no-obstacle")
no_obs = True if (no_obs != None) else False

draw_network = True if my_args.get_attribute("draw-network") != None else False

radius = 0.0015
for i_file in range(nb_out):
    myContacts = ContactNetwork()
    myContacts.read_vtp(file_path, nb_procs, no_obs, zmin=radius+EPS)
    print("Nb of considered contacts =", myContacts.nbContacts)
    myContacts.fabric_anisotropy(nbSamples = nbSamples)

if draw_network:
    myContacts.draw_network()
