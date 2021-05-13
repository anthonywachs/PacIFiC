from math import *
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
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
            return 0.
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

    def read_vtp(self, filePath):
        """This function reads the VTP files storing contact informations written by Grains3D.

        For now, this function only reads line 7 of these files, which stores the coordinates of the centers of the colliding particles (or the contact point if the contact occurs with a wall).

        The data are simply a long line of floats separated by spaces, and they correspond to x, y, z coordinates of each contact points.
        """
        myFile = open(filePath,'r')
        line_nb = 0
        for line in myFile:
            if line_nb == 6:
                coord = np.array((line.strip(' \n')).split(" "))
                break
            line_nb += 1
        coord = coord.astype(np.float)
        self.nbContacts = int(len(coord)/6)
        for i in range(0,self.nbContacts):
            self.add_contact(Contact(Point(coord[6*i], coord[6*i+1],
            coord[6*i+2]), Point(coord[6*i+3], coord[6*i+4], coord[6*i+5])))

    def draw_network(self):
        """This function builds on Contact.draw() to output the contact network in a 3D plot. We use the function set_axes_equal provided by Karlo (see at the very top of this file).
        """
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        for i in range(self.nbContacts):
            self.allContacts[i].draw(ax)
        set_axes_equal(ax)
        plt.show()

    def fabric_anisotropy(self, nbSamples = 51):
        """This function computes and draw the probability density function P_n(theta) = Nc(theta)/Nc, where Nc(theta) is the number of contact which orientation is between theta and theta+dtheta, and Nc is the total number of contacts. See, e.g., D. Cantor et al., "Rheology and structure of polydisperse three-dimensional packings of spheres", Phys. Rev. E., 2018.

        This function allows to tune through the argument nbSamples, since dtheta = pi/nbSamples.
        """
        nbTheta = np.zeros((1,nbSamples))
        nbRelevantContacts = 0
        for i in range(self.nbContacts):
            theta = self.allContacts[i].get_xz_angle()
            # Because the contacts with walls add huge spikes in the horizontal
            # and vertical directions, we try to sort them out with this if
            # statement. The idea is that if a contact line is exactly
            # horizontal or vertical, it is very likely that it is a contact
            # with a wall.
            # Interestingly, the spikes in the normal direction are very
            # sensitive to the value of EPS.
            if (abs(theta) > EPS and abs(theta - pi) > EPS and
            abs(theta - pi/2) > EPS):
                nbTheta[0,int(nbSamples*theta/pi)] += 1
                nbRelevantContacts += 1
        nbTheta /= nbRelevantContacts
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(np.array([(i*pi/nbSamples)%(2*nbSamples) for i in range(int(2*nbSamples+1))]), np.transpose(np.hstack((nbTheta, nbTheta, np.array([[nbTheta[0,0]]])))))
        plt.show()

# Example
file_path = "test-cases/fabric_anisotropy_test2.vtp"
myContacts = ContactNetwork()
myContacts.read_vtp(file_path)
myContacts.fabric_anisotropy()
myContacts.draw_network()
