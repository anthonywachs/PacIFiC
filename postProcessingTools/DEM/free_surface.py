import os
import sys
import numpy as np
import skimage
from skimage import *
from scipy import ndimage
import matplotlib.pyplot as plt
from cmdLineArguments import *

sys.path.append('./')

# Some variable that could change from one application to another
WIDTH = 20  # the width of the channel is 20cm
MIN_INTENSITY_DEFAULT = 0.05
pvpython_path = ('/home/damien/softwares/ParaView-5.8.0-RC1-MPI-Linux-Python3.7'
    +'-64bit/bin/pvpython')
generate_frames_path = '$PP_HOME/DEM/generate_first_last_frames.py'

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

def find_L0_from_file(rootfolderpath):
    insert_path=rootfolderpath+'/Grains/Init/insert.xml'
    insert_file = open(insert_path,'r')
    gate = False
    L0=""
    for line in insert_file:
        if gate == True and "Centre" in line:
            c=22
            while line[c]!='\"':
                L0+=line[c]
                c+=1
            L0=float(L0)*100-0.25
        if gate == False and "Gate" in line:
            gate = True
    return(L0)

def compute_vertical_intensity(frame):
    nbr, nbc = np.shape(frame)[:2]
    intensity = np.zeros(nbc)
    c0 = 0
    cf = 0
    for col in range(0,nbc):
        for row in range(0,nbr):
            intensity[col] += frame[row,col,0]
        if c0 == 0:
            if intensity[col] != 0:
                c0 = col
        elif cf == 0:
            if intensity[col] == 0:
                cf = col
    intensity /= (255*nbr)
    return intensity, c0, cf

def compute_horizontal_intensity(frame):
    nbr, nbc = np.shape(frame)[:2]
    intensity = np.zeros(nbr)
    r0 = 0
    rf = 0
    for row in range(0,nbr):
        for col in range(0,nbc):
            intensity[row] += frame[row,col,0]
        if r0 == 0:
            if intensity[row] != 0:
                r0 = row
        elif rf == 0:
            if intensity[row] == 0:
                rf = row
    intensity /= (255*nbr)
    return intensity, r0, rf


def set_origin_axes(first_frame):
    first_frame = filters.roberts(first_frame[:,:,0])
    nbr,nbc=np.shape(first_frame)
    negative_x=True
    c=0
    while negative_x and c<nbc:
        c+=1
        y_origin=np.argmax(first_frame[:,c])
        if y_origin!=0:
            negative_x=False
    y_origin=np.argmax(first_frame[:,c+10])
    x_origin=np.argmax(first_frame[y_origin-1,:])+1
    return(x_origin,y_origin)

def set_unit_length(first_frame,x0,y0,length,direction="side"):
    first_frame = first_frame[:y0,x0:,0]
    nbr,nbc=np.shape(first_frame)
    if direction == "side":
        r=nbr-10
        c=1
        while first_frame[r,c]>0:
            c+=1
        pixel_length=length/c
        return(pixel_length,c)
    elif direction == "top":
        r=nbr-1
        c=10
        while first_frame[r,c]>0:
            r-=1
        pixel_length=length/(nbr-r)
        return(pixel_length,nbr-r)

def analyze_first_frame(first_frame_side_path,L0):
    my_first_frame=io.imread(first_frame_side_path)
    (x0,y0)=set_origin_axes(my_first_frame)
    pixel_length,nb_white_pixels_t0=set_unit_length(my_first_frame,x0,y0,L0)

    my_initial_first_frame = my_first_frame
    my_first_frame = filters.roberts(my_first_frame[:y0,x0:,0])

    nbr = np.shape(my_first_frame)[0]
    x_t0 = []
    y_t0 = []
    for c in range(3,nb_white_pixels_t0-3):
        max=np.argmax(my_first_frame[:,c])
        if max<2:
            max=nbr
        x_t0.append(c)
        y_t0.append(max)

    H0=(y0-np.mean(y_t0))*pixel_length
    return(x0, y0, pixel_length, H0, my_initial_first_frame, x_t0, y_t0)

def analyze_first_frame_top(first_frame_top_path):
    my_first_frame=io.imread(first_frame_top_path)
    v_intensity,c0,cf = compute_vertical_intensity(my_first_frame)
    h_intensity,r0,rf = compute_horizontal_intensity(my_first_frame)
    pixel_length_top = width/(rf-r0)
    L0_top = (cf-c0)*pixel_length_top
    return(c0,cf,r0,rf,pixel_length_top,L0_top)

def analyze_last_frame_top(last_frame_top_path, pixel_length, c0, r0, rf,
    rootfolderpath):
    my_last_frame=io.imread(last_frame_top_path)[r0:rf,c0:,0]
    nbr,nbc = np.shape(my_last_frame)
    intensity = np.zeros(nbc)
    for col in range(0,nbc):
        for row in range(0,nbr):
            intensity[col] += my_last_frame[row,col]
        intensity[col] /= (nbr*255)
    col = 0
    L_infinity = 0
    while L_infinity == 0 and col<nbc:
        if intensity[col]<runout_intensity:
            L_infinity = col*pixel_length
        col += 1

    f, (ax1) = plt.subplots(1,1)
    x = np.linspace(0,nbc*pixel_length,nbc)
    ax1.imshow(my_last_frame,cmap='gray',extent=[0,nbc*pixel_length,0,width])
    ax1.plot(x, intensity*width,color='red')
    ax1.plot([L_infinity,L_infinity], [0,width],color='red',linestyle="--")
    ax1.set_xlabel("x (cm)")
    ax1.set_ylabel("y (cm)")
    plt.suptitle(rootfolderpath)
    return L_infinity


def analyze_last_frame(last_frame_side_path,pixel_length,x0,y0):
    image = io.imread(last_frame_side_path)
    initial_image = image
    nbrow,nbcol,nbch=np.shape(image)

    image = image[:y0,x0:,:]
    image = filters.roberts(image[:,:,0])

    nbr,nbc=np.shape(image)
    x=[]
    y=[]
    max_of_max=nbr
    runout=0
    for c in range(0,nbc):
        max=np.argmax(image[:,c])
        if max<2:
            max=nbr
        if max<max_of_max:
            max_of_max=max
        if c>0 and runout==0 and max==nbr:
            runout = c+x0
        x.append(c)
        y.append(max)

    if runout==0:
        print("Error: image does not show full extent of the avalanch")
    else:
        L_infinity = (runout-x0)*pixel_length

    y = smooth(y, 20)

    #correct the off-shoots due to the smoothing
    c=nbc-1
    while y[c]!=nbr and c>=0:
        y[c]=nbr
        c-=1
    c=0
    while y[c]<max_of_max:
        y[c]=max_of_max
        c+=1

    return(x,y,L_infinity,initial_image,runout)

def plot_free_surfaces(initial_image,my_initial_first_frame,x_t0,y_t0,x,y,
                        folderpath,pixel_length,runout):
    f, (ax1,ax2) = plt.subplots(1,2)
    nbr,nbc = np.shape(my_initial_first_frame)[:2]
    xl0 = -x0*pixel_length
    xr0 = nbc*pixel_length + xl0
    yl0 = -(nbr-y0)*pixel_length
    yr0 = nbr*pixel_length + yl0
    extent_t0 = [xl0,xr0,yl0,yr0]
    ax1.imshow(my_initial_first_frame,cmap='gray',extent=extent_t0)
    ax1.plot((np.array(x_t0))*pixel_length,np.ones(len(x_t0))*H0,color='red')
    ax1.set_xlabel("x (cm)")
    ax1.set_ylabel("y (cm)")
    #
    ax2.imshow(initial_image,cmap='gray',extent=extent_t0)
    ax2.plot(np.array(x[:])*pixel_length,yr0
        -np.array(y[:])*pixel_length,color='red')
    ax2.set_xlabel("x (cm)")
    ax2.set_ylabel("y (cm)")
    f.suptitle(folderpath)
    return()


### Start of the proper program
my_args = cmdLineArgs()
my_args.add_cmd_arg("run_paraview","run_pv")
my_args.add_cmd_arg("width","width")
my_args.add_cmd_arg("runout_intensity","runout_intensity")
my_args.read_cmd_args()
run_pv = True if my_args.get_attribute("run_paraview")=="yes" else False
width = float(my_args.get_attribute("width")) if my_args.get_attribute("width")!=None else WIDTH
runout_intensity =  float(my_args.get_attribute("runout_intensity")) if my_args.get_attribute("runout_intensity")!=None else MIN_INTENSITY_DEFAULT

output_file_path1 = "/home/damien/phd/dem/physics/fig3.10/damien.csv"
output_file1 = open(output_file_path1,"w")

output_file_path2 = "/home/damien/phd/dem/physics/fig3.11/damien_Hinf.csv"
output_file2 = open(output_file_path2,"w")

print("Parameters of this application")
print("Cannel width:",width,"cm")
print("Intensity threshold for run-out length:",100*runout_intensity,"%")
print("\n---------------------------------------------------------------\n")

aspect_ratios = []
runouts = []
free_surfaces = []
for rootfolderpath in my_args.get_attributes("root_path"):
    print("Entering "+rootfolderpath+"...")
    if rootfolderpath[-1]=='/':
        rootfolderpath.strip('/')
    folderpath = rootfolderpath + '/Grains/Simu/'

    L0 = find_L0_from_file(rootfolderpath)

    cmd = pvpython_path + ' ' + generate_frames_path + ' ' + folderpath

    if run_pv:
        os.system(cmd)

    first_frame_side_path=(folderpath+'initial_frame_side.png')
    last_frame_side_path=(folderpath+'last_frame_side.png')
    first_frame_top_path=(folderpath+'initial_frame_top.png')
    last_frame_top_path=(folderpath+'last_frame_top.png')

    c0_top, cf_top, r0_top, rf_top, pixel_length_top, L0_top =\
    analyze_first_frame_top(first_frame_top_path)

    L_infinity = analyze_last_frame_top(last_frame_top_path,
        pixel_length_top,c0_top,r0_top,rf_top,rootfolderpath)

    print("L_infinity from top frame =", L_infinity)

    print("L0 from file:",L0)
    print("L0 from top frame:",L0_top)

    (x0,y0,pixel_length,H0,my_initial_first_frame,
        x_t0,y_t0) = analyze_first_frame(first_frame_side_path,L0)
    (x,y,L_infinity_wrong,initial_image,runout) = analyze_last_frame(
                                last_frame_side_path, pixel_length,x0,y0)

    print("H0 =",H0)
    print("L0 =",L0)
    print("L_infinity =", L_infinity)
    print("L_infinity from the side (over estimated) =", L_infinity_wrong)

    plot_free_surfaces(initial_image, my_initial_first_frame,x_t0,y_t0,x,y,
    rootfolderpath,pixel_length,runout)

    aspect_ratios.append(H0/L0)
    runouts.append(L_infinity)

    # Often some crosses are leaning against the wall, falsely increasing the
    # final height of the pile. We consider the final height 3 times the
    # particle width away from the wall, since the image is blurry and the
    # y function is smoothed out: we need this margin.
    nb_pixels_per_particle_width = round(3 * 0.3/pixel_length)
    print("nb_pixels_per_particle_width =",nb_pixels_per_particle_width)
    x=np.array(x[:runout])
    y=y0-np.array(y[:runout])
    y_hinf=y[nb_pixels_per_particle_width:runout]
    H_infinity = np.max(y_hinf)*pixel_length
    print("H_infinity =",H_infinity)

    aspect_ratio = H0 / L0
    nd_runout = (L_infinity - L0) / L0
    nd_height = H0 / H_infinity

    print("aspect ratio a="+str(aspect_ratio)+", (L_inf-L0)/L0="+str(nd_runout)
        +", H_inf/H0="+str(H0/H_infinity))
    output_file1.write(str(aspect_ratio)+", "+str(nd_runout)+"\n")
    output_file2.write(str(aspect_ratio)+", "+str(nd_height)+"\n")


    # free_surfaces.append(((x*pixel_length-L0)/(L_infinity-L0),y/np.max(y),
    #                 rootfolderpath))
    free_surfaces.append(((x*pixel_length)/(H_infinity),(y*pixel_length)/H_infinity,rootfolderpath))
    print("\n---------------------------------------------------------------\n")

output_file1.close()
output_file2.close()

fig = plt.figure()

for case in free_surfaces:
    color = my_args.get_attribute("color",case[2])
    linestyle = my_args.get_attribute("linestyle",case[2])
    linewidth = my_args.get_attribute("linewidth",case[2])
    label = my_args.get_attribute("label",case[2])
    plt.plot(case[0],case[1],label=label,c=color,ls=linestyle,lw=linewidth)
# plt.xlabel(r'$\frac{L-L_0}{L_\infty-L_0}$')
# plt.ylabel(r'$\frac{h}{H_\infty}$')
plt.xlabel(r'$\frac{x}{H_\infty}$')
plt.ylabel(r'$\frac{h}{H_\infty}$')
if my_args.is_legend():
    plt.legend(loc='best')
plt.show()
