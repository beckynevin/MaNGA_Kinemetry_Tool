#!/usr/bin/env python
# coding: utf-8

# In[170]:


# Step 1: Go to Marvin Github and import some tutorial stuff
from marvin.tools.maps import Maps
import matplotlib.pyplot as plt
import numpy as np
import marvin
from marvin import config, marvindb
import math
from bresenham import bresenham
import numpy.ma as ma
from scipy.optimize import curve_fit
from scipy import asarray as ar,exp
import cv2
get_ipython().run_line_magic('matplotlib', 'inline')
# Make sure you have collaborator access - if not go set it up following these instructinos:
# https://sdss-marvin.readthedocs.io/en/stable/installation.html
config.access = 'collab'
# Choose the data release you would like to use (could also use MPL)
config.setRelease('MPL-8')
print(config.access)

# I had to re-log in when I first ran this code:
config.login(refresh=True)
print(config.token)

#Strong bars:
#8444-12703
#8984-12704
#8135-6103
#9196-12701
#9886-12705

#Weak bars (some of these might not have enough signal to do anything with their velocity fields. Let me know if so and I’ll find some others):
#8987-12702
#8567-12703
#8318-12704
#8988-12705
#8332-12701

# Check to see if this map loads remotely, if so you are good to go- 
plateids = ['8332-12701']# this is my favorite galaxy
Maps(plateids[0]).download()

maps = [Maps(plateid) for plateid in plateids]
# .datamodel shows all the options for various map extensions
#print(Maps('7443-12703').datamodel)

print(maps[0].header['REFF'], 'Reff')
print(maps[0].header['ECOOELL'], 'ellipticity')
print(maps[0].header['ECOOPA'], 'PA')



svel_maps = [mapz['stellar_vel'].value for mapz in maps]# this is in case you want to feed it a list
svel_map = svel_maps[0]
svel_map = ma.masked_where(svel_map==0, svel_map)

svel_e_maps = [1/np.sqrt(mapz['stellar_vel'].ivar) for mapz in maps]# this is in case you want to feed it a list
svel_e_map = svel_e_maps[0]

# sigma_star_corr = sqrt( STELLAR_SIGMA2 - STELLAR_SIGMACORR2 )
ssig_maps = [np.sqrt(mapz['stellar_sigma'].value**2-mapz['stellar_sigmacorr'].value**2) for mapz in maps]# this is in case you want to feed it a list
ssig_map = ssig_maps[0]
ssig_map = ma.masked_where(ssig_map==0, ssig_map)


ssig_e_maps = [1/np.sqrt(mapz['stellar_sigma'].ivar) for mapz in maps]# this is in case you want to feed it a list
ssig_e_map = ssig_e_maps[0]

plt.clf()
plt.imshow(svel_map, cmap='RdBu_r')
plt.title(r'$V_*$')
plt.colorbar()
plt.show()

plt.clf()
plt.imshow(svel_e_map, cmap='RdBu_r')
plt.title(r'Error $V_*$')
plt.colorbar()
plt.show()

plt.clf()
plt.imshow(ssig_map, vmax=200, cmap='plasma')
plt.title(r'$\sigma_*$')
plt.colorbar()
plt.show()

plt.clf()
plt.imshow(ssig_e_map,vmax=200, cmap='plasma')
plt.title(r'Error $\sigma_*$')
plt.colorbar()
plt.show()


# In[161]:


'''
STILL TO DO AND THIS IS HUGE - 
YOU NEED TO CHECK THE MANGA IMAGING PARAMETERS LIKE R_EFF AND ELLIPTICITY AGAINST THOSE YOU DERIVE
FROM YOUR IMAGING PIPELINE.
'''


# In[171]:


print(np.shape(svel_map))



# First get effective radius, PA, and ellipticity from MaNGA, then run GALFIT on em

#ECOOPA	Position angle used for the semi-major axis polar coordinate calculations
#ECOOELL	Ellipticity (1-b/a) used for the semi-major axis polar coordinate calculations

R_e_MaNGA = maps[0].header['REFF']
ellip_MaNGA = maps[0].header['ECOOELL']
PA_MaNGA = maps[0].header['ECOOPA']


'''img_params=extract_GALFIT_parameters(view, myr, im, run)
#return PA_img, size_a, arcs_totes, inc, pixelscale, r_e
if img_params[0]==0:#this means the GALFIT file doesn't exist
    STOP
#    continue



try:
    #I have found the PA from GALFIT to be better than the PA from statmorph 
    #which is commented out below
    PA_imag=img_params[0]
except TypeError:
     STOP'''
#   continue


epsilon = ellip_MaNGA#1-img_params[1]
r_e = R_e_MaNGA#img_params[2]

#print('deltapos', deltapos)
#continue

#epsilon = float(kin_cube[0].header['ellip']), again the epsilon from GALFIT is better
#r_e = float(kin_cube[0].header['REFF'])
#PA_img = float(kin_cube[0].header['PA_img']), again, the PA from GALFIT seems to be better

#axial ratio (q) = b/a = semi-minor/semi-major
#ellipticity (ε) = 1 – q



#from the r-band image you need to run statmorph and get both the ellipticity and the r-band effective radius
#we already measured this from the r-band image in the kinematic program so I'd like to extract those two 
#parameters from the header of the stellar_kinematics file



size=np.shape(svel_map)[0]

# Prepare kin_cube, which has four dimensions - vel, vel_e, sig, and sig_e

kin_cube=[[svel_map],[svel_e_map],[ssig_map],[ssig_e_map]]
coords=map_to_coords(kin_cube, size)

rad = radon_python_mod(svel_map, 30, 30, r_e, 1, 'yes', 'no')#was 30,30
if rad[8]==1:#this means we gotta expand the grid
    rad = radon_python_mod(vel, 30, 30, r_e, 2, 'no')



input_kinemetry(plateids[0], coords[0], coords[1], coords[2], coords[3], coords[4], coords[5] , size/2-rad[0][1], size/2-rad[0][0],  'no', size)


# In[ ]:





# In[ ]:





# In[ ]:





# In[145]:


def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]
def gauss(x,a,x0,sigma,offset):
    return a*exp(-(x-x0)**2/(2*sigma**2))+offset

def ndim_grid(start,stop):
    # Set number of dimensions
    ndims = len(start)

    # List of ranges across all dimensions
    L = [np.arange(start[i],stop[i]) for i in range(ndims)]

    # Finally use meshgrid to form all combinations corresponding to all 
    # dimensions and stack them as M x ndims array
    return np.hstack((np.meshgrid(*L))).swapaxes(0,1).reshape(ndims,-1).T



def extract_GALFIT_parameters(view, myr, im, run):
    """
    This retrives the GALFIT predictors 
    """
    output='../LAURA_Sims/GALFIT_folder/out_'+str(run)+'_'+str(view)+'_'+str(myr)+'.fits'
    try:
        out=pyfits.open(output)
    except FileNotFoundError:
        #print('NO GALFIT FILEEEE')
        STOP
        
        return 0, 0, 0

    

    try:
        try:
            attempt=float(out[2].header['1_MAG'][0:5])
            attempt_2=float(out[2].header['2_MAG'][0:5])
        except ValueError:
            return 0, 0, 0
        if float(out[2].header['1_MAG'][0:5]) < float(out[2].header['2_MAG'][0:5]):
        #this means the 1st one is brighter
            inc=float(out[2].header['1_AR'][0:5])
            r_e=float(out[2].header['1_RE'][0:5])
            PA_img=(float(out[2].header['1_PA'][:5]))
            
        else:
            inc=float(out[2].header['2_AR'][0:5])
            r_e=float(out[2].header['2_RE'][0:5])
            PA_img=(float(out[2].header['2_PA'][:5]))
            
    except KeyError or ValueError:#if there is no #2
        try:
            
            inc=float(out[2].header['1_AR'][0:5])
            r_e=float(out[2].header['1_RE'][0:5])
            PA_img=(float(out[2].header['1_PA'][:5]))
        except ValueError:
            return 0, 0, 0
        
    return PA_img, inc, r_e

def map_to_coords(map_cube, size):
    """
    Converts coordinates to list form (in order to feed through kinemetry).
    """
    x_list=[]
    y_list=[]
    vel_list=[]
    vel_e_list=[]
    sig_list=[]
    sig_e_list=[]
    
    vel_dimension=np.reshape(map_cube[0],(size,size))

    vel_e_dimension=np.reshape(map_cube[1],(size,size))
    sig_dimension=np.reshape(map_cube[2],(size,size))
    sig_e_dimension=np.reshape(map_cube[3],(size,size))
    
    for i in range(size):
        for j in range(size):
            try:
                value = vel_dimension[i,j]
                if str(value) == '--':
                    continue
                vel_list.append(vel_dimension[i,j])
                x_list.append(i)
                y_list.append(j)
                vel_e_list.append(vel_e_dimension[i,j])
                sig_list.append(sig_dimension[i,j])
                sig_e_list.append(sig_e_dimension[i,j])
            except IndexError:
                continue
            
    
    
    return x_list, y_list, vel_list, vel_e_list, sig_list, sig_e_list


def radon_python_mod(vel_field, n_p, n_theta, r_ap, factor, plot, trouble):
    
    
    """
    This section performs the radon transform, from Stark et al. 2018.
    
    It is a long calculation, because it first calculates the Absolute Radon Transform (R_AB) of the velocity field, 
    then it iterates for multiple choices of the kinematic center. It determines the kinematic center by minimizing
    the asymmetry, A, of the Radon profile calculated from R_AB using a centroiding method.
    
    If the kinematic center is on the edge of the search grid of spaxels, it throws a flag and the code will rerun
    this function after expanding the grid.
    
    """
    #It first converts the x and y coordinates into p and theta coordinates (circular)
    #p is rho, which is the distance of the point on the velocity field from the kinematic center
    p_list = np.linspace(-int(np.shape(vel_field)[0]/2)+5,int(np.shape(vel_field)[0]/2)-5, int(np.shape(vel_field)[0]/2)+1)#was 5
    p_list = np.linspace(-int(np.shape(vel_field)[0]/2)+5, int(np.shape(vel_field)[0]/2)-5,n_p)#was 20
   
    #theta is the angle from the negative y axis CCW to the point on the velocity map.
    theta_list = np.linspace(0, 180, n_theta)#was 10
    
    #It searches over a grid of coordinates around the photometric center to find the 'kinematic center',
    #so it creates a 3x3 grid (from -3 to 3 including 0) in 2D
    box_list=list(ndim_grid([-3, -3],[4,4]))
    #If the kinematic center is not found on the first iteration, it expands the dimensions of the grid
    #by a factor of 2 upon rerunning.
    box_list=[factor*x for x in box_list]
    
    
    
    #Here I create a X and Y meshgrid type list of these points, since this is a rough approximation of the 
    #full grid of points in order to later plot what is happening.
    X_list=[]
    Y_list=[]

    for j in range(len(box_list)):

        X_list.append(int(np.shape(vel_field)[0]/2+box_list[j][0]))#-10+box_list[b][0])
        Y_list.append(int(np.shape(vel_field)[1]/2+box_list[j][1]))#-10+box_list[b][1])
        
    #print('X_list', X_list)
    #print('Y_list', Y_list)
    
    #creates empty lists that will be populated with R_AB, and all other derived quantities from this for every 
    #rho, theta point.
    A_list=[]
    A_e_list=[]
    A_2_list=[]
    R_AB_list=[]
    theta_hat_list=[]
    theta_hat_e_list=[]
    
    
    #First run it for just the center one, index b=24 in order to normalize relative to the center
    #later on in the calculation of A.
    R_AB=[]
    b=24
    
    for i in range(len(p_list)):
        for j in range(len(theta_list)):
            #
            X = int(p_list[i]*math.cos(math.radians(theta_list[j]))+np.shape(vel_field)[0]/2-0.5+box_list[b][0])#-10+box_list[b][0])
            Y = int(p_list[i]*math.sin(math.radians(theta_list[j]))+np.shape(vel_field)[1]/2-0.5+box_list[b][1])#-10+box_list[b][1])
            

            '''We have an X and a Y and a theta (slope) so we should be able to get the intercept'''
            '''And then two more points on either end'''


            '''But we only want to calculate for things that are on the circle'''

            try:
                #if this point exists in the velocity field then you can continue
                test_value = vel_field[X,Y]
            except IndexError:
                R_AB.append(-1000)
                continue
            if np.isnan(vel_field[X,Y]):
                R_AB.append(-1000)
                STOP2
                continue

            if str(vel_field[X,Y]) == '--':
                R_AB.append(-1000)
                continue
            #calculate the slope of the line segment from the kinematic center (in this case the photometric center)
            #to the given point
            deltay = Y - np.shape(vel_field)[1]/2
            deltax = X - np.shape(vel_field)[0]/2
            #draw a line perpendicular to this; the radon transform will be calculated along this line
            slope_p = math.tan(math.radians(theta_list[j]+90))#-deltax/deltay
            #draw a line from the point to where it intersects the bottom left of the map, which is the new origin
            intercept = Y - slope_p*X



            if slope_p > 1000:
                #vertical, so calculate along one value of X for the entire length of y
                x_min = X
                x_max = X
                y_min = 0
                y_max = np.shape(vel_field)[0]
            else:
                x_min = 0
                x_max = np.shape(vel_field)[0]
                y_min = intercept
                y_max = intercept+slope_p*x_max

            #This neat line draws a line through a given set of coordinates
            bres_list = list(bresenham(int(x_min), int(y_min), int(x_max), int(y_max)))

            #to calculate the absolute Bounded Radon Transform, do this for all points that are within r_e/2 of the center
            #of the line:
            vel_append=[]
            for k in range(len(bres_list)):
                if bres_list[k][0] < 0 or bres_list[k][1] < 0:
                    continue
                if np.sqrt((bres_list[k][0]-X)**2+(bres_list[k][1]-Y)**2) > r_e/2:
                    continue
               
                try:
                    vel_append.append(vel_field[bres_list[k][1], bres_list[k][0]])
                    
                except IndexError:
                    continue

            #clean it up, no masked values in here!
            vel_append_clean=[]
            for k in range(len(vel_append)):
                if ma.is_masked(vel_append[k]):
                    continue
                else:
                    vel_append_clean.append(vel_append[k])


            #finally, create R_AB by summing all of these velocity differences
            if vel_append_clean:
                inside=vel_append_clean-np.mean(vel_append_clean)
                R_AB.append(np.sum(np.abs(inside)))
            else:
                R_AB.append(-1000)


    R_AB=ma.masked_where(np.isnan(R_AB), R_AB)
    R_AB = ma.masked_where(R_AB==-1000, R_AB)

    
    R_AB_array = np.reshape(R_AB, (len(p_list),len(theta_list)))

    
    '''Try to plot R_AB values'''
    plt.clf()
    plt.imshow(R_AB_array)
    locs, labels = plt.yticks()
    #p_list = np.linspace(-int(np.shape(vel_field)[0]/2)+5, int(np.shape(vel_field)[0]/2)-5,n_p)#was 20
   
    ylabels = [p_list[0], 0, p_list[-1]]
    ylocs = [locs[0], locs[-1]/2, locs[-1]] 
    plt.yticks(ylocs, ylabels)
    
    locs, labels = plt.xticks()
    #p_list = np.linspace(-int(np.shape(vel_field)[0]/2)+5, int(np.shape(vel_field)[0]/2)-5,n_p)#was 20
   
    xlabels = [theta_list[0], theta_list[-1]/2, theta_list[-1]]
    xlocs = [locs[0], locs[-1]/2, locs[-1]] 
    plt.xticks(xlocs, xlabels)
    plt.colorbar()
    plt.show()

    #Now, extract the R_AB value at each rho value across all theta values.
    #This creates the Radon profile; the estimated value of theta hat that minimizes R_AB at each value of rho.
    #We minimize R_AB because we want the theta value at each rho that 

    #these are the estimated values of theta that are best fit for a value of rho
    theta_hat=[]
    theta_hat_e=[]

    for l in range(len(p_list)):

        marginalized = R_AB_array[l,:]
        marginalized = ma.masked_where(marginalized<1e-4, marginalized)
        count = len([i for i in marginalized if i > 1e-3])
        #count up how many elements are in the row of R_AB --> if it is less than 6, don't measure it
        #because it will cause an error when trying to fit a Gaussian, k = n+1
        if count < 6:
            theta_hat.append(0)
            theta_hat_e.append(0)
            continue

        if ma.is_masked(marginalized)==True:
            theta_hat.append(0)
            theta_hat_e.append(0)
            continue

        try:
            #initially, try to fit a negative gaussian to determine theta hat
            popt,pcov = curve_fit(gauss,theta_list,marginalized,p0=[-abs(np.min(marginalized)-np.max(marginalized)),theta_list[list(marginalized).index(np.min(marginalized))],20,np.mean(marginalized)])
            append_value = popt[1]
            
            
            
        except RuntimeError or OptimizeError: 
            theta_hat.append(0)
            theta_hat_e.append(0)

            continue
            
            

        #if the code fits a positive Gaussian, try a new guess for the correct theta hat, the second smallest value       
        if (popt[0]>0):
            try:
                popt,pcov = curve_fit(gauss,theta_list,marginalized,p0=[-abs(np.min(marginalized)-np.max(marginalized)),theta_list[list(marginalized).index(sorted(marginalized)[1])],20,np.mean(marginalized)])
                append_value = popt[1]
                if popt[0] > 0:
                    #if this doesn't work, quit and move on

                    theta_hat.append(0)
                    theta_hat_e.append(0)

                    continue

                
                

            except RuntimeError or OptimizeError :
                theta_hat.append(0)
                theta_hat_e.append(0)

                continue
        
        #sometimes, it is necessary to shift the xs because the peak is at 180, which is right at the edge
        if ((popt[1] - 3*np.sqrt(pcov[1][1])) < 0) or ((popt[1] + 3*np.sqrt(pcov[1][1])) > 180):
            theta_list_shifted = theta_list+find_nearest(theta_list,90)
            index_1 = list(theta_list).index(find_nearest(theta_list,90))
            new_marginalized = np.concatenate((marginalized[index_1:], marginalized[:index_1]))

            try:
                popt,pcov = curve_fit(gauss,theta_list_shifted,new_marginalized,p0=[-abs(np.min(new_marginalized)-np.max(new_marginalized)),theta_list_shifted[list(new_marginalized).index(np.min(new_marginalized))],20,np.mean(new_marginalized)])
                

                if popt[0] > 0:
                    theta_hat.append(0)
                    theta_hat_e.append(0)

                    continue
                    
                    
                if popt[1] > 180:
                    append_value = popt[1]-180
                else:
                    append_value = popt[1]

            except RuntimeError or OptimizeError :
                theta_hat.append(0)
                theta_hat_e.append(0)

                continue
            
        theta_hat.append(append_value)
        theta_hat_e.append(np.sqrt(pcov[1][1]))

    if trouble=='yes':
        plt.clf()
        plt.scatter(p_list, theta_hat)
        plt.errorbar(p_list, theta_hat, yerr=theta_hat_e)
        plt.show()
        
    #now to calculate A, it is necessary to sum the values of theta hat at a mirror image of themselves
    #across p=0
    delta_theta_sum=[]
    delta_theta_sum_e=[]


    for l in range(int(len(p_list)/2)):
        if (theta_hat[0+l]==0) or (theta_hat[-1-l]==0) or (abs(theta_hat[0+l])>180) or (abs(theta_hat[-1-l]) > 180):

            delta_theta_sum.append(0)
            delta_theta_sum_e.append(0)
        else:
            if abs(theta_hat[0+l]-theta_hat[-1-l]) > 90:
                #because the maximum you can be apart is 90
                inside = 180 - abs(theta_hat[0+l]-theta_hat[-1-l])
            else:
                inside = abs(theta_hat[0+l]-theta_hat[-1-l])

            delta_theta_sum.append(inside)
            #I would also like to have an error estimate on this quantity:
            delta_theta_sum_e.append(np.sqrt((theta_hat_e[0+l])**2+
                                         (theta_hat_e[-1-l])**2))
    


    delta_theta_sum_masked=ma.masked_where(np.array(delta_theta_sum)==0, delta_theta_sum)

    

    OG_weight=ma.count(delta_theta_sum_masked)
    
    
    plt.clf()
    fig=plt.figure()
    ax0 = fig.add_subplot(111)
    if plot=='yes':
        
        im0 = ax0.imshow(vel_field, cmap='RdBu_r')
    
    
    #Okay now do this for all the other positions in box_list
    A_list=[]
    A_e_list=[]
    
    A_2_list=[]
    R_AB_list=[]
    theta_hat_list=[]
    theta_hat_e_list=[]

    
            
    A_list=[]
    A_e_list=[]
    A_2_list=[]
    R_AB_list=[]
    theta_hat_list=[]
    theta_hat_e_list=[]
    
    X_list_real=[]
    Y_list_real=[]
    

    for b in range(len(box_list)):
        R_AB=[]
        X_list_real.append(np.shape(vel_field)[0]/2-0.5+box_list[b][0])
        Y_list_real.append(np.shape(vel_field)[1]/2-0.5+box_list[b][1])
        for i in range(len(p_list)):
            for j in range(len(theta_list)):
                #
                X = int(p_list[i]*math.cos(math.radians(theta_list[j]))+np.shape(vel_field)[0]/2-0.5+box_list[b][0])#-10+box_list[b][0])
                Y = int(p_list[i]*math.sin(math.radians(theta_list[j]))+np.shape(vel_field)[1]/2-0.5+box_list[b][1])#-10+box_list[b][1])
                
                try:
                    #if this point exists in the velocity field then you can continue
                    test_value = vel_field[X,Y]
                except IndexError:
                    R_AB.append(-1000)
                    continue
                if np.isnan(vel_field[X,Y]):
                    R_AB.append(-1000)
                    STOP2
                    continue

                if str(vel_field[X,Y]) == '--':
                    R_AB.append(-1000)
                    continue
                #calculate the slope of the line segment from the kinematic center (in this case the photometric center)
                #to the given point
                deltay = Y - np.shape(vel_field)[1]/2
                deltax = X - np.shape(vel_field)[0]/2
                #draw a line perpendicular to this; the radon transform will be calculated along this line
                slope_p = math.tan(math.radians(theta_list[j]+90))#-deltax/deltay
                #draw a line from the point to where it intersects the bottom left of the map, which is the new origin
                intercept = Y - slope_p*X



                if slope_p > 1000:
                    #vertical, so calculate along one value of X for the entire length of y
                    x_min = X
                    x_max = X
                    y_min = 0
                    y_max = np.shape(vel_field)[0]
                else:
                    x_min = 0
                    x_max = np.shape(vel_field)[0]
                    y_min = intercept
                    y_max = intercept+slope_p*x_max

                #This neat line draws a line through a given set of coordinates
                bres_list = list(bresenham(int(x_min), int(y_min), int(x_max), int(y_max)))

                
                #to calculate the absolute Bounded Radon Transform, do this for all points that are within r_e/2 of the center
                #of the line:
                vel_append=[]
                for k in range(len(bres_list)):
                    if bres_list[k][0] < 0 or bres_list[k][1] < 0:
                        continue
                    if np.sqrt((bres_list[k][0]-X)**2+(bres_list[k][1]-Y)**2) > r_e/2:
                        continue

                    try:

                        vel_append.append(vel_field[bres_list[k][1], bres_list[k][0]])
                        #vel_new[bres_list[j][1], bres_list[j][0]] = vel_field[bres_list[j][1], bres_list[j][0]]
                    except IndexError:
                        continue

                #clean it up, no masked values in here!
                vel_append_clean=[]
                for k in range(len(vel_append)):
                    if ma.is_masked(vel_append[k]):
                        continue
                    else:
                        vel_append_clean.append(vel_append[k])


                #finally, create R_AB by summing all of these velocity differences
                if vel_append_clean:
                    inside=vel_append_clean-np.mean(vel_append_clean)
                    R_AB.append(np.sum(np.abs(inside)))
                else:
                    R_AB.append(-1000)


        R_AB=ma.masked_where(np.isnan(R_AB), R_AB)
        R_AB = ma.masked_where(R_AB==-1000, R_AB)


        R_AB_array = np.reshape(R_AB, (len(p_list),len(theta_list)))


        #Now, extract the R_AB value at each rho value across all theta values.
        #This creates the Radon profile; the estimated value of theta hat that minimizes R_AB at each value of rho.
        #We minimize R_AB because we want the theta value at each rho that 

        #these are the estimated values of theta that are best fit for a value of rho
        theta_hat=[]
        theta_hat_e=[]

        for l in range(len(p_list)):

            marginalized = R_AB_array[l,:]
            marginalized = ma.masked_where(marginalized<1e-4, marginalized)
            count = len([i for i in marginalized if i > 1e-3])
            #count up how many elements are in the row of R_AB --> if it is less than 6, don't measure it
            #because it will cause an error when trying to fit a Gaussian, k = n+1
            if count < 6:
                theta_hat.append(0)
                theta_hat_e.append(0)
                continue

            if ma.is_masked(marginalized)==True:
                theta_hat.append(0)
                theta_hat_e.append(0)
                continue

            try:
                popt,pcov = curve_fit(gauss,theta_list,marginalized,p0=[-abs(np.min(marginalized)-np.max(marginalized)),theta_list[list(marginalized).index(np.min(marginalized))],20,np.mean(marginalized)])
                append_value = popt[1]

                

            #sometimes, it is necessary to shift the xs because the peak is at 180, which is right at the cutoff
            except RuntimeError or OptimizeError: 
                theta_hat.append(0)
                theta_hat_e.append(0)

                continue
                


                    
            if (popt[0]>0):
                try:
                    popt,pcov = curve_fit(gauss,theta_list,marginalized,p0=[-abs(np.min(marginalized)-np.max(marginalized)),theta_list[list(marginalized).index(sorted(marginalized)[1])],20,np.mean(marginalized)])
                    append_value = popt[1]
                    if popt[0] > 0:


                        theta_hat.append(0)
                        theta_hat_e.append(0)

                        continue


                except RuntimeError or OptimizeError :
                    theta_hat.append(0)
                    theta_hat_e.append(0)

                    continue
            if ((popt[1] - 3*np.sqrt(pcov[1][1])) < 0) or ((popt[1] + 3*np.sqrt(pcov[1][1])) > 180):
                theta_list_shifted = theta_list+find_nearest(theta_list,90)
                index_1 = list(theta_list).index(find_nearest(theta_list,90))
                new_marginalized = np.concatenate((marginalized[index_1:], marginalized[:index_1]))

                try:
                    popt,pcov = curve_fit(gauss,theta_list_shifted,new_marginalized,p0=[-abs(np.min(new_marginalized)-np.max(new_marginalized)),theta_list_shifted[list(new_marginalized).index(np.min(new_marginalized))],20,np.mean(new_marginalized)])
                    
                    if popt[0] > 0:
                        theta_hat.append(0)
                        theta_hat_e.append(0)

                        continue

                    
                    if popt[1] > 180:
                        append_value = popt[1]-180
                    else:
                        append_value = popt[1]

                except RuntimeError or OptimizeError :
                    theta_hat.append(0)
                    theta_hat_e.append(0)

                    continue

            theta_hat.append(append_value)
            theta_hat_e.append(np.sqrt(pcov[1][1]))



 

        delta_theta_sum=[]
        delta_theta_sum_e=[]
        
        
        for l in range(int(len(p_list)/2)):
            if (theta_hat[0+l]==0) or (theta_hat[-1-l]==0) or (abs(theta_hat[0+l])>180) or (abs(theta_hat[-1-l]) > 180):
                
                delta_theta_sum.append(0)
                delta_theta_sum_e.append(0)
            else:
                
                
                if abs(theta_hat[0+l]-theta_hat[-1-l]) > 90:
                    #because the maximum you can be apart is 90
                    inside = 180 - abs(theta_hat[0+l]-theta_hat[-1-l])
                else:
                    inside = abs(theta_hat[0+l]-theta_hat[-1-l])
                
                delta_theta_sum.append(inside)
                delta_theta_sum_e.append(np.sqrt((theta_hat_e[0+l])**2+
                                             (theta_hat_e[-1-l])**2))
        
            
            
            
        theta_hat_mask = ma.masked_where(theta_hat==0, theta_hat)
        theta_hat_mask = ma.masked_where(abs(theta_hat_mask) >180, theta_hat_mask)
        
        theta_hat_list.append(theta_hat_mask)
        theta_hat_e_list.append(theta_hat_e)
        
        delta_theta_sum_masked=ma.masked_where(np.array(delta_theta_sum)==0, delta_theta_sum)

        #A is weighted by the A value of the center of the map
        A = ((ma.sum(delta_theta_sum)/(ma.count(delta_theta_sum_masked)**2)))*OG_weight
       
        A_percent_e = []
        for l in range(len(delta_theta_sum_e)):
            if delta_theta_sum[l] != 0 :
                A_percent_e.append((delta_theta_sum_e[l])**2)
        A_abs_error = np.sqrt(np.sum(A_percent_e))#so this is on the numerator of the A error only
        A_error = (A_abs_error/ma.sum(delta_theta_sum))*A
       
        
        
        A_list.append(A)
        
        
        A_e_list.append(A_error)

        #Also calculates the other types of asymmetry:
        #A_2
        delta_theta_sum=ma.masked_where(np.array(delta_theta_sum)==0, delta_theta_sum)
        delta_theta_sum_e=ma.masked_where(np.array(delta_theta_sum_e)==0, delta_theta_sum_e)
        
        
        A_2 = ma.sum(delta_theta_sum/delta_theta_sum_e)
        
        
        
        A_2_list.append(A_2)
        
        R_AB_list.append(R_AB_array)
        
    
      
    
    
    A_list = np.array(A_list)
    A_list = ma.masked_where(np.isnan(A_list), A_list)
    
    A_list_array = np.reshape(np.array(A_list), (7,7))
    A_e_list_array = np.reshape(np.array(A_e_list), (7,7))

    
    #first_number = 33-3*factor#33-3*factor = 27 if factor = 2
    
    #print('first number', first_number, factor)
    first_number = X_list_real[0]
    x = factor*np.linspace(0, np.shape(A_list_array)[0]-1,np.shape(A_list_array)[0])+first_number
    #y = (factor*6-factor*np.linspace(0, np.shape(A_list_array)[1]-1, np.shape(A_list_array)[1]))+first_number
    y = (factor*np.linspace(0, np.shape(A_list_array)[1]-1, np.shape(A_list_array)[1]))+first_number
    
    print('matching x and y with X and Y real', x, X_list_real)
    
    if min(X_list_real) != min(x):
        #this is improper scaling
        print('x and y', x, y)
    
        print('this is what it should match X_list', X_list_real, Y_list_real)
    
        STOP
    
    
    x, y = np.meshgrid(x, y)
    
    
    
    
    ax0.imshow(A_list_array, cmap='plasma')
    cs = ax0.contour(x,y, A_list_array, 5, colors='orange')#was 5
    p = cs.levels#collections[0].get_paths()[0]
    
    if plot=='yes':
        im1 = ax0.scatter(np.flip(Y_list_real, axis=0),X_list_real,c=A_list, s=30, zorder=100)

        plt.colorbar(im0)
        
    #Now, snatch the last level and use it to blur everything else out and find the center
    #now fill with zeros and do your best to turn into a binary bitmask'''
    if trouble=='yes':
        print(p[0], p)
    ret,thresh = cv2.threshold(A_list_array,p[-1]/2,100,cv2.THRESH_BINARY_INV)
    
    #ret, thresh = cv2.adaptiveThreshold(A_list_array, 50, cv2.ADAPTIVE_THRESH_GAUSSIAN_C, cv2.THRESH_BINARY_INV, 115, 1)
    M = cv2.moments(thresh)
    
    if trouble=='yes':
        plt.clf()
        plt.imshow(A_list_array)
        plt.colorbar()
        plt.show()

        print('moments', M)
        print('ret', ret)

       
    
        plt.clf()
        plt.imshow(thresh)
        plt.colorbar()
        plt.show()

    # calculate x,y coordinate of center
    cX = int(M["m10"] / M["m00"])
    cY = int(M["m01"] / M["m00"])

    
    
    
    
    #transform to flip in the x coordinate around the center coordinate, which changes based on the factor
    cX_flip = 6-cX
    
    x_match = cX_flip+first_number#33-3*factor
    y_match = (cY)+first_number#33-3*factor
    
    
    
    
    
    if plot=='yes':
        

        print('this is the cX and cY', x_match, y_match)


        plt.scatter( x_match, y_match, marker='x', color='red', zorder=105)
        plt.axhline(y=np.shape(vel_field)[0]/2, color='red')
        plt.axvline(x=np.shape(vel_field)[0]/2, color='red')

        plt.show()
    
    
    
    
        plt.clf()
        plt.imshow(thresh)
        plt.scatter(cX,cY, marker='x', color='red')

        x = np.linspace(0, np.shape(A_list_array)[0]-1,np.shape(A_list_array)[0])
        y = (np.linspace(0, np.shape(A_list_array)[1]-1, np.shape(A_list_array)[1]))
        x, y = np.meshgrid(x,y)

        plt.contour(x,y, A_list_array, 5, colors='orange')#was 5
        plt.xlim([6,0])

        plt.show()
    
    
    
    
    
    '''now find the min index of this'''
    for ii in range(len(X_list)):
        if X_list[ii]==find_nearest(np.array(X_list),x_match) and Y_list[ii]==find_nearest(np.array(np.flip(Y_list, axis=0)),y_match):
            min_index = ii
        
    
    
    
    if factor ==1:
        if abs(box_list[min_index][0])>2 or abs(box_list[min_index][1])>2:
            #then we are on the edge
            expand=1
        else:
            expand=0
    else:
        if abs(box_list[min_index][0])>5 or abs(box_list[min_index][1])>5:
            min_index=24
        expand=0
        
    
    
    if plot=='yes':
        plt.clf()
        plt.imshow(R_AB_list[min_index], cmap='summer')
        locs, labels = plt.yticks()
        #p_list = np.linspace(-int(np.shape(vel_field)[0]/2)+5, int(np.shape(vel_field)[0]/2)-5,n_p)#was 20

        ylabels = [p_list[0], 0, p_list[-1]]
        ylocs = [locs[0], locs[-1]/2, locs[-1]] 
        plt.yticks(ylocs, ylabels)

        locs, labels = plt.xticks()
        #p_list = np.linspace(-int(np.shape(vel_field)[0]/2)+5, int(np.shape(vel_field)[0]/2)-5,n_p)#was 20

        xlabels = [theta_list[0], theta_list[-1]/2, theta_list[-1]]
        xlocs = [locs[0], locs[-1]/2, locs[-1]] 
        plt.xticks(xlocs, xlabels)
        plt.colorbar()
        plt.show()
    
   
    return box_list[min_index],R_AB_list[min_index], A_list[min_index],  A_2_list[min_index],p_list, theta_list, theta_hat_list[min_index], theta_hat_e_list[min_index], expand

def input_kinemetry(plateid, x_list_after, y_list_after, vel_list_after, vel_e_list_after, sig_list_after, sig_e_list_after,  img_x, img_y,  plot, size):
    """
    Creates an input file for kinemetry (.txt file) that has columns of x and y coords, velocity, velocity dispersion,
    and the errors on both of these values.
    """
    file2=open('kinemetry_input_txt/kinemetry_input_'+str(plateid)+'.txt','w')
    file2.write('#'+'\t'+'XBIN'+'\t'+'YBIN'+'\t'+'VEL'+'\t'+'ER_VEL'+'\t'+'SIG'+'\t'+'ER_SIG'+'\n')
    
    #These are the coordinates of the kinematic center of the galaxy
    middle_x=img_x
    middle_y=img_y
    
    #for some reason, kinemetry wants each row indexed
    counter=1
    
    for i in range(len(x_list_after)):
        if np.isnan(vel_list_after[i]):#mask these values by not including them in the file
            #file2.write(str(counter)+'\t'+str((middle_x-x_list_after[len(x_list_after)-1-i]))+'\t'+str((middle_y-y_list_after[len(x_list_after)-1-i]))+'\t')
            #file2.write('--'+'\t'+'--'+'\t'+'--'+'\t'+'--'+'\n')
        
            continue
        if vel_e_list_after[i]==0.0:#these are also list indices to mask
            #file2.write(str(counter)+'\t'+str((middle_x-x_list_after[len(x_list_after)-1-i]))+'\t'+str((middle_y-y_list_after[len(x_list_after)-1-i]))+'\t')
            #file2.write('--'+'\t'+'--'+'\t'+'--'+'\t'+'--'+'\n')
        
            continue
        file2.write(str(counter)+'\t'+str((middle_x-x_list_after[len(x_list_after)-1-i]))+'\t'+str((middle_y-y_list_after[len(x_list_after)-1-i]))+'\t')
        file2.write(str(vel_list_after[i])+'\t'+str(vel_e_list_after[i])+'\t'+str(sig_list_after[i])+'\t'+str(sig_e_list_after[i])+'\n')
        counter +=1

    file2.close()
    
    #optional, you can plot it to check
    if plot=='yes':
        file_check='kinemetry_input_txt/kinemetry_input_'+str(plateid)+'.txt'
        stel_vel=np.zeros((size,size))
    
        with open(file_check, 'r') as f:
            data = f.readlines()


            x_list=[]
            y_list=[]
            vel_list=[]
            

            for line in data:
                words = line.split()
                #print('words', words)

                if words[1]=='XBIN':
                    continue
                else:
                    x_list.append(float(words[1]))
                    y_list.append(float(words[2]))
                    vel_list.append(float(words[3]))
                    
                    stel_vel[int(float(words[1])-size/2),int(float(words[2])-size/2)] = float(words[3])
                    
        plt.clf()
        im0=plt.imshow(stel_vel, cmap='RdBu_r')
        plt.scatter(np.shape(stel_vel)[0]/2, np.shape(stel_vel)[1]/2, color='red', marker='x')
        plt.axhline(y=np.shape(stel_vel)[0]/2, color='red')
        plt.axvline(x=np.shape(stel_vel)[0]/2, color='red')
    
        plt.colorbar(im0)
        plt.show()
print('compiled')


# In[ ]:




