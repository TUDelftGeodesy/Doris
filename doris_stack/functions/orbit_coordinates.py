# The three functions here are used to transform coordinates and times in the original xml files to a usefull format for doris.
# Therefore these functions are generally used to read metadata of sentinel files.
import math
import numpy as np

def lph2xyz(line,pixel,container,norm_orbit_line,centroid_lon,centroid_lat,height):

    # initialization
    MAXITER=500
    CRITERPOS=1.0e-16
    # WGS84 Elliposid:
    ellipsoid=[6378137.0 , 6356752.3141]
    num_points = np.array([line]).shape[0]
    xyz       = np.zeros((num_points,3))
    #$$$ parameter of the image
    tr1  = float(container['rangeTimePix'][0])/2 # one way range time [sec]
    RSR  = float(container['rangeRSR'][0])*2 # one way in [HZ]
    centerphi = float(centroid_lat)
    centerlambda= float(centroid_lon)
    SOL=299792458
    # $$$ reference surace: WGS84

    ell_a  = ellipsoid[0]                  # semimajor of the ellipsoid
    ell_b  = ellipsoid[1]                  # semiminor of the ellipsoid
    ell_e2 = (ell_a**2-ell_b**2)/ell_a**2    # squared first eccentricity(derived)
    # $$$ ell_e2b=(ell_a^2-ell_b^2)/ell_b^2;  % squared second eccentricity(derived)

    # $$$ [lat,long,h] of scene center to [x,y,z]
    h            = height        # this is only for initial values

    centerphi    = centerphi*np.pi/180
    centerlambda = centerlambda*np.pi/180

    Ncenter      = ell_a/np.sqrt(1-ell_e2*(np.sin(centerphi)**2)) # eccentricity
    scenecenterx = (Ncenter+h)*np.cos(centerphi)*np.cos(centerlambda)
    scenecentery = (Ncenter+h)*np.cos(centerphi)*np.sin(centerlambda)
    scenecenterz = (Ncenter+h-ell_e2*Ncenter)*np.sin(centerphi)

    for n in range(0,num_points):   # loop through points

        posonellx = scenecenterx
        posonelly = scenecentery
        posonellz = scenecenterz
        ratime = tr1 + (pixel-1.0)/RSR

    #get position and velocity of the satellite
        possatx = norm_orbit_line[n,1]
        possaty = norm_orbit_line[n,2]
        possatz = norm_orbit_line[n,3]
        velsatx = norm_orbit_line[n,4]
        velsaty = norm_orbit_line[n,5]
        velsatz = norm_orbit_line[n,6]

        equationset=np.zeros((3,1))
        partialsxyz=np.zeros((3,3))
        for iter in range(1, MAXITER+1):

            #update equations and slove system
            dsat_Px = posonellx - possatx    #vector of 'satellite to P on ellipsoid'
            dsat_Py = posonelly - possaty
            dsat_Pz = posonellz - possatz

            equationset[0,0] = -(velsatx*dsat_Px+velsaty*dsat_Py+velsatz* dsat_Pz)

            equationset[1,0] = -(dsat_Px*dsat_Px+dsat_Py*dsat_Py+dsat_Pz*dsat_Pz-(SOL*ratime)**2)

            equationset[2,0] = -((posonellx*posonellx+posonelly*posonelly)/((ell_a+height)**2)+(posonellz/(ell_b+height))**2-1.0)

            partialsxyz[0,0] = velsatx
            partialsxyz[0,1] = velsaty
            partialsxyz[0,2] = velsatz
            partialsxyz[1,0] = 2*dsat_Px
            partialsxyz[1,1] = 2*dsat_Py
            partialsxyz[1,2] = 2*dsat_Pz
            partialsxyz[2,0] = (2*posonellx)/((ell_a+height)**2)
            partialsxyz[2,1] = (2*posonelly)/((ell_a+height)**2)
            partialsxyz[2,2] = (2*posonellz)/((ell_b+height)**2)

        # solve system [NOTE] orbit_normalized, otherwise close to
        # singular
            solpos = np.linalg.solve(partialsxyz,equationset)

            solx = solpos[0,0]
            soly = solpos[1,0]
            solz = solpos[2,0]

        # update solution
            posonellx = posonellx + solx
            posonelly = posonelly + soly
            posonellz = posonellz + solz

        # check convergence
            if (abs(solx)<CRITERPOS and abs(soly)<CRITERPOS and abs(solz)<CRITERPOS):
                break
            elif(iter>=MAXITER):
                MAXITER=MAXITER+1
        # final solution: array of XYZ coordinates
        xyz[n,:]=np.array([posonellx, posonelly, posonellz]).copy()

        return xyz


def xyz2ell(position):

    ellipsoid=[6378137.0 , 6356752.3141]
    ell_a = ellipsoid[0]
    ell_b = ellipsoid[1]

    ell_e2  = (ell_a**2-ell_b**2)/ell_a**2    #squared first eccentricity(derived)
    ell_e2b = (ell_a**2-ell_b**2)/ell_b**2    # squared second eccentricity(derived)

    posx = position[:,0]
    posy = position[:,1]
    posz = position[:,2]

    r    = math.sqrt(posx**2 + posy**2)
    mu   = math.atan2(posz*ell_a, r*ell_b)

    sin3 = (math.sin(mu))**3
    cos3 = (math.cos(mu))**3
    phi  = math.atan2((posz + ell_e2b * ell_b * sin3),(r - ell_e2 * ell_a* cos3))

    Radar_lambda = math.atan2(posy,posx)
    N = ell_a / math.sqrt(1 - ell_e2 * (math.sin(phi))**2) #for every point no
                                                # approx with scene.center
    height = (r/math.cos(phi)) - N

    phi_lam_height = np.zeros(3)
    phi_lam_height[0] = phi*180/math.pi
    phi_lam_height[1] = Radar_lambda*180/math.pi
    phi_lam_height[2] = height

    return phi_lam_height

def intrp_orbit(line,container,burst_number):

    intrpOrder = 'spline'
    orbit_time       = np.zeros(len(container['orbitTime']),dtype=np.float64)
    orbit_x = np.zeros(len(container['orbitTime']),dtype=np.float64)
    orbit_y = np.zeros(len(container['orbitTime']),dtype=np.float64)
    orbit_z = np.zeros(len(container['orbitTime']),dtype=np.float64)
    for row in range(len(container['orbitTime'])):

        #orbit_time_position=precorb[row]
        orbit_time[row]    =hms2sec(container['orbitTime'][row].split('T')[1])
        orbit_x[row]       =float( container['orbitX'][row])
        orbit_y[row]       =float( container['orbitY'][row])
        orbit_z[row]       =float( container['orbitZ'][row])

    # compute normalization factors
    px = orbit_time # time

    f  = min(px)
    g  = (max(px)-min(px))
    px = (px-f)/g
    # polyDegree
    polyDegree = 2

    coef_x1 = (np.polyfit(px,orbit_x,polyDegree));
    a = coef_x1[2]
    b = coef_x1[1]
    c = coef_x1[0]
    coef_x = [c/(g**2), b/g-(2*c*f)/(g**2), a-b*f/g+c*(f/g)**2]

    coef_y1 = (np.polyfit(px,orbit_y,polyDegree))
    a = coef_y1[2]
    b = coef_y1[1]
    c = coef_y1[0]
    coef_y = [c/(g**2), b/g-(2*c*f)/(g**2), a-b*f/g+c*(f/g)**2]

    coef_z1 = (np.polyfit(px,orbit_z,polyDegree));
    a = coef_z1[2]
    b = coef_z1[1]
    c = coef_z1[0]
    coef_z = [c/(g**2), b/g-(2*c*f)/(g**2), a-b*f/g+c*(f/g)**2]

    vel_x = np.polyval(np.polyder(coef_x),orbit_time)
    vel_y = np.polyval(np.polyder(coef_y),orbit_time)
    vel_z = np.polyval(np.polyder(coef_z),orbit_time)

    acc_x = np.kron(np.ones(len(container['orbitTime'])),np.polyder(np.polyder(coef_x)))
    acc_y = np.kron(np.ones(len(container['orbitTime'])),np.polyder(np.polyder(coef_y)))
    acc_z = np.kron(np.ones(len(container['orbitTime'])),np.polyder(np.polyder(coef_z)))
    #print('acc_x.shape=',acc_x.shape)

    # interpolated orbit
    norm_orbit = np.array([orbit_time, orbit_x,orbit_y,orbit_z,vel_x,  vel_y,  vel_z,acc_x,  acc_y,  acc_z]);

    # interpolated orbit for l_aztime
    PRF       = float(container['azimuthPRF'][0])   #[Hz]     # Pulse repeition frequency
    Taz_start = hms2sec(container['azimuthTimeStart'][burst_number].split('T')[1])

    ta1          = Taz_start  # start time in UTC [sec]

    l_aztime     = (line-1)/PRF + ta1
    pos_orbit_x = np.interp(l_aztime,orbit_time,orbit_x)
    pos_orbit_y = np.interp(l_aztime,orbit_time,orbit_y)
    pos_orbit_z = np.interp(l_aztime,orbit_time,orbit_z)

    vel_orbit_x = np.interp(l_aztime,orbit_time,vel_x)
    vel_orbit_y = np.interp(l_aztime,orbit_time,vel_y)
    vel_orbit_z = np.interp(l_aztime,orbit_time,vel_z)

    #acc = np.interp(orbit_time,[acc_x, acc_y, acc_z],    l_aztime,intrpOrder)
    acc_orbit_x = np.interp(l_aztime,orbit_time,acc_x)
    acc_orbit_y = np.interp(l_aztime,orbit_time,acc_y)
    acc_orbit_z = np.interp(l_aztime,orbit_time,acc_z)


    norm_orbit_line = np.array([l_aztime, pos_orbit_x,pos_orbit_y,pos_orbit_z,vel_orbit_x,vel_orbit_y,vel_orbit_z,acc_orbit_x,acc_orbit_y,acc_orbit_z])

    return norm_orbit.transpose(),norm_orbit_line.reshape(1,-1,order='F')

####################################################################################################

def hms2sec(hmsString,convertFlag='float'):
    # input hmsString syntax: XX:XX:XX.xxxxxx
    secString = int(hmsString[0:2])*3600 + \
        int(hmsString[3:5])*60 + \
        float(hmsString[6:])
    if convertFlag == 'int' :
        return int(secString)
    elif convertFlag == 'float' :
        return float(secString)
    else:
        return int(secString)
