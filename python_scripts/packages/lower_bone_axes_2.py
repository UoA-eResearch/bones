"""
Modified version of bone_axes.py created by Desney Geybe 
Made for the pelvis, femur, tibia-fibula. They are functions used to align the 
bones to their joint coordinate system (Wu, et al., 2002). 

Author: Alicia Lie 
"""
import numpy as np
import scipy
from numpy.linalg import norm
# from scipy.linalg import lstsq

from packages.calculate import Calculate
from packages.transform import Align, Transform
from gias2.common import geoprimitives as gp

#%%#########################################################################%%#
class PelvisAxes:
    ## Align eigenvectors with global coordinate system #######################
    @staticmethod
    def align_with_coordinate_system(data):
        # Define eigenvectors
        u = Calculate.eigenvectors(data)
        com = np.mean(data, axis=0)

        # Define the rotation and translation 4x4 matrices
        u_t = np.transpose(u)
        r = np.array([[u_t[0, 0], u_t[0, 1], u_t[0, 2], 0],
                      [u_t[1, 0], u_t[1, 1], u_t[1, 2], 0],
                      [u_t[2, 0], u_t[2, 1], u_t[2, 2], 0],
                      [0, 0, 0, 1]])

        t = np.array([[1, 0, 0, -com[0]],
                      [0, 1, 0, -com[1]],
                      [0, 0, 1, -com[2]],
                      [0, 0, 0, 1]])

        # Define the combined transformation and apply it
        tm = np.matmul(r, t)
        data = Transform.by_tm(data, tm)

        # Align first eigenvector with z-axis (instead of x) this is for the pelvis
        data, tm_temp = Transform.by_degrees(data, "y", 90)
        tm = tm_temp.dot(tm)

        return tm

    @staticmethod
    def set_origin(data, Pel_Landmarks):
        '''
        sets origin to left hip joint of pelvis
        '''
        vertices = Pel_Landmarks.vertices[:]
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        Pel_Landmarks.vertices = data
        # move origin
        origin = Pel_Landmarks.landmarks["left_hip_joint"]["coords"]
        data, tm_1 = Transform.by_translation(data, -origin)
        Pel_Landmarks.landmarks["left_hip_joint"]["coords"] = [0,0,0]
        for ld in Pel_Landmarks.landmarks:
            if ld == "left_hip_joint" or ld == 'lengths' or ld == 'angles':
                continue
            Pel_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Pel_Landmarks.landmarks[ld]["coords"], tm_1)
        tm = tm_1.dot(tm)
        Pel_Landmarks.vertices = vertices[:]
        return tm
    
    @staticmethod
    def set_origin_lasis(data, Pel_Landmarks):
        '''
        sets origin to left hip joint of pelvis
        '''
        vertices = Pel_Landmarks.vertices[:]
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        Pel_Landmarks.vertices = data
        # move origin
        origin = Pel_Landmarks.landmarks["ASIS_L"]["coords"]
        data, tm_1 = Transform.by_translation(data, -origin)
        Pel_Landmarks.landmarks["ASIS_L"]["coords"] = [0,0,0]
        for ld in Pel_Landmarks.landmarks:
            if ld == "ASIS_L" or ld == 'lengths' or ld == 'angles':
                continue
            Pel_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Pel_Landmarks.landmarks[ld]["coords"], tm_1)
        tm = tm_1.dot(tm)
        Pel_Landmarks.vertices = vertices[:]
        return tm
    
    @staticmethod
    def set_origin_midasis(data, Pel_Landmarks):
        '''
        sets origin to left hip joint of pelvis
        '''
        vertices = Pel_Landmarks.vertices[:]
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        Pel_Landmarks.vertices = data
        # move origin
        origin = Pel_Landmarks.landmarks["ASIS_mid"]["coords"]
        data, tm_1 = Transform.by_translation(data, -origin)
        Pel_Landmarks.landmarks["ASIS_mid"]["coords"] = [0,0,0]
        for ld in Pel_Landmarks.landmarks:
            if ld == "ASIS_mid" or ld == 'lengths' or ld == 'angles':
                continue
            Pel_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Pel_Landmarks.landmarks[ld]["coords"], tm_1)
        tm = tm_1.dot(tm)
        Pel_Landmarks.vertices = vertices[:]
        return tm

    ## Define the anteroposterior axis of the pelvis and align with the x axis ###
    @staticmethod
    def align_x_axis(data, Pel_Landmarks):

        # Restore landmark vertices
        vertices = Pel_Landmarks.vertices[:]

        # Initialise variables and find optimal rotation
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

        # Get location of landmarks
        Pel_Landmarks.vertices = data

        # Define posterior->anterior(x) axis
        ASIS_mid = Pel_Landmarks.landmarks["ASIS_mid"]["coords"]
        PSIS_mid = Pel_Landmarks.landmarks["PSIS_mid"]["coords"]
        AP_vector = ASIS_mid - PSIS_mid

        # Rotation about z-axis that aligns anterior posterior axis with x-axis
        vec_xy = np.array([AP_vector[0], AP_vector[1], 0])
        vec_xy = vec_xy / norm(vec_xy)
        gamma = np.arccos(np.dot(vec_xy, [1, 0, 0]))
        if vec_xy[1] > 0:
            gamma = -gamma
            # Apply rotation about the z axis
        data, tm_temp = Transform.by_radians(data, "z", gamma)
        for ld in Pel_Landmarks.landmarks:
            if ld == 'lengths' or ld == 'angles':
                continue
            Pel_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Pel_Landmarks.landmarks[ld]["coords"], tm_temp)
        tm = tm.dot(tm_temp)
        Pel_Landmarks.vertices = vertices[:]
        return tm

    @staticmethod
    def align_z_axis_x(data, Pel_Landmarks):
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        vertices = Pel_Landmarks.vertices[:]
        # Initialise variables and find optimal rotation
        # Get location of landmarks
        Pel_Landmarks.vertices = data

        # Define posterior->anterior(x) axis
        l_ASIS = Pel_Landmarks.landmarks["ASIS_L"]["coords"]
        r_ASIS = Pel_Landmarks.landmarks["ASIS_R"]["coords"]

        AP_vector = r_ASIS - l_ASIS
        # Rotation about y-axis that aligns asis axis with z-axis
        vec_yz = np.array([0, AP_vector[1], AP_vector[2]])
        vec_yz = vec_yz / norm(vec_yz)
        beta = np.arccos(np.dot(vec_yz, [0, 0, 1]))
        if vec_yz[0] > 0:
            beta = -beta
            # Apply rotation about the z axis
        data, tm_temp = Transform.by_radians(data, "x", beta)
        for ld in Pel_Landmarks.landmarks:
            if  ld == 'lengths' or ld == 'angles':
                continue
            Pel_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Pel_Landmarks.landmarks[ld]["coords"], tm_temp)
        tm = tm.dot(tm_temp)

        Pel_Landmarks.vertices = vertices[:]

        return tm
    @staticmethod
    def align_z_axis_y(data, Pel_Landmarks):

        vertices = Pel_Landmarks.vertices[:]
        # Initialise variables and find optimal rotation
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        # Get location of landmarks
        Pel_Landmarks.vertices = data

        # Define posterior->anterior(x) axis
        l_ASIS = Pel_Landmarks.landmarks["ASIS_L"]["coords"]
        r_ASIS = Pel_Landmarks.landmarks["ASIS_R"]["coords"]

        AP_vector = r_ASIS - l_ASIS
        # Rotation about y-axis that aligns asis axis with z-axis
        vec_xz = np.array([AP_vector[0], 0, AP_vector[2]])
        vec_xz = vec_xz / norm(vec_xz)
        beta = np.arccos(np.dot(vec_xz, [0, 0, 1]))
        if vec_xz[0] > 0:
            beta = -beta
            # Apply rotation about the z axis
        data, tm_temp = Transform.by_radians(data, "y", beta)
        for ld in Pel_Landmarks.landmarks:
            if ld == 'lengths' or ld == 'angles':
                continue
            Pel_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Pel_Landmarks.landmarks[ld]["coords"], tm_temp)
        tm = tm.dot(tm_temp)

        Pel_Landmarks.vertices = vertices[:]

        return tm

class FemurAxes:
    ## Align eigenvectors with global coordinate system #######################
    @staticmethod
    def align_with_coordinate_system(data):
        # Define eigenvectors
        u = Calculate.eigenvectors(data)
        com = np.mean(data, axis=0)
        
        # Define the rotation and translation 4x4 matrices
        u_t = np.transpose(u)
        r = np.array([[u_t[0,0], u_t[0,1], u_t[0,2],        0], 
                      [u_t[1,0], u_t[1,1], u_t[1,2],        0],
                      [u_t[2,0], u_t[2,1], u_t[2,2],        0],
                      [       0,        0,        0,        1]])
        
        t = np.array([[       1,        0,        0,  -com[0]],
                      [       0,        1,        0,  -com[1]],
                      [       0,        0,        1,  -com[2]],
                      [       0,        0,        0,        1]])
        
        # Define the combined transformation and apply it
        tm = np.matmul(r, t)
        data = Transform.by_tm(data, tm)
        
        # Align first eigenvector with y axis instead of x
        data, tm_temp = Transform.by_degrees(data, "z", 90)
        tm = tm_temp.dot(tm)
        
        return tm

    @staticmethod
    def set_origin(data, Fem_Landmarks):
        '''
        sets origin to condylar midpoint of femur
        '''
        vertices = Fem_Landmarks.vertices[:]
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        Fem_Landmarks.vertices = data
        # move origin
        origin = Fem_Landmarks.landmarks["con_mid"]["coords"]
        data, tm_1 = Transform.by_translation(data, -origin)
        Fem_Landmarks.landmarks["con_mid"]["coords"] = [0,0,0]
        for ld in Fem_Landmarks.landmarks:
            if ld == "con_mid" or ld == 'lengths' or ld == 'angles':
                continue
            Fem_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Fem_Landmarks.landmarks[ld]["coords"],tm_1)
        tm = tm_1.dot(tm)
        Fem_Landmarks.vertices = vertices[:]
        return tm

    @staticmethod
    def align_long_axis_z(data, Fem_Landmarks):
        vertices = Fem_Landmarks.vertices[:]
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        Fem_Landmarks.vertices = data
        head_centre = Fem_Landmarks.landmarks["head_centre"]["coords"]
        vector = (head_centre - [0,0,0])
        vec_xy = np.array([vector[0], vector[1], 0])
        vec_xy = vec_xy / norm(vec_xy)
        alpha = np.arccos(np.dot(vec_xy, [0, 1, 0]))
        if vec_xy[0] < 0:
            alpha = -alpha
        data, tm_temp = Transform.by_radians(data, "z", alpha)
        tm = tm_temp.dot(tm)
        for ld in Fem_Landmarks.landmarks:
            if ld == "con_mid" or ld == 'lengths' or ld == 'angles':
                continue
            Fem_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Fem_Landmarks.landmarks[ld]["coords"], tm_temp)
        Fem_Landmarks.vertices = vertices[:]
        return tm

    @staticmethod
    def align_long_axis_x(data, Fem_Landmarks):
        vertices = Fem_Landmarks.vertices[:]
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])

        Fem_Landmarks.vertices = data
        head_centre = Fem_Landmarks.landmarks["head_centre"]["coords"]
        vector = (head_centre - [0,0,0])
        # Find x-axis rotation that makes the vector have a x value of 0
        vec_yz = np.array([0, vector[1], vector[2]])
        vec_yz = vec_yz / norm(vec_yz)
        beta = np.arccos(np.dot(vec_yz, [0, 1, 0]))
        if vec_yz[2] > 0:
            beta = -beta
        data, tm_temp = Transform.by_radians(data, "x", beta)
        tm = tm_temp.dot(tm)
        for ld in Fem_Landmarks.landmarks:
            if ld == "con_mid" or ld == 'lengths' or ld == 'angles':
                continue
            Fem_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Fem_Landmarks.landmarks[ld]["coords"], tm_temp)
        Fem_Landmarks.vertices = vertices[:]
        return tm


    @staticmethod
    def align_z_axis(data, Fem_Landmarks):
        '''
        align the z axis of the coordinate system to the medial and lateral condylar axis in the y plane
        '''

        lat_con = Fem_Landmarks.landmarks["lat_con"]["coords"]
        med_con = Fem_Landmarks.landmarks["med_con"]["coords"]
        # Define medial->lateral (z) axis (lat to med if left femur)
        vertices = Fem_Landmarks.vertices[:]
        tm = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
        Fem_Landmarks.vertices = data
        ML_vector = med_con - lat_con
    
        # Rotation about y-axis that aligns mediolateral axis with z-axis
        vec_xz = np.array([ML_vector[0], 0, ML_vector[2]])
        vec_xz = vec_xz/norm(vec_xz)
        gamma = np.arccos(np.dot(vec_xz, [0,0,1]))
        if vec_xz[0] > 0:
            gamma = -gamma        
        
        # Apply rotation about the long axis(y)
        data, tm_temp = Transform.by_radians(data, "y", gamma)
        tm = tm_temp.dot(tm)
        for ld in Fem_Landmarks.landmarks:
            if ld == "con_mid" or ld == 'lengths' or ld == 'angles':
                continue
            Fem_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Fem_Landmarks.landmarks[ld]["coords"], tm_temp)
        Fem_Landmarks.vertices = vertices[:]
        return tm

class TibFibAxes:
    ## Align eigenvectors with global coordinate system #######################
    @staticmethod
    def align_with_coordinate_system(data):
        # Define eigenvectors
        u = Calculate.eigenvectors(data)
        com = np.mean(data, axis=0)
        
        # Define the rotation and translation 4x4 matrices
        u_t = np.transpose(u)
        r = np.array([[u_t[0,0], u_t[0,1], u_t[0,2],        0], 
                      [u_t[1,0], u_t[1,1], u_t[1,2],        0],
                      [u_t[2,0], u_t[2,1], u_t[2,2],        0],
                      [       0,        0,        0,        1]])
        
        t = np.array([[       1,        0,        0,  -com[0]],
                      [       0,        1,        0,  -com[1]],
                      [       0,        0,        1,  -com[2]],
                      [       0,        0,        0,        1]])
        
        # Define the combined transformation and apply it
        tm = np.matmul(r, t)
        data = Transform.by_tm(data, tm)
        
        # Align first eigenvector with y axis instead of z
        data, tm_temp = Transform.by_degrees(data, "z", 90) 
        tm = tm_temp.dot(tm)

        return tm

    @staticmethod
    def set_origin(data, Tibfib_Landmarks):
        '''
        sets origin to condylar midpoint of femur
        '''
        vertices = Tibfib_Landmarks.vertices[:]
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        Tibfib_Landmarks.vertices = data
        # move origin
        origin = Tibfib_Landmarks.landmarks["int_mal"]["coords"]
        data, tm_1 = Transform.by_translation(data, -origin)
        tm = tm_1.dot(tm)
        Tibfib_Landmarks.landmarks["int_mal"]["coords"] = [0, 0, 0]
        for ld in Tibfib_Landmarks.landmarks:
            if ld == "int_mal" or ld == 'lengths' or ld == 'angles':
                continue
            Tibfib_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Tibfib_Landmarks.landmarks[ld]["coords"], tm_1)

        Tibfib_Landmarks.vertices = vertices[:]
        return tm

    ## Define origin and align the mediolateral axis(z axis) ##################
    @staticmethod
    def align_ML_axis(data, Tibfib_Landmarks):
        # Use the MM and LM to align the z axis
        # Define medial->lateral (z) axis (lat to med if left femur) 
        vertices = Tibfib_Landmarks.vertices[:]
        tm = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])

        Tibfib_Landmarks.vertices = data
        lat_malleolus = Tibfib_Landmarks.landmarks["lateral_malleolus"]["coords"]
        med_malleolus = Tibfib_Landmarks.landmarks["medial_malleolus"]["coords"]

        ML_vector = med_malleolus - lat_malleolus

        # Rotation about y-axis that aligns mediolateral axis with z-axis
        vec_xz = np.array([ML_vector[0], 0, ML_vector[2]])
        vec_xz = vec_xz/norm(vec_xz)
        gamma = np.arccos(np.dot(vec_xz, [0,0,1]))
        if vec_xz[0] > 0:
            gamma = -gamma

            # Apply rotation about the long axis(y)
        data, tm_temp = Transform.by_radians(data, "y", gamma)
        tm = tm_temp.dot(tm)

        for ld in Tibfib_Landmarks.landmarks:
            if ld == "int_mal" or ld == 'lengths' or ld == 'angles':
                continue
            Tibfib_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Tibfib_Landmarks.landmarks[ld]["coords"], tm_temp)

        Tibfib_Landmarks.vertices = vertices[:]

        return tm

    @staticmethod
    def align_long_axis_z(data, Tibfib_Landmarks):
    #align y axis between the midpoint of the malleoli and the midpoint of the condyles i.e. up the shaft of the tibia
        # Create unit vector through COM of the two points
        vertices = Tibfib_Landmarks.vertices[:]
        tm = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
        Tibfib_Landmarks.vertices = data
        int_con = Tibfib_Landmarks.landmarks["con_mid"]["coords"]
        vector = (int_con - [0,0,0])
        vec_xy = np.array([vector[0],vector[1],0])
        vec_xy = vec_xy/norm(vec_xy)
        alpha = np.arccos(np.dot(vec_xy, [0,1,0]))
        if vec_xy[0] < 0:
            alpha = -alpha
        data, tm_temp = Transform.by_radians(data, "z", alpha)
        tm = tm_temp.dot(tm)

        for ld in Tibfib_Landmarks.landmarks:
            if ld == "int_mal" or ld == 'lengths' or ld == 'angles':
                continue
            Tibfib_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Tibfib_Landmarks.landmarks[ld]["coords"], tm_temp)
        Tibfib_Landmarks.vertices = vertices[:]
        return tm

    @staticmethod
    def align_long_axis_x(data,Tibfib_Landmarks):
        vertices = Tibfib_Landmarks.vertices[:]
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        Tibfib_Landmarks.vertices = data
        int_con = Tibfib_Landmarks.landmarks["con_mid"]["coords"]
        vector = (int_con - [0,0,0])
        #Find x-axis rotation that makes the vector have a z value of 0
        vec_yz = np.array([0,vector[1],vector[2]])
        vec_yz = vec_yz/norm(vec_yz)
        beta = np.arccos(np.dot(vec_yz, [0,1,0]))
        if vec_yz[2] > 0:
            beta = -beta
        data, tm_temp = Transform.by_radians(data, "x", beta)
        tm = tm_temp.dot(tm)

        for ld in Tibfib_Landmarks.landmarks:
            if ld == "int_mal" or ld == 'lengths' or ld == 'angles':
                continue
            Tibfib_Landmarks.landmarks[ld]["coords"] = Transform.by_tm(Tibfib_Landmarks.landmarks[ld]["coords"],
                                                                       tm_temp)
        
        Tibfib_Landmarks.vertices = vertices[:]
             
        return tm
