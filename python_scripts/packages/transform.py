# -*- coding: utf-8 -*-
"""
Created by Desney Greybe, 2018. 
Modified by Alicia Lie

"""

import numpy as np

from numpy.linalg import inv, norm

from packages.calculate import Calculate

#%%#########################################################################%%#
class Transform:    
    ## Apply a transformation matrix to an mx3 data array #####################
    @staticmethod
    def by_tm(data, tm):
        if len(np.shape(data)) == 1:
            data = np.append(data, 1)
            data = np.transpose(data)
            data = np.matmul(tm, data)
            data = np.transpose(data)
            data = data[:3]
        else:
            ones = np.ones((np.shape(data)[0],1))
            data = np.append(data, ones, axis=1)
            data = np.transpose(data)
            data = np.matmul(tm, data)
            data = np.transpose(data)    
            data = data[:,0:3]
        
        return data
    
    ## Apply rotation about a coordinate system axis ##########################
    @staticmethod
    def by_degrees(data, axis, deg):
        # Define rotation matrix
        theta = np.radians(deg)
        cos_theta, sin_theta = np.cos(theta), np.sin(theta)
        if axis == "x":
            tm = np.array([[1,0,0,0],[0,cos_theta,-sin_theta,0],[0,sin_theta,cos_theta,0],[0,0,0,1]])
        if axis == "y":
            tm = np.array([[cos_theta,0,sin_theta,0],[0,1,0,0],[-sin_theta,0,cos_theta,0],[0,0,0,1]])
        if axis == "z":
            tm = np.array([[cos_theta,-sin_theta,0,0],[sin_theta,cos_theta,0,0],[0,0,1,0],[0,0,0,1]])
        
        # Apply rotation
        data = Transform.by_tm(data, tm)
                
        return data, tm

    @staticmethod
    def by_radians(data, axis, rad):
        # Define rotation matrix
        theta = rad
        cos_theta, sin_theta = np.cos(theta), np.sin(theta)
        if axis == "x":
            tm = np.array([[1,0,0,0],[0,cos_theta,-sin_theta,0],[0,sin_theta,cos_theta,0],[0,0,0,1]])
        if axis == "y":
            tm = np.array([[cos_theta,0,sin_theta,0],[0,1,0,0],[-sin_theta,0,cos_theta,0],[0,0,0,1]])
        if axis == "z":
            tm = np.array([[cos_theta,-sin_theta,0,0],[sin_theta,cos_theta,0,0],[0,0,1,0],[0,0,0,1]])
        
        # Apply rotation
        data = Transform.by_tm(data, tm)
                
        return data, tm

    ## Apply rotation (in degrees) about a coordinate system axis #############
    @staticmethod
    def by_translation(data, t):        
        # Define translation matrix
        tm = np.array([[1,0,0,t[0]],[0,1,0,t[1]],[0,0,1,t[2]],[0,0,0,1]])

        # Apply translation
        data = Transform.by_tm(data, tm)
                
        return data, tm

    ## Reflect mesh in the specified plane ####################################
    @staticmethod
    def by_reflection(data, plane):        
        # Define transformation matrix        
        if plane == "xy":    # z-axis
            tm = np.array([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,1]])
        if plane == "xz":    # y-axis
            tm = np.array([[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,1]])
        if plane == "yz":    # x-axis
            tm = np.array([[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    
        # Apply transformation
        data = Transform.by_tm(data, tm)

        return data, tm

    @staticmethod
    def reflect(self, plane):
        try:
            self.vertices, tm = Transform.by_reflection(self.vertices, plane)

        except:
            #            print "~ Unique vertices not yet identified ~"
            pass

        try:
            # as reflecting about x plane, swap v1 and v2 in order to have faces and normals facing outwards
            v0 = self.v0
            v1 = self.v1
            self.v0, tm = Transform.by_reflection(v1, plane)
            self.v1, tm = Transform.by_reflection(v0, plane)
            self.v2, tm = Transform.by_reflection(self.v2, plane)
            self.normals, tm = Transform.by_reflection(self.normals, plane)

        except:
            pass

        return tm
    ## Scale data by isometrically using a scale factor #######################
    @staticmethod
    def by_scale_factor(data, scale_factor):
        # Define transformation matrix        
        tm = np.array([[scale_factor,0,0,0],[0,scale_factor,0,0],[0,0,scale_factor,0],[0,0,0,1]]) 

        # Apply transformation
        data = Transform.by_tm(data, tm)
                
        return data, tm

    ## Rotate data about an arbitrary axis ####################################
    @staticmethod
    def about_vector(data, point, vector, rad):
        # Align vector with z-axis
        tm_1 = Align.vector_with_axis(point, vector, "z")
            
        # Perform transformations
        data = Transform.by_tm(data, tm_1)
        data, tm_2 = Transform.by_radians(data, "z", rad)
        data = Transform.by_tm(data, inv(tm_1))
    
        # Calculate combined transformation        
        tm = inv(tm_1).dot(tm_2).dot(tm_1)
        
        return data, tm
        
#%%#########################################################################%%#
class Align:    
    ## Apply a transformation matrix to an mx3 data array #####################
    @staticmethod
    def with_coordinate_system(data, centre_method):
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
        
        # Align first eigenvector with z-axis (instead of x)
        data, tm_temp = Transform.by_degrees(data, "y", 90)
        tm = tm_temp.dot(tm)
        
        # Correct alignment in z-axis
        if centre_method != "com":
            prox = np.amin(data, axis=0)[2]
            dist = np.amax(data, axis=0)[2]        
            if centre_method == "proximal":
                t = np.array([0,0,-prox])
            if centre_method == "distal":
                t = np.array([0,0,-dist])
            if centre_method == "midpoint":
                mid = prox + (dist-prox)/2
                t = np.array([0,0,-mid])
            else:
                print ("!! Alignment method is incorrect !!")
            
            data, tm_temp = Transform.by_translation(data, t)
            tm = tm_temp.dot(tm)
        
        return data, tm
    
    ## Transformation that aligns vector with coordinate system axis ##########    
    @staticmethod
    def vector_with_axis(point, vector, axis):
        # Create translation matrix
        #point, tm_1 = Transform.by_translation(point, -point)
        
        # Find y-axis rotation to bring long axis into zy-plane
        # vec_xz = np.array([vector[0],0,vector[2]])
        # vec_xz = vec_xz/norm(vec_xz)
        # alpha = np.arccos(np.dot(vec_xz, [0,0,1]))
        # if vec_xz[0] > 0:
        #     alpha = -alpha
        # vector, tm_2 = Transform.by_radians(vector, "y", alpha)

        # Find z-axis rotation to make the vector have a x value of 0
        vec_xy = np.array([vector[0],vector[1],0])
        vec_xy = vec_xy/norm(vec_xy)
        alpha = np.arccos(np.dot(vec_xy, [0,1,0]))
        if vec_xy[0] > 0:
            alpha = -alpha
        vector, tm_2 = Transform.by_radians(vector, "z", alpha)
        
        # Find x-axis rotation that makes the vector have a z value of 0 #OH DIDNT LOOK AT THIS
        vec_yz = np.array([0,vector[1],vector[2]])
        vec_yz = vec_yz/norm(vec_yz)
        beta = np.arccos(np.dot(vec_yz, [0,1,0]))
        if vec_yz[2] > 0:
            beta = -beta
        vector, tm_3 = Transform.by_radians(vector, "x", beta)

        # Combine tranformations
        tm = tm_3.dot(tm_2)#.dot(tm_1)

        # If not aligned with z-axis
        # if axis != "y":
        #     if axis == "x":
        #         vector, tm_4 = Transform.by_degrees(vector, "y", 90)
        #     elif axis == "z":
        #         vector, tm_4 = Transform.by_degrees(vector, "x", -90)
        #     tm = tm_4.dot(tm)
        
        return tm

#%%#########################################################################%%#