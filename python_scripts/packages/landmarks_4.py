# -*- coding: utf-8 -*-
# landmarks_v4.py
"""
Created on Fri Sep  3 13:01:35 2021
Classes and functions for finding anatomical landmarks and bone measurements from the pelvis, femur,
and tibia/fibula. Used in find_landmarks.py

Requires lists prepared from cmgui as a result of node selection. 
Works for 'left' bones only

x = posterior --> anterior 
y = distal --> proximal 
z = medial --> lateral 

@author: laura Carman
"""
# %%#########################################################################%%#
import numpy as np
import scipy
from numpy.linalg import norm
from packages.cylinder_fit import fit as cfit
from packages.lower_bone_axes_2 import FemurAxes, PelvisAxes, TibFibAxes
from packages.cmiss import cm_lists
from gias2.common import geoprimitives as gp
from packages.transform import Transform
import scipy.linalg
import open3d as o3d


def calc_distances(p0, points):
    # Calculates the distance of each point in 3D space from the origin/p0
    squared_dist = np.sum((p0 - points) ** 2, axis=1)
    return np.sqrt(squared_dist)


def _norm(v):
    return scipy.divide(v, scipy.sqrt((scipy.array(v) ** 2.0).sum()))


def _mag(v):
    return scipy.sqrt((scipy.array(v) ** 2.0).sum())


# %%#########################################################################%%#
class PelvisLandmarks(STL):
    ## Initiate class #########################################################
    def __init__(self, file_name):
        self.file_name = file_name.split(".")[0]
        self.vertices = None
        self.vertex_lists = {}
        self.landmarks = {}

    ## Load stl mesh vertices #################################################
    def load_vertices(self, file_path):
        mesh = o3d.io.read_triangle_mesh(file_path)
        self.vertices = mesh.vertices

        return

        ## Load vertex lists for femur landmarks #################################

    def load_landmark_lists(self, lists_path):
        # Distal portions 
        self.vertex_lists["pubic_symphysis_R"] = cm_lists.node("pubic_symphysis_R", lists_path) - 1
        self.vertex_lists["pubic_symphysis_L"] = cm_lists.node("pubic_symphysis_L", lists_path) - 1
        # self.vertex_lists["pubic_tubercle_L"] = cm_lists.node("pubic_tubercle_L", lists_path) - 1
        self.vertex_lists["right_hip_joint"] = cm_lists.node("right_hip_joint", lists_path) - 1
        self.vertex_lists["left_hip_joint"] = cm_lists.node("left_hip_joint", lists_path) - 1

        # Proximal portions
        self.vertex_lists["ASIS_L"] = cm_lists.node("ASIS_L", lists_path) - 1
        self.vertex_lists["ASIS_R"] = cm_lists.node("ASIS_R", lists_path) - 1
        self.vertex_lists["PSIS_L"] = cm_lists.node("PSIS_L", lists_path) - 1
        self.vertex_lists["PSIS_R"] = cm_lists.node("PSIS_R", lists_path) - 1

        return

        ## Tranform vertices into the joint coordinate system #############

    def align_JCS(self):
        tm_1 = PelvisAxes.set_origin(self.vertices,self)
        self.transform(tm_1)
        tm_2 = PelvisAxes.align_x_axis(self.vertices, self)
        self.transform(tm_2)
        tm_4 = PelvisAxes.align_z_axis_x(self.vertices, self)
        self.transform(tm_4)
        tm_6 = PelvisAxes.align_z_axis_y(self.vertices, self)
        self.transform(tm_6)

        return tm_1.dot(tm_2).dot(tm_4).dot(tm_6)
    
    def align_com(self):
        tm_2 = PelvisAxes.align_x_axis(self.vertices, self)
        self.transform(tm_2)
        tm_4 = PelvisAxes.align_z_axis_x(self.vertices, self)
        self.transform(tm_4)
        tm_6 = PelvisAxes.align_z_axis_y(self.vertices, self)
        self.transform(tm_6)

        return tm_2.dot(tm_4).dot(tm_6)
    
    
    def align_LASIS(self):
        tm_1 = PelvisAxes.set_origin_lasis(self.vertices,self)
        self.transform(tm_1)
        tm_2 = PelvisAxes.align_x_axis(self.vertices, self)
        self.transform(tm_2)
        tm_4 = PelvisAxes.align_z_axis_x(self.vertices, self)
        self.transform(tm_4)
        tm_6 = PelvisAxes.align_z_axis_y(self.vertices, self)
        self.transform(tm_6)
        return tm_1
    
    def align_midasis(self):
        tm_1 = PelvisAxes.set_origin_midasis(self.vertices,self)
        self.transform(tm_1)
        tm_2 = PelvisAxes.align_x_axis(self.vertices, self)
        self.transform(tm_2)
        tm_4 = PelvisAxes.align_z_axis_x(self.vertices, self)
        self.transform(tm_4)
        tm_6 = PelvisAxes.align_z_axis_y(self.vertices, self)
        self.transform(tm_6)
        return tm_1

        ## Define the anterior superior iliac spine of the left and right pelvis ############################

    def ASIS(self):
        # Find the point that is the most +ve in the x direction for the left and right side
        y_centre = np.mean(self.vertices[:, 1])
        z_centre = np.mean(self.vertices[:, 2])
        prox_vertices = self.vertices[self.vertices[:, 1] >= y_centre]
        right_vertices = prox_vertices[prox_vertices[:, 2] >= z_centre]
        left_vertices = prox_vertices[prox_vertices[:, 2] <= z_centre]
        i = right_vertices.argmax(axis=0)[0]
        self.landmarks["ASIS_R"] = {}
        self.landmarks["ASIS_R"]["coords"] = right_vertices[i, :]

        # Find the point that is the most +ve in the x direction for the right side
        j = left_vertices.argmax(axis=0)[0]
        self.landmarks["ASIS_L"] = {}
        self.landmarks["ASIS_L"]["coords"] = left_vertices[j, :]

        mid_point = (self.landmarks["ASIS_R"]["coords"] + self.landmarks["ASIS_L"]["coords"]) / 2

        self.landmarks["ASIS_mid"] = {}
        self.landmarks["ASIS_mid"]["coords"] = mid_point
        return

    ## Define the posterior superior iliac spine of the left and right pelvis ############################
    def PSIS(self):
        # Find the point that is the most -ve in the x direction for the left and right side
        y_centre = np.mean(self.vertices[:, 1])
        z_centre = np.mean(self.vertices[:, 2])
        prox_vertices = self.vertices[self.vertices[:, 1] >= y_centre]
        right_vertices = prox_vertices[prox_vertices[:, 2] >= z_centre]
        left_vertices = prox_vertices[prox_vertices[:, 2] <= z_centre]
        i = right_vertices.argmin(axis=0)[0]
        self.landmarks["PSIS_R"] = {}
        self.landmarks["PSIS_R"]["coords"] = right_vertices[i, :]

        # Find the point that is the most +ve in the x direction for the right side
        i = left_vertices.argmin(axis=0)[0]
        self.landmarks["PSIS_L"] = {}
        self.landmarks["PSIS_L"]["coords"] = left_vertices[i, :]

        mid_point = (self.landmarks["PSIS_R"]["coords"] + self.landmarks["PSIS_L"]["coords"]) / 2

        self.landmarks["PSIS_mid"] = {}
        self.landmarks["PSIS_mid"]["coords"] = mid_point
        return

    ## Define the centre of the hip joint of the right pelvis ############################
    def hip_joint(self):
        # Find the point that is in the middle of the hip joint if modelled by a sphere
        # Extract relevant vertices
        hip_joint = self.vertices[self.vertex_lists["right_hip_joint"]]

        # Centre gives the coords of the centre point and rad is the radius of that sphere
        centre, rad = gp.fitSphereAnalytic(hip_joint)

        self.landmarks["right_hip_joint"] = {}
        self.landmarks["right_hip_joint"]["radius"] = rad
        self.landmarks["right_hip_joint"]["coords"] = centre

        left_hip_joint = self.vertices[self.vertex_lists["left_hip_joint"]]

        # Centre gives the coords of the centre point and rad is the radius of that sphere
        centre, rad = gp.fitSphereAnalytic(left_hip_joint)

        self.landmarks["left_hip_joint"] = {}
        self.landmarks["left_hip_joint"]["radius"] = rad
        self.landmarks["left_hip_joint"]["coords"] = centre
        return

    ## Define the pubic symphysis as well as the pubic tubercles of the pelvis ############################
    def PS(self):
        # Find the point that is most +ve in the z axis for the left tubercle
        i = self.vertices[self.vertex_lists["pubic_symphysis_L"]].argmax(axis=0)[2]
        self.landmarks["PT_L"] = {}  # Re select the nodes###
        self.landmarks["PT_L"]["ID"] = self.vertex_lists["pubic_symphysis_L"][i]
        self.landmarks["PT_L"]["coords"] = self.vertices[self.vertex_lists["pubic_symphysis_L"][i]]

        # Find the point that is most -ve in the z axis for the right pubic turbercle
        i = self.vertices[self.vertex_lists["pubic_symphysis_R"]].argmin(axis=0)[2]
        self.landmarks["PT_R"] = {}
        self.landmarks["PT_R"]["ID"] = self.vertex_lists["pubic_symphysis_R"][i]
        self.landmarks["PT_R"]["coords"] = self.vertices[self.vertex_lists["pubic_symphysis_R"][i]]

        # The pubic symphysis is the point in the middle of the two pubic turbercles 
        self.landmarks["PS"] = {}
        self.landmarks["PS"]["coords"] = (self.landmarks["PT_L"]["coords"] + self.landmarks["PT_R"]["coords"]) / 2

        return

    def length_measurements(self):
        squared_dist = np.sum((self.landmarks["ASIS_R"]["coords"] - self.landmarks["ASIS_L"]["coords"]) ** 2, axis=0)
        length = np.sqrt(squared_dist)

        self.landmarks["lengths"] = {}
        self.landmarks["lengths"]["ASIS_width"] = length

        squared_dist = np.sum((self.landmarks["PSIS_R"]["coords"] - self.landmarks["PSIS_L"]["coords"]) ** 2, axis=0)
        length = np.sqrt(squared_dist)

        self.landmarks["lengths"]["PSIS_width"] = length

        squared_dist = np.sum((self.landmarks["ASIS_mid"]["coords"] - self.landmarks["PSIS_mid"]["coords"]) ** 2,
                              axis=0)
        length = np.sqrt(squared_dist)

        self.landmarks["lengths"]["depth"] = length
        return


# %%#########################################################################%%#
class FemurLandmarks(STL):
    ## Initiate class #########################################################
    def __init__(self, file_name):
        self.file_name = file_name.split(".")[0]
        self.vertices = None
        self.vertex_lists = {}
        self.landmarks = {}

        STL.__init__(self, self.file_name)

    ## Load stl mesh vertices #################################################
    def load_vertices(self, file_path):
        self.load_stl(file_path)
        self.unique_vertices()

        return

    ## Load vertex lists for femur landmarks #################################
    def load_landmark_lists(self, lists_path):
        '''
        Load vertex list for femur landmarks, lists obtained using cm-gui
        '''
        # Distal portions
        self.vertex_lists["lateral_epicondyle"] = cm_lists.node("lateral_epicondyle", lists_path) - 1
        self.vertex_lists["medial_epicondyle"] = cm_lists.node("medial_epicondyle", lists_path) - 1
        self.vertex_lists["lateral_condyle"] = cm_lists.node("lateral_condyle", lists_path) - 1
        self.vertex_lists["medial_condyle"] = cm_lists.node("medial_condyle", lists_path) - 1
        self.vertex_lists["femur_int_notch"] = cm_lists.node("femur_int_notch", lists_path) - 1
        self.vertex_lists["femur_pat_surf"] = cm_lists.node("femur_pat_surf", lists_path) - 1
        self.vertex_lists["lat_bot_con"] = cm_lists.node("lat_bot_con", lists_path) - 1
        self.vertex_lists["med_bot_con"] = cm_lists.node("med_bot_con", lists_path) - 1
        self.vertex_lists["left_condyle"] = cm_lists.node("left_condyle", lists_path) - 1
        self.vertex_lists["right_condyle"] = cm_lists.node("right_condyle", lists_path) - 1
        self.vertex_lists["lat_bot_con"] = cm_lists.node("lat_bot_con", lists_path) - 1
        self.vertex_lists["med_bot_con"] = cm_lists.node("med_bot_con", lists_path) - 1

        # Proximal portions
        self.vertex_lists["femur_head"] = cm_lists.node("femur_head", lists_path) - 1
        # self.vertex_lists["greater_trochanter"]= cm_lists.node("greater_trochanter", lists_path) - 1
        self.vertex_lists["femur_neck"] = cm_lists.node("femur_neck", lists_path) - 1
        self.vertex_lists["femur_prox_shaft"] = cm_lists.node("femur_prox_shaft", lists_path) - 1
        self.vertex_lists["fem_bot"] = cm_lists.node("fem_bot", lists_path) - 1
        # self.vertex_lists["femur_top"] = cm_lists.node("femur_top", lists_path) - 1

        return

    def align_angles(self):
        '''
        align femur to ISB recommendations with origin at condylar midpoint
        y axis between condylar midpoint and head centre
        z axis between condyles  
        x axis cross product
        '''
        tm_1 = FemurAxes.set_origin(self.vertices,self)
        self.transform(tm_1)
        tm_2 = FemurAxes.align_z_axis(self.vertices, self)
        self.transform(tm_2)
        tm_4 = FemurAxes.align_long_axis_z(self.vertices, self)
        self.transform(tm_4)
        tm_6 = FemurAxes.align_long_axis_x(self.vertices, self)
        self.transform(tm_6)

        return tm_1.dot(tm_2).dot(tm_4).dot(tm_6)

    def head_centre(self):
        '''
        Find the head and neck centre of the femur by a sphere fit to the femoral head.
        
        '''
        # Extract relevant vertices
        femur_head = self.vertices[self.vertex_lists["femur_head"]]

        # centre gives the coords of the centre point and rad is the radius of that sphere
        centre, rad = gp.fitSphereAnalytic(femur_head)
        self.landmarks["head_centre"] = {}
        self.landmarks["head_centre"]["coords"] = centre
        self.landmarks["head_centre"]["radius"] = rad
        return

    def neck_centre(self):
        '''
        find the centre of the femoral neck and the centre of the proximal femoral
        shaft using cylinder fit
        '''
        # neck centre
        neck_cen = self.vertices[self.vertex_lists["femur_neck"]]

        # centre gives the coords of the centre point and rad is the radius of that sphere
        W, neck_centre, rad, err = cfit(neck_cen)
        self.landmarks["neck_centre"] = {}
        self.landmarks["neck_centre"]["coords"] = neck_centre
        self.landmarks["neck_centre"]["radius"] = rad

        # shaft centre
        # find the centre of the femoral neck
        shaft_cen = self.vertices[self.vertex_lists["femur_prox_shaft"]]

        # centre gives the coords of the centre point and rad is the radius of that sphere
        W, shaft_centre, rad, err = cfit(shaft_cen)
        # neck_centre = self.vertices[self.vertex_lists["femur_neck"]].mean(axis=0)
        self.landmarks["shaft_centre"] = {}
        self.landmarks["shaft_centre"]["coords"] = shaft_centre
        return

    def great_trochant(self):
        '''
        Define the greater trochanter of the femur by most dorsal (-x) and most 
        proximal (+y) by calculating distance from origin. Find the most lateral 
        point of the greater trochanter (furthest point in both -z and x from origin)
        '''
        y_centre = np.mean(self.vertices[:, 1])
        z_centre = np.mean(self.vertices[:, 2])
        threshold = (np.amax(self.vertices[:, 2]) - np.min(self.vertices[:, 2])) * 0.2
        prox_vertices = self.vertices[self.vertices[:, 1] >= y_centre]
        prox_lat_vertices = prox_vertices[prox_vertices[:, 2] <= z_centre - threshold]

        # greater trochanter proximal
        dist = 0
        for n in range(len(prox_lat_vertices)):
            # x = prox_lat_vertices[n,0]
            y = prox_lat_vertices[n, 1]
            z = prox_lat_vertices[n, 2]
            length = np.sqrt(y ** 2 + z ** 2)
            if length > dist:
                dist = length
                i = n
                # i = prox_lat_vertices.argmax(axis=0)[1]
        self.landmarks["great_trochant"] = {}
        # self.landmarks["great_trochant"]["ID"] = prox_lat_vertices[i]
        self.landmarks["great_trochant"]["coords"] = prox_lat_vertices[i, :]

        # greater trochanter lateral
        dist = 0
        for n in range(len(prox_lat_vertices)):
            x = prox_lat_vertices[n, 0]
            z = prox_lat_vertices[n, 2]
            length = np.sqrt(0.5 * x ** 2 + z ** 2)
            if length > dist:
                dist = length
                i = n

                # i = prox_lat_vertices.argmin(axis=0)[2]
        self.landmarks["lat_great_trochant"] = {}
        # self.landmarks["lat_great_trochant"]["ID"] = prox_lat_vertices[i]
        self.landmarks["lat_great_trochant"]["coords"] = prox_lat_vertices[i, :]
        return

    def epicondyles(self):
        '''
        align z axis to epicondyles to find most medial and lateral points (true epicondyles),
        transformation disregarded. Find the epicondylar width
        '''
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        gamma = None
        c = 0
        self.landmarks["lat_epicon"] = {}
        self.landmarks["med_epicon"] = {}
        while gamma != 0 and c < 10:
            y_centre = np.mean(self.vertices[:, 1])
            dist_vertices = self.vertices[self.vertices[:, 1] <= y_centre]
            # data = self.vertices
            i = dist_vertices.argmin(axis=0)[2]
            j = dist_vertices.argmax(axis=0)[2]

            lat_epi = dist_vertices[i,:]
            med_epi = dist_vertices[j,:]

            self.landmarks["lat_epicon"]["coords"] = lat_epi
            self.landmarks["med_epicon"]["coords"] = med_epi

            ML_vector = med_epi - lat_epi

            # Rotation about y-axis that aligns mediolateral axis with z-axis
            vec_xz = np.array([ML_vector[0], 0, ML_vector[2]])
            vec_xz = vec_xz / norm(vec_xz)
            gamma = np.arccos(np.dot(vec_xz, [0, 0, 1]))
            if vec_xz[0] > 0:
                gamma = -gamma

                # Apply rotation about the long axis(y)
            self.vertices, tm_temp = Transform.by_radians(self.vertices, "y", gamma)
            tm = tm_temp.dot(tm)
            c = c + 1
        tm_inv = np.linalg.inv(tm)
        self.vertices = Transform.by_tm(self.vertices, tm_inv)
        self.landmarks["lat_epicon"]["coords"] = Transform.by_tm(self.landmarks["lat_epicon"]["coords"],
                                                                        tm_inv)
        self.landmarks["med_epicon"]["coords"] = Transform.by_tm(self.landmarks["med_epicon"]["coords"],
                                                                       tm_inv)

        # epicondylar width
        squared_dist = np.sum((self.landmarks["med_epicon"]["coords"] - self.landmarks["lat_epicon"]["coords"]) ** 2,
                              axis=0)
        width = np.sqrt(squared_dist)

        self.landmarks["lengths"] = {}
        self.landmarks["lengths"]["epicon_width"] = width

        # Find the epicondylar midpoint
        mid_point = (self.landmarks["med_epicon"]["coords"] + self.landmarks["lat_epicon"]["coords"]) / 2

        self.landmarks["epicon_mid"] = {}
        self.landmarks["epicon_mid"]["coords"] = mid_point
        return

    def femoral_length(self):
        '''
        calculate the length of the femur from the greater trochanter landmark 
        to the lateral epicondyle landmark
        '''
        squared_dist = np.sum(
            (self.landmarks["great_trochant"]["coords"] - self.landmarks["lat_epicon"]["coords"]) ** 2, axis=0)
        length = np.sqrt(squared_dist)

        self.landmarks["lengths"]["femoral_length"] = length
        return

    def posterior_condyles(self):
        '''
        find posterior condyles by aligning them to the z axis to find the most posterior points,
        disregards the transformation
        '''
        tm = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
        gamma = None
        c = 0
        self.landmarks["lat_con_post"] = {}
        self.landmarks["med_con_post"] = {}
        while gamma != 0 and c < 10:
            #data = self.vertices
            y_centre = np.mean(self.vertices[:,1])
            z_centre = np.mean(self.vertices[:,2])
            dist_vertices = self.vertices[self.vertices[:,1] <= y_centre]
            dist_lat_vertices = dist_vertices[dist_vertices[:,2]<=z_centre]
            dist_med_vertices = dist_vertices[dist_vertices[:,2]>=z_centre]
            i = dist_lat_vertices.argmin(axis=0)[0]
            j = dist_med_vertices.argmin(axis=0)[0]
            lat_con = dist_lat_vertices[i,:]

            #self.landmarks["lat_con_post"]["ID"] = self.vertex_lists["lateral_condyle"][i]
            self.landmarks["lat_con_post"]["coords"] = dist_lat_vertices[i,:]
            med_con = dist_med_vertices[j,:]

            #self.landmarks["med_con_post"]["ID"] = self.vertex_lists["medial_condyle"][j]
            self.landmarks["med_con_post"]["coords"] = dist_med_vertices[j,:]
            ML_vector = med_con - lat_con

            # Rotation about y-axis that aligns mediolateral axis with z-axis
            vec_xz = np.array([ML_vector[0], 0, ML_vector[2]])
            vec_xz = vec_xz/norm(vec_xz)
            gamma = np.arccos(np.dot(vec_xz, [0,0,1]))
            if vec_xz[0] > 0:
                gamma = -gamma

                # Apply rotation about the long axis(y)
            self.vertices, tm_temp = Transform.by_radians(self.vertices, "y", gamma)
            tm = tm_temp.dot(tm)
            c = c + 1
        tm_inv = np.linalg.inv(tm)
        self.vertices = Transform.by_tm(self.vertices, tm_inv)
        self.landmarks["lat_con_post"]["coords"]= Transform.by_tm(self.landmarks["lat_con_post"]["coords"], tm_inv)
        self.landmarks["med_con_post"]["coords"]= Transform.by_tm(self.landmarks["med_con_post"]["coords"], tm_inv)
        return

    def condyles(self):
        '''
        define the condyles of the femur using a cylinder fit to each condyle
        '''
        # fit cylinder to condyles
        lat_con = self.vertices[self.vertex_lists["left_condyle"]]

        # centre gives the coords of the centre point and rad is the radius of that cylinder
        W, centre, rad, err = cfit(lat_con)
        self.landmarks["lat_con"] = {}
        self.landmarks["lat_con"]["coords"] = centre
        self.landmarks["lat_con"]["radius"] = rad

        # fit cylinder to condyles
        med_con = self.vertices[self.vertex_lists["right_condyle"]]

        # centre gives the coords of the centre point and rad is the radius of that cylinder
        W, centre, rad, err = cfit(med_con)

        self.landmarks["med_con"] = {}
        self.landmarks["med_con"]["coords"] = centre
        self.landmarks["med_con"]["radius"] = rad

        # find condyle midpoint
        mid_point = (self.landmarks["med_con"]["coords"] + self.landmarks["lat_con"]["coords"]) / 2

        self.landmarks["con_mid"] = {}
        self.landmarks["con_mid"]["coords"] = mid_point
        return

    def bot_condyles(self):
        '''
        find the most distal points of the medial and lateral condyles of the femur
        '''
        y_centre = np.mean(self.vertices[:,1])
        z_centre = np.mean(self.vertices[:,2])
        threshold = (np.amax(self.vertices[:,2]) - np.min(self.vertices[:,2])) * 0.2
        dist_vertices = self.vertices[self.vertices[:,1] <= y_centre]
        dist_lat_vertices = dist_vertices[dist_vertices[:,2]<=z_centre-threshold]
        dist_med_vertices = dist_vertices[dist_vertices[:,2]>=z_centre-threshold]
        i = dist_lat_vertices.argmin(axis=0)[1]
        j = dist_med_vertices.argmin(axis=0)[1]
        #lateral condyle
        #i = self.vertices[self.vertex_lists["lat_bot_con"]].argmin(axis=0)[1]
        self.landmarks["lat_bot_con"] = {}
        #self.landmarks["lat_bot_con"]["ID"] = self.vertex_lists["lat_bot_con"][i]
        #self.landmarks["lat_bot_con"]["coords"] = self.vertices[self.vertex_lists["lat_bot_con"][i]]
        self.landmarks["lat_bot_con"]["coords"] = dist_lat_vertices[i,:]

        #medial condyle
        #i = self.vertices[self.vertex_lists["med_bot_con"]].argmin(axis=0)[1]
        self.landmarks["med_bot_con"] = {}
        #self.landmarks["med_bot_con"]["ID"] = self.vertex_lists["med_bot_con"][i]
        #self.landmarks["med_bot_con"]["coords"] = self.vertices[self.vertex_lists["med_bot_con"][i]]
        self.landmarks["med_bot_con"]["coords"] = dist_med_vertices[j,:]

        return

    def int_fossa(self):
        '''
        finds the patellar surface (most anterior node of the femoral notch) and 
        the intercondylar fossa (most posterior vertex of the int notch)
        '''
        # patellar surface
        i = self.vertices[self.vertex_lists["femur_int_notch"]].argmax(axis=0)[0]
        self.landmarks["pat_surf"] = {}
        self.landmarks["pat_surf"]["ID"] = self.vertex_lists["femur_int_notch"][i]
        self.landmarks["pat_surf"]["coords"] = self.vertices[self.vertex_lists["femur_int_notch"][i]]
        # intercondylar fossa
        i = self.vertices[self.vertex_lists["femur_int_notch"]].argmin(axis=0)[0]
        self.landmarks["int_fossa"] = {}
        self.landmarks["int_fossa"]["ID"] = self.vertex_lists["femur_int_notch"][i]
        self.landmarks["int_fossa"]["coords"] = self.vertices[self.vertex_lists["femur_int_notch"][i]]
        # notch midpoint
        mid_point = (self.vertices[self.landmarks["pat_surf"]["ID"]] + self.vertices[
            self.landmarks["int_fossa"]["ID"]]) / 2
        self.landmarks["notch_mid"] = {}
        self.landmarks["notch_mid"]["coords"] = mid_point
        return

    def anteversion_angle(self):
        '''
        Calculation of the anteversion in 2d looking down the y axis. Uses the
        landmarks of the posterior condyles, head centre, and neck centre.
        '''
        # 2d anteversion angle
        lcp = self.landmarks['lat_con_post']['coords']
        mcp = self.landmarks['med_con_post']['coords']
        lc = self.landmarks['lat_epicon']['coords']
        mc = self.landmarks['med_epicon']['coords']
        lcc = self.landmarks['lat_con']['coords']
        mcc = self.landmarks['med_con']['coords']
        hjc = self.landmarks['head_centre']['coords']
        nc = self.landmarks['neck_centre']['coords']
        post_vec = _norm(mcp - lcp)
        epi_vec = _norm(mc-lc)
        con_vec = _norm(mcc-lcc)
        if post_vec[2] < 0:
            post_vec = post_vec * -1
        if epi_vec[2] < 0:
            epi_vec = epi_vec * -1
        if con_vec[2] < 0:
            con_vec = con_vec * -1
        neck_vec = _norm(hjc - nc)
        if neck_vec[2] < 0:
            neck_vec = neck_vec * -1
        # calculate anteversion angle
        angle = np.rad2deg(np.arctan2(post_vec[2] * neck_vec[0] - post_vec[0] * neck_vec[2],
                                      post_vec[0] * neck_vec[0] + post_vec[2] * neck_vec[2]))
        self.landmarks["angles"] = {}
        self.landmarks["angles"]["AA_2d"] = angle

        angle2 = np.rad2deg(np.arctan2(epi_vec[2] * neck_vec[0] - epi_vec[0] * neck_vec[2],
                                      epi_vec[0] * neck_vec[0] + epi_vec[2] * neck_vec[2]))
        self.landmarks['angles']['AA_2d_epi'] = angle2

        angle3 = np.rad2deg(np.arctan2(con_vec[2] * neck_vec[0] - con_vec[0] * neck_vec[2],
                                      con_vec[0] * neck_vec[0] + con_vec[2] * neck_vec[2]))
        self.landmarks['angles']['AA_2d_con'] = angle3

        # fit plane to bottom nodes
        p = gp.fitPlaneLS(self.vertices[self.vertex_lists["fem_bot"]])
        if abs(p.X[2]) > abs(p.Y[2]):
            vec1 = p.X
        else:
            vec1 = p.Y
        if vec1[2] < 0:
            vec1 = vec1 * -1
        angle4 = np.rad2deg(np.arctan2(vec1[2] * neck_vec[0] - vec1[0] * neck_vec[2],
                                      vec1[0] * neck_vec[0] + vec1[2] * neck_vec[2]))
        self.landmarks['angles']['AA_2d_plane'] = angle4

        # 3d anteversion angle
        angle_3d = np.rad2deg(np.arccos(np.dot(neck_vec, post_vec) / (_mag(neck_vec) * _mag(post_vec))))
        self.landmarks['angles']['AA_3d'] = angle_3d
        return

    def neck_shaft_angle(self):
        '''
        calculation of neck shaft angle in 2d looking down the x axis. Uses the
        landmarks of the head centre, neck centre, proximal shaft centre, and 
        condylar midpoint.
        '''
        con_mid = self.landmarks['con_mid']['coords']
        shaft_mid = self.landmarks['shaft_centre']['coords']
        hjc = self.landmarks['head_centre']['coords']
        nc = self.landmarks['neck_centre']['coords']
        shaft_vec = _norm(con_mid - shaft_mid)
        if shaft_vec[1] > 0:
            shaft_vec = shaft_vec * -1
        neck_vec = _norm(hjc - nc)
        if neck_vec[1] < 0:
            neck_vec = neck_vec * -1
        angle = np.rad2deg(np.arctan2(shaft_vec[2] * neck_vec[1] - shaft_vec[1] * neck_vec[2],
                                      shaft_vec[1] * neck_vec[1] + shaft_vec[2] * neck_vec[2]))
        self.landmarks["angles"]["NSA_2d"] = angle

        # 3d neck shaft angle
        angle_3d = np.rad2deg(np.arccos(np.dot(neck_vec, shaft_vec) / (_mag(neck_vec) * _mag(shaft_vec))))
        self.landmarks['angles']['NSA_3d'] = angle_3d
        return

    def mLDFA(self):
        '''
        calculation of mechanical lateral distal femoral angle as the angle between
        the y axis and the bottom axis of the femur.
        '''
        lbc = self.landmarks['lat_bot_con']['coords']
        mbc = self.landmarks['med_bot_con']['coords']
        bot_vec = _norm(lbc - mbc)
        if bot_vec[2]>0:
            bot_vec = bot_vec*-1
        angle = np.rad2deg(np.arctan2(0 * bot_vec[1] - 1 * bot_vec[2], 1 * bot_vec[1] + 0 * bot_vec[2]))
        self.landmarks['angles']['mLDFA_2d'] = angle

        # 3d mLDFA
        angle_3d = np.rad2deg(np.arccos(np.dot(bot_vec, [0, 1, 0]) / (_mag(bot_vec) * _mag([0, 1, 0]))))
        self.landmarks['angles']['mLDFA_3d'] = angle_3d

        # fit plane to bottom nodes
        p = gp.fitPlaneLS(self.vertices[self.vertex_lists["fem_bot"]])
        if abs(p.X[2]) > abs(p.Y[2]):
            vec1 = p.X
        else:
            vec1 = p.Y
        if vec1[2] > 0:
            vec1 = vec1 * -1
        angle = np.rad2deg(np.arctan2(0 * vec1[1] - 1 * vec1[2], 1 * vec1[1] + 0 * vec1[2]))
        self.landmarks['angles']['mLDFA_plane'] = angle

        return

    def bicondylar_angle(self):
        '''
        calculation of the biconylar angle as the angle between the shaft axis of 
        the femur and the vertical axis.
        '''
        lbc = self.landmarks['lat_bot_con']['coords']
        mbc = self.landmarks['med_bot_con']['coords']
        con_mid = self.landmarks['con_mid']['coords']
        shaft_mid = self.landmarks['shaft_centre']['coords']
        bot_vec = (lbc - mbc)
        up_vec = ([0, bot_vec[2], -bot_vec[1]])
        shaft_vec = _norm(con_mid - shaft_mid)
        angle = np.rad2deg(np.arctan2(up_vec[2] * shaft_vec[1] - up_vec[1] * shaft_vec[2],
                                      up_vec[1] * shaft_vec[1] + up_vec[2] * shaft_vec[2]))
        self.landmarks['angles']['BA_2d'] = angle

        # 3d bicondylar angle
        angle_3d = np.rad2deg(np.arccos(np.dot(up_vec, shaft_vec) / (_mag(up_vec) * _mag(shaft_vec))))
        self.landmarks['angles']['BA_3d'] = angle_3d
        return


# %%#########################################################################%%#
class TibFibLandmarks(STL):
    ## Initiate class #########################################################
    def __init__(self, file_name):
        self.file_name = file_name.split(".")[0]
        self.vertices = None
        self.vertex_lists = {}
        self.landmarks = {}

        STL.__init__(self, self.file_name)

    ## Load stl mesh vertices #################################################
    def load_vertices(self, file_path):
        self.load_stl(file_path)
        self.unique_vertices()

        return

    ## Load vertex lists for tibia/fibula landmarks #################################
    def load_landmark_lists(self, lists_path):
        # Distal portions
        # self.vertex_lists["medial_malleolus"] = cm_lists.node("medial_malleolus", lists_path) - 1
        # self.vertex_lists["lateral_malleolus"] = cm_lists.node("lateral_malleolus", lists_path) - 1
        self.vertex_lists["malleoli"] = cm_lists.node("malleoli", lists_path) - 1
        self.vertex_lists["fib_notch1"] = cm_lists.node("fib_notch1", lists_path) - 1
        self.vertex_lists["fib_notch2"] = cm_lists.node("fib_notch2", lists_path) - 1

        # Proximal portions
        self.vertex_lists["lateral_condyle"] = cm_lists.node("lateral_condyle", lists_path) - 1
        self.vertex_lists["posterior_lateral_condyle"] = cm_lists.node("posterior_lateral_condyle", lists_path) - 1
        self.vertex_lists["posterior_medial_condyle"] = cm_lists.node("posterior_medial_condyle", lists_path) - 1
        self.vertex_lists["anterior_medial_condyle"] = cm_lists.node("anterior_medial_condyle", lists_path) - 1
        self.vertex_lists["medial_condyle"] = cm_lists.node("medial_condyle", lists_path) - 1
        self.vertex_lists["tibial_tuberosity"] = cm_lists.node("tibial_tuberosity", lists_path) - 1
        self.vertex_lists["right_condyle"] = cm_lists.node("right_condyle", lists_path) - 1
        self.vertex_lists["left_condyle"] = cm_lists.node("left_condyle", lists_path) - 1
        self.vertex_lists["lat_plat"] = cm_lists.node("lat_plat", lists_path) - 1
        self.vertex_lists["med_plat"] = cm_lists.node("med_plat", lists_path) - 1
        self.vertex_lists["tib_plat"] = cm_lists.node("tib_plat", lists_path) - 1
        self.vertex_lists["tib_dist"] = cm_lists.node("tib_dist", lists_path) - 1
        return

    ## Tranform vertices into the joint coordinate system #############
    def align_JCS(self):
        tm_1 = TibFibAxes.set_origin(self.vertices, self)
        self.transform(tm_1)
        tm_2 = TibFibAxes.align_ML_axis(self.vertices, self)
        self.transform(tm_2)
        tm_4 = TibFibAxes.align_long_axis_z(self.vertices, self)
        self.transform(tm_4)
        tm_6 = TibFibAxes.align_long_axis_x(self.vertices, self)
        self.transform(tm_6)
        return tm_1.dot(tm_2).dot(tm_4).dot(tm_6)

    def malleoli(self):
        '''
        Find the medial and lateral malleoli of the tibfib and the midpoint
        '''
        tm = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
        gamma = None
        c = 0
        self.landmarks["lateral_malleolus"] = {}
        self.landmarks["medial_malleolus"] = {}

        while gamma != 0 and c<10:
            i = self.vertices[self.vertex_lists["malleoli"]].argmin(axis=0)[2]
            j = self.vertices[self.vertex_lists["malleoli"]].argmax(axis=0)[2]

            lat_mal = self.vertices[self.vertex_lists["malleoli"][i]]
            med_mal = self.vertices[self.vertex_lists["malleoli"][j]]

            self.landmarks["lateral_malleolus"]["coords"] = lat_mal
            self.landmarks["medial_malleolus"]["coords"] = med_mal

            ML_vector = med_mal - lat_mal

            # Rotation about y-axis that aligns mediolateral axis with z-axis
            vec_xz = np.array([ML_vector[0], 0, ML_vector[2]])
            vec_xz = vec_xz / norm(vec_xz)
            gamma = np.arccos(np.dot(vec_xz, [0, 0, 1]))
            if vec_xz[0] > 0:
                gamma = -gamma

                # Apply rotation about the long axis(y)
            self.vertices, tm_temp = Transform.by_radians(self.vertices, "y", gamma)
            tm = tm_temp.dot(tm)
            c = c + 1
        tm_inv = np.linalg.inv(tm)
        self.vertices = Transform.by_tm(self.vertices, tm_inv)
        self.landmarks["lateral_malleolus"]["coords"] = Transform.by_tm(self.landmarks["lateral_malleolus"]["coords"], tm_inv)
        self.landmarks["medial_malleolus"]["coords"] = Transform.by_tm(self.landmarks["medial_malleolus"]["coords"], tm_inv)

        self.landmarks["int_mal"] = {}
        #lat_mal = [self.landmarks["lateral_malleolus"]["coords"][0],self.landmarks["medial_malleolus"]["coords"][1],self.landmarks["lateral_malleolus"]["coords"][2]]
        self.landmarks["int_mal"]["coords"] = (self.landmarks["medial_malleolus"]["coords"] +
                                               self.landmarks["lateral_malleolus"]["coords"]) / 2
        squared_dist = np.sum(
            (self.landmarks["medial_malleolus"]["coords"] - self.landmarks["lateral_malleolus"]["coords"]) ** 2, axis=0)
        width = np.sqrt(squared_dist)

        self.landmarks["lengths"]["malleolar_width"] = width

        squared_dist = np.sum(
            (self.landmarks["lateral_condyle"]["coords"] - self.landmarks["lateral_malleolus"]["coords"]) ** 2, axis=0)
        length = np.sqrt(squared_dist)

        self.landmarks["lengths"]["tibial_length"] = length
        return

    def cyl_condyles(self):
        ##find medial and lateral condyle by cylinder fit
        med_con = self.vertices[self.vertex_lists["right_condyle"]]
        lat_con = self.vertices[self.vertex_lists["left_condyle"]]
        # centre gives the coords of the centre point and rad is the radius of that cylinder
        # centre, rad = gp.fitSphereAnalytic(med_con)
        W, centre, rad, err = cfit(med_con)

        self.landmarks["medial_condyle"] = {}
        self.landmarks['med_con'] = {}
        self.landmarks["med_con"]["coords"] = centre
        self.landmarks["med_con"]["radius"] = rad

        W, centre, rad, err = cfit(lat_con)

        self.landmarks['lateral_condyle'] = {}
        self.landmarks["lat_con"] = {}
        self.landmarks["lat_con"]["coords"] = centre
        self.landmarks["lat_con"]["radius"] = rad

        # align condylar axis to the z axis
        ML_vector = self.landmarks['med_con']['coords'] - self.landmarks['lat_con']['coords']

        # Rotation about y-axis that aligns mediolateral axis with z-axis
        vec_xz = np.array([ML_vector[0], 0, ML_vector[2]])
        vec_xz = vec_xz / norm(vec_xz)
        gamma = np.arccos(np.dot(vec_xz, [0, 0, 1]))
        if vec_xz[0] > 0:
            gamma = -gamma
            # Apply rotation about the long axis(y) and find the inverse transformation for the coords
        self.vertices, tm_temp = Transform.by_radians(self.vertices, "y", gamma)
        tm_inv = np.linalg.inv(tm_temp)

        ## find the lateral and medial condyle as the most lateral and medial points
        i = self.vertices[self.vertex_lists["lateral_condyle"]].argmin(axis=0)[2]
        j = self.vertices[self.vertex_lists["medial_condyle"]].argmax(axis=0)[2]

        self.landmarks["lateral_condyle"]["coords"] = Transform.by_tm(
            self.vertices[self.vertex_lists["lateral_condyle"][i]], tm_inv)
        self.landmarks["medial_condyle"]["coords"] = Transform.by_tm(
            self.vertices[self.vertex_lists["medial_condyle"][j]], tm_inv)

        #transform back
        self.vertices = Transform.by_tm(self.vertices, tm_inv)

        # Find the epicondylar midpoint
        mid_point = (self.landmarks["medial_condyle"]["coords"] + self.landmarks["lateral_condyle"]["coords"]) / 2

        self.landmarks["con_mid"] = {}
        self.landmarks["con_mid"]["coords"] = mid_point

        # condylar width
        squared_dist = np.sum(
            (self.landmarks["medial_condyle"]["coords"] - self.landmarks["lateral_condyle"]["coords"]) ** 2, axis=0)
        width = np.sqrt(squared_dist)
        self.landmarks['lengths'] = {}
        self.landmarks["lengths"]["con_width"] = width
        return
    def cyl_condyles_align(self):
        '''
        cylinder fit to condyles of the tibia to find med_con and lat_con. Align 
        the z axis to the cylinder fit condyles and find the most lateral and medial
        points for the medial and lateral condyle. Also in this orientation finds the 
        tibial tuberosity, the posterior lateral and medial condyles, 
        '''


        # align condylar axis to the z axis
        ML_vector = self.landmarks['med_con']['coords'] - self.landmarks['lat_con']['coords']

        # Rotation about y-axis that aligns mediolateral axis with z-axis
        vec_xz = np.array([ML_vector[0], 0, ML_vector[2]])
        vec_xz = vec_xz / norm(vec_xz)
        gamma = np.arccos(np.dot(vec_xz, [0, 0, 1]))
        if vec_xz[0] > 0:
            gamma = -gamma
            # Apply rotation about the long axis(y) and find the inverse transformation for the coords
        self.vertices, tm_temp = Transform.by_radians(self.vertices, "y", gamma)
        tm_inv = np.linalg.inv(tm_temp)



        ## find the tibial tuberosity
        k = self.vertices[self.vertex_lists["tibial_tuberosity"]].argmax(axis=0)[0]
        self.landmarks["tib_tub"] = {}
        self.landmarks["tib_tub"]["coords"] = Transform.by_tm(self.vertices[self.vertex_lists["tibial_tuberosity"][k]],
                                                              tm_inv)
        # self.landmarks["tib_tub"]["coords"]= Transform.by_tm(self.landmarks["tib_tub"]["coords"], tm_inv)

        ## find the posterior lateral and medial condyles
        self.landmarks["PLC"] = {}
        self.landmarks["PMC"] = {}
        l = self.vertices[self.vertex_lists["posterior_lateral_condyle"]].argmin(axis=0)[0]
        m = self.vertices[self.vertex_lists["posterior_medial_condyle"]].argmin(axis=0)[0]
        self.landmarks["PLC"]["coords"] = Transform.by_tm(
            self.vertices[self.vertex_lists["posterior_lateral_condyle"][l]], tm_inv)
        self.landmarks["PMC"]["coords"] = Transform.by_tm(
            self.vertices[self.vertex_lists["posterior_medial_condyle"][m]], tm_inv)
        # Find the postcondylar midpoint
        mid_point = (self.landmarks["PMC"]["coords"] + self.landmarks["PLC"]["coords"]) / 2

        self.landmarks["PC_mid"] = {}
        self.landmarks["PC_mid"]["coords"] = mid_point

        ## find the slope of the tibial plateau by finding the anterior condyle and the posterior condyle
        k = self.vertices[self.vertex_lists["tib_plat"]].argmax(axis=0)[0]
        l = self.vertices[self.vertex_lists['tib_plat']].argmin(axis=0)[0]

        self.landmarks["AMC"] = {}
        self.landmarks["AMC"]["coords"] = Transform.by_tm(self.vertices[self.vertex_lists["tib_plat"][k]], tm_inv)
        self.landmarks["PC"] = {}
        self.landmarks["PC"]["coords"] = Transform.by_tm(self.vertices[self.vertex_lists['tib_plat'][l]], tm_inv)

        # find tibal slope
        vec = _norm(self.vertices[self.vertex_lists["tib_plat"][k]] - self.vertices[self.vertex_lists['tib_plat'][l]])
        if vec[0] < 0:
            vec = vec * -1
        angle = np.rad2deg(np.arctan2(vec[1] * 1 - vec[0] * 0, 1 * vec[0] + 0 * vec[1]))
        self.landmarks['angles'] = {}
        self.landmarks['angles']['tib_slope'] = angle

        ## find mMPTA
        # mMTPA con is between the lateral and medial cylinder fit condyle
        con_vec = _norm(Transform.by_tm(self.landmarks['med_con']['coords'], tm_temp) - Transform.by_tm(
            self.landmarks['lat_con']['coords'], tm_temp))
        if con_vec[2] < 0:
            con_vec = con_vec * -1
        angle_con = np.rad2deg(np.arctan2(0 * con_vec[1] - -1 * con_vec[2], -1 * con_vec[1] + 0 * con_vec[2]))
        self.landmarks['angles']['mMPTA_con'] = angle_con

        # mMTPA found by finding the most medial and lateral points of the tibial plateau nodes
        mc = self.vertices[
            self.vertex_lists["tib_plat"][self.vertices[self.vertex_lists["tib_plat"]].argmin(axis=0)[2]]]
        lc = self.vertices[
            self.vertex_lists["tib_plat"][self.vertices[self.vertex_lists["tib_plat"]].argmax(axis=0)[2]]]

        self.landmarks['med_plat'] = {}
        self.landmarks['lat_plat'] = {}
        self.landmarks["med_plat"]["coords"] = Transform.by_tm(mc, tm_inv)
        self.landmarks["lat_plat"]["coords"] = Transform.by_tm(lc, tm_inv)

        vec = _norm(mc - lc)
        if vec[2] < 0:
            vec = vec * -1
        # calculate mMPTA
        angle = np.rad2deg(np.arctan2(0 * vec[1] - -1 * vec[2], -1 * vec[1] + 0 * vec[2]))
        self.landmarks["angles"]["mMPTA_2d"] = angle

        # 3d mMPTA
        angle_3d = np.rad2deg(np.arccos(np.dot(vec, [0, -1, 0]) / (_mag(vec) * _mag([0, -1, 0]))))
        self.landmarks['angles']['mMPTA_3d'] = angle_3d



        ###########transform back##############################################
        self.vertices = Transform.by_tm(self.vertices, tm_inv)
        p = gp.fitPlaneLS(self.vertices[self.vertex_lists["tib_plat"]])
        # find mMPTA plane
        vec_p = p.Y
        if p.Y[2] < 0:
            vec_p = vec_p * -1
        angle = np.rad2deg(np.arctan2(0 * vec_p[1] - -1 * vec_p[2], -1 * vec_p[1] + 0 * vec_p[2]))
        self.landmarks['angles']['mMPTA_plane'] = angle

        # find tibial slope plane
        vector = p.X
        if p.X[0] < 0:
            vector = vector * -1
        angle2 = np.rad2deg(np.arctan2(vector[1] * 1 - vector[0] * 0, 1 * vector[0] + 0 * vector[1]))
        self.landmarks["angles"]["tib_slope_plane"] = angle2






        return

    def tibial_torsion(self):
        '''
        Find the tibial torsion as the angle between the posterior tibial condyles
        and the malleoli. In 2d looking down the y axis
        '''
        # 2d tibial torsion
        pmc = self.landmarks['PMC']['coords']
        plc = self.landmarks['PLC']['coords']
        mc = self.landmarks['medial_condyle']['coords']
        lc = self.landmarks['lateral_condyle']['coords']
        mcc = self.landmarks['med_con']['coords']
        lcc = self.landmarks['lat_con']['coords']
        mm = self.landmarks['medial_malleolus']['coords']
        lm = self.landmarks['lateral_malleolus']['coords']
        prox_vec = _norm(pmc - plc)
        prox_vec_con = _norm(mc - lc)
        prox_vec_cyl = _norm(mcc - lcc)
        if prox_vec[2] < 0:
            prox_vec = prox_vec * -1
        if prox_vec_con[2] < 0:
            prox_vec_con = prox_vec_con * -1
        if prox_vec_cyl[2] < 0:
            prox_vec_cyl = prox_vec_cyl * -1
        dist_vec = _norm(mm - lm)
        if dist_vec[2] < 0:
            dist_vec = dist_vec * -1
        # calculate tibial torsion
        angle = np.rad2deg(np.arctan2(prox_vec[2] * dist_vec[0] - prox_vec[0] * dist_vec[2],
                                      prox_vec[0] * dist_vec[0] + prox_vec[2] * dist_vec[2]))
        self.landmarks["angles"]["TT_2d"] = angle

        angle_con = np.rad2deg(np.arctan2(prox_vec_con[2] * dist_vec[0] - prox_vec_con[0] * dist_vec[2],
                                          prox_vec_con[0] * dist_vec[0] + prox_vec_con[2] * dist_vec[2]))
        self.landmarks["angles"]["TT_2d_con"] = angle_con

        angle_cyl = np.rad2deg(np.arctan2(prox_vec_cyl[2] * dist_vec[0] - prox_vec_cyl[0] * dist_vec[2],
                                          prox_vec_cyl[0] * dist_vec[0] + prox_vec_cyl[2] * dist_vec[2]))
        self.landmarks["angles"]["TT_2d_cyl"] = angle_cyl

        # 3d tibial torsion
        angle_3d = np.rad2deg(np.arccos(np.dot(prox_vec, dist_vec) / (_mag(prox_vec) * _mag(dist_vec))))
        self.landmarks['angles']['TT_3d'] = angle_3d

        # tibial torsion but fit planes to the top and bottom of the tibfib
        p = gp.fitPlaneLS(self.vertices[self.vertex_lists["tib_plat"]])
        vec1 = p.Y
        if vec1[2] < 0:
            vec1 = vec1 * -1

        angle2 = np.rad2deg(
            np.arctan2(vec1[2] * dist_vec[0] - vec1[0] * dist_vec[2], vec1[0] * dist_vec[0] + vec1[2] * dist_vec[2]))
        self.landmarks["angles"]["TT_plane"] = angle2

        return

    def fib_notch_angle(self):
        '''
        find the fibular notches and the angle between the mediolateral axis
        and the line connecting the two notches
        '''

        # Define the most +ve in the z axis in fib notch 1
        i = self.vertices[self.vertex_lists["fib_notch1"]].argmin(axis=0)[2]

        # Define the most +ve in the z axis in fib notch 2
        j = self.vertices[self.vertex_lists["fib_notch2"]].argmin(axis=0)[2]

        self.landmarks["fib_notch1"] = {}
        self.landmarks["fib_notch1"]["ID"] = self.vertex_lists["fib_notch1"][i]
        self.landmarks["fib_notch1"]["coords"] = self.vertices[self.vertex_lists["fib_notch1"][i]]

        self.landmarks["fib_notch2"] = {}
        self.landmarks["fib_notch2"]["ID"] = self.vertex_lists["fib_notch2"][j]
        self.landmarks["fib_notch2"]["coords"] = self.vertices[self.vertex_lists["fib_notch2"][j]]

        # Find the angle between the mediolateral axis and the line connecting the two notches 
        vector = _norm(self.landmarks['fib_notch2']['coords'] - self.landmarks['fib_notch1']['coords'])
        angle = np.rad2deg(np.arctan2(vector[2] * 0 - vector[0] * 1, 0 * vector[0] + 1 * vector[2]))
        self.landmarks["angles"]["fib_notch_angle"] = angle
        return

# %%#
