"""
functions for find landmarks and align code
author: Laura Carman 22/06/22
"""
import os
# os.add_dll_directory('OpenSim 4.5\\bin\\')

import numpy as np
import scipy
from packages.transform import Transform
from numpy.linalg import norm
from gias3.common import geoprimitives as gp
from packages.cylinder_fit import fit as cfit
from gias3.musculoskeletal import model_alignment
from gias3.musculoskeletal.bonemodels.modelcore import ACSCartesian
from gias3.common import transform3D, math
from scipy.spatial import cKDTree
from gias3.mesh import vtktools

def _norm(v):
    return v / np.linalg.norm(v)  # scipy.divide(v, scipy.sqrt((scipy.array(v) ** 2.0).sum()))


class PelvisPLY:
    def __init__(self, file_name):
        self.file_name = file_name.split(".")[0]
        self.vertices = None
        self.elements = None
        self.faces = None
        self.vertex_lists = {}
        self.landmarks = {}
        self.original_landmarks = {}
        self.acs = ACSCartesian((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))

    def load_vertices(self, file_path):
        mesh = vtktools.loadpoly(file_path)
        self.vertices = mesh.v
        self.faces = mesh.f

        return

    def load_landmark_lists(self, lists_path):
        self.vertex_lists["RHJC"] = node("RHJC", lists_path)
        self.vertex_lists["LHJC"] = node("LHJC", lists_path)
        self.vertex_lists["LASIS"] = node("LASIS", lists_path)
        self.vertex_lists["RASIS"] = node("RASIS", lists_path)
        self.vertex_lists["LPSIS"] = node("LPSIS", lists_path)
        self.vertex_lists["RPSIS"] = node("RPSIS", lists_path)
        return

    def save_ply(self, save_path):
        save_name = os.path.join(save_path, self.file_name + '.ply')
        mesh = vtktools.Writer(filename=save_name, v=self.vertices, f=self.faces)
        mesh.write()
        return

    def find_landmarks(self):
        # ASIS
        '''
        Find the point that is the most +ve in the x direction for the left and right side
        '''
        rasis = self.vertices[self.vertex_lists["RASIS"]]
        lasis = self.vertices[self.vertex_lists["LASIS"]]
        i = rasis.argmax(axis=0)[0]
        j = lasis.argmax(axis=0)[0]

        self.landmarks["RASIS"] = {}
        self.landmarks["RASIS"]["ID"] = self.vertex_lists["RASIS"][i]
        self.landmarks["RASIS"]["coords"] = rasis[i, :]

        self.landmarks["LASIS"] = {}
        self.landmarks["LASIS"]["ID"] = self.vertex_lists["LASIS"][j]
        self.landmarks["LASIS"]["coords"] = lasis[j, :]

        mid_point = (self.landmarks["RASIS"]["coords"] + self.landmarks["LASIS"]["coords"]) / 2

        self.landmarks["ASIS_mid"] = {}
        self.landmarks["ASIS_mid"]["coords"] = mid_point

        # PSIS
        '''
        Find the point that is the most -ve in the x direction for the left and right side
        '''
        rpsis = self.vertices[self.vertex_lists["RPSIS"]]
        lpsis = self.vertices[self.vertex_lists["LPSIS"]]
        i = rpsis.argmin(axis=0)[0]
        j = lpsis.argmin(axis=0)[0]

        self.landmarks["RPSIS"] = {}
        self.landmarks["RPSIS"]["ID"] = self.vertex_lists["RPSIS"][i]
        self.landmarks["RPSIS"]["coords"] = rpsis[i, :]

        self.landmarks["LPSIS"] = {}
        self.landmarks["LPSIS"]["ID"] = self.vertex_lists["LPSIS"][j]
        self.landmarks["LPSIS"]["coords"] = lpsis[j, :]

        mid_point = (self.landmarks["RPSIS"]["coords"] + self.landmarks["LPSIS"]["coords"]) / 2

        self.landmarks["PSIS_mid"] = {}
        self.landmarks["PSIS_mid"]["coords"] = mid_point

        # hip joint centre
        '''
        Find the centre of a sphere fit to the acetabulum of the pelvis
        '''
        hip_joint = self.vertices[self.vertex_lists["RHJC"]]

        # Centre gives the coords of the centre point and rad is the radius of that sphere
        centre, rad = gp.fitSphereAnalytic(hip_joint)
        self.landmarks["lengths"] = {}
        self.landmarks["RHJC"] = {}
        self.landmarks["RHJC"]["radius"] = rad
        self.landmarks["lengths"]["RHJ_diameter"] = 2 * rad
        self.landmarks["RHJC"]["coords"] = centre

        left_hip_joint = self.vertices[self.vertex_lists["LHJC"]]

        # Centre gives the coords of the centre point and rad is the radius of that sphere
        centre, rad = gp.fitSphereAnalytic(left_hip_joint)

        self.landmarks["LHJC"] = {}
        self.landmarks["LHJC"]["radius"] = rad
        self.landmarks["lengths"]["LHJ_diameter"] = 2 * rad
        self.landmarks["LHJC"]["coords"] = centre

        return

    def length_measurements(self):
        # ASIS width
        squared_dist = np.sum((self.landmarks["RASIS"]["coords"] - self.landmarks["LASIS"]["coords"]) ** 2, axis=0)
        length = np.sqrt(squared_dist)
        self.landmarks["lengths"]["ASIS_width"] = length

        # PSIS width
        squared_dist = np.sum((self.landmarks["RPSIS"]["coords"] - self.landmarks["LPSIS"]["coords"]) ** 2, axis=0)
        length = np.sqrt(squared_dist)
        self.landmarks["lengths"]["PSIS_width"] = length

        # Pelvis depth
        squared_dist = np.sum((self.landmarks["ASIS_mid"]["coords"] - self.landmarks["PSIS_mid"]["coords"]) ** 2,
                              axis=0)
        length = np.sqrt(squared_dist)
        self.landmarks["lengths"]["pelvis_depth"] = length
        return

    def init_align(self, pel_landmarks):
        # Define eigenvectors
        com = np.mean(self.vertices, axis=0)
        cov_mat = np.cov((self.vertices - com).T)
        u, s, v = np.linalg.svd(cov_mat)
        com = np.mean(self.vertices, axis=0)

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
        self.vertices = Transform.by_tm(self.vertices, tm)

        # Align first eigenvector with z-axis (instead of x) this is for the pelvis
        self.vertices, tm_temp = Transform.by_degrees(self.vertices, "y", 90)
        tm = tm_temp.dot(tm)

        # align x axis
        LASIS = np.asarray(self.vertices)[pel_landmarks[0], :]
        RASIS = np.asarray(self.vertices)[pel_landmarks[1], :]
        LPSIS = np.asarray(self.vertices)[pel_landmarks[2], :]
        RPSIS = np.asarray(self.vertices)[pel_landmarks[3], :]
        # Define posterior->anterior(x) axis
        ASIS_mid = (LASIS + RASIS) / 2
        PSIS_mid = (LPSIS + RPSIS) / 2
        AP_vector = ASIS_mid - PSIS_mid

        # Rotation about z-axis that aligns anterior posterior axis with x-axis
        vec_xy = np.array([AP_vector[0], AP_vector[1], 0])
        vec_xy = vec_xy / norm(vec_xy)
        gamma = np.arccos(np.dot(vec_xy, [1, 0, 0]))
        if vec_xy[1] > 0:
            gamma = -gamma
            # Apply rotation about the z axis
        self.vertices, tm_temp = Transform.by_radians(self.vertices, "z", gamma)
        tm = tm_temp.dot(tm)

        if np.asarray(self.vertices)[pel_landmarks[0], :][1] < 0:
            self.vertices, tm_temp = Transform.by_degrees(self.vertices, "x", 180)
            tm = tm_temp.dot(tm)

        return tm

    def init_align_ISB(self, pel_landmarks):
        '''
        initial alignment of the pelvis to ISB system with manually selected landmarks
        '''

        tm_1 = self.init_align(pel_landmarks)

        LASIS = np.asarray(self.vertices)[pel_landmarks[0], :]
        RASIS = np.asarray(self.vertices)[pel_landmarks[1], :]
        LPSIS = np.asarray(self.vertices)[pel_landmarks[2], :]
        RPSIS = np.asarray(self.vertices)[pel_landmarks[3], :]
        self.acs.update(*model_alignment.createPelvisACSISB(LASIS,
                                                            RASIS,
                                                            LPSIS,
                                                            RPSIS))
        u = np.array([self.acs.o, self.acs.o + self.acs.x, self.acs.o + self.acs.y, self.acs.o + self.acs.z])
        ut = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1]])

        tm_2 = transform3D.directAffine(u, ut)
        self.vertices = transform3D.transformAffine(self.vertices, tm_2)
        tm = tm_2.dot(tm_1)
        return tm

    def align_ISB(self):
        '''
        ISB alignment based on the found landmarks
        '''
        self.acs.update(*model_alignment.createPelvisACSISB(self.landmarks['LASIS']['coords'],
                                                            self.landmarks['RASIS']['coords'],
                                                            self.landmarks['LPSIS']['coords'],
                                                            self.landmarks['RPSIS']['coords']))

        u = np.array([self.acs.o, self.acs.o + self.acs.x, self.acs.o + self.acs.y, self.acs.o + self.acs.z])
        ut = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1]])

        tm = transform3D.directAffine(u, ut)
        self.vertices = transform3D.transformAffine(self.vertices, tm)

        return tm

    def update_acs(self):
        '''
        update anatomical coordinate system of the pelvis
        '''
        self.acs.update(*model_alignment.createPelvisACSISB(self.landmarks['LASIS']['coords'],
                                                            self.landmarks['RASIS']['coords'],
                                                            self.landmarks['LPSIS']['coords'],
                                                            self.landmarks['RPSIS']['coords']))
        return

    def back_transform_landmarks(self, tm):
        '''
        used if both landmarks in the original and new coordinate system are required
        '''
        tm_inv = np.linalg.inv(tm)
        for ld in self.landmarks:
            if ld == 'lengths' or ld == 'angles':
                continue
            self.original_landmarks[ld] = {}
            self.original_landmarks[ld]["coords"] = list(
                transform3D.transformAffine(np.asarray([self.landmarks[ld]["coords"]]), tm_inv).flatten())
            try:
                self.original_landmarks[ld]["ID"] = self.landmarks[ld]['ID']
            except KeyError:
                continue
        return

    def transform_landmarks(self, tm):
        '''
        transforms landmarks to the updated coordinate system
        '''
        for ld in self.landmarks:
            if ld == 'lengths' or ld == 'angles':
                continue
            self.landmarks[ld]["coords"] = transform3D.transformAffine(np.asarray([self.landmarks[ld]["coords"]]), tm)[
                0]
        return


class FemurPLY:
    def __init__(self, file_name, side):
        self.file_name = file_name.split(".")[0]
        self.vertices = None
        self.elements = None
        self.vertex_lists = {}
        self.landmarks = {}
        self.original_landmarks = {}
        self.acs = ACSCartesian((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        self.side = side

    def load_vertices(self, file_path):
        mesh = vtktools.loadpoly(file_path)
        self.vertices = mesh.v
        self.faces = mesh.f
        return

    def load_landmark_lists(self, lists_path):
        self.vertex_lists["fem_head"] = node("fem_head", lists_path)
        self.vertex_lists["dist_fem"] = node("dist_fem", lists_path)
        self.vertex_lists["great_troch"] = node("great_troch", lists_path)
        self.vertex_lists["fem_neck"] = node("fem_neck", lists_path)
        self.vertex_lists["lat_con"] = node("lat_con", lists_path)
        self.vertex_lists["med_con"] = node("med_con", lists_path)
        self.vertex_lists["fem_prox_shaft"] = node("fem_prox_shaft", lists_path)
        self.vertex_lists['fem_knee_surface'] = node('fem_knee_surface', lists_path)
        self.vertex_lists['lat_epi'] = node('lat_epi', lists_path)
        self.vertex_lists['med_epi'] = node('med_epi', lists_path)
        self.vertex_lists['lat_con_cyl'] = node('lat_con_cyl', lists_path)
        self.vertex_lists['med_con_cyl'] = node('med_con_cyl', lists_path)
        self.vertex_lists['con_cyl'] = node('con_cyl', lists_path)

        return

    def save_ply(self, save_path):
        save_name = os.path.join(save_path, self.file_name + '.ply')
        mesh = vtktools.Writer(filename=save_name, v=self.vertices, f=self.faces)
        mesh.write()
        return

    def find_landmarks(self):
        # head centre
        '''
        Find the head centre of the femur by a sphere fit to the femoral head.
        '''
        femur_head = self.vertices[self.vertex_lists["fem_head"]]

        # centre gives the coords of the centre point and rad is the radius of that sphere
        centre, rad = gp.fitSphereAnalytic(femur_head)
        self.landmarks["FHC"] = {}
        self.landmarks["FHC"]["coords"] = centre
        self.landmarks["FHC"]["radius"] = rad
        self.landmarks["lengths"] = {}
        self.landmarks["lengths"]["FHC_diameter"] = 2 * rad

        # neck centre
        '''
        find the centre of the femoral neck and the centre of the proximal femoral
        shaft using cylinder fit
        '''
        neck_cen = self.vertices[self.vertex_lists["fem_neck"]]

        # centre gives the coords of the centre point and rad is the radius of that cylinder
        W, neck_centre, rad, err = cfit(neck_cen)
        self.landmarks["neck_centre"] = {}
        self.landmarks["neck_centre"]["coords"] = neck_centre
        self.landmarks["neck_centre"]["radius"] = rad

        # shaft centre
        shaft_cen = self.vertices[self.vertex_lists["fem_prox_shaft"]]

        # centre gives the coords of the centre point and rad is the radius of that cylinder
        W, shaft_centre, rad, err = cfit(shaft_cen)
        self.landmarks["shaft_centre"] = {}
        self.landmarks["shaft_centre"]["coords"] = shaft_centre

        '''
        find femoral distal landmarks by
        1. fits cylinder to each femoral condyle
        2. aligns z axis with centre of cylinders (knee flexion axis)
        3. calculates landmarks
        '''

        con_cyl = self.vertices[self.vertex_lists["con_cyl"]]

        # fit cylinders
        W, con_centre, rad, err = cfit(con_cyl, [(0, 0)])

        # align z axis to cylinder
        con_vec = W
        if con_vec[2] < 0:
            con_vec *= -1

        # find medial epicondyle by closest distance to cylinder axis
        if self.side == 'left':
            con_vec_med = con_vec
            con_vec_line_med = gp.Line3D(con_vec_med, con_centre)
        else:
            con_vec_med = con_vec * -1
            con_vec_line_med = gp.Line3D(con_vec_med, con_centre)
        med_epi = self.vertices[self.vertex_lists['med_con']]
        j = 0
        dist_old = 10.0
        for p in range(len(med_epi)):
            dist_new = con_vec_line_med.calcDistanceFromPoint(med_epi[p, :])
            if dist_new < dist_old:
                dist_old = dist_new
                j = p

        # find lateral epicondyle by closest distance to cylinder axis
        if self.side == 'left':
            con_vec_lat = con_vec * -1
            con_vec_line_lat = gp.Line3D(con_vec_lat, con_centre)
        else:
            con_vec_lat = con_vec
            con_vec_line_lat = gp.Line3D(con_vec_lat, con_centre)
        lat_epi = self.vertices[self.vertex_lists['lat_con']]
        i = 0
        dist_old = 10.0
        for p in range(len(lat_epi)):
            dist_new = con_vec_line_lat.calcDistanceFromPoint(lat_epi[p, :])
            if dist_new < dist_old:
                dist_old = dist_new
                i = p

        self.landmarks['lat_epicon'] = {}
        self.landmarks["lat_epicon"]["coords"] = lat_epi[i, :]  # lat_epi[i,:]
        self.landmarks["lat_epicon"]["ID"] = self.vertex_lists['lat_con'][i]
        self.landmarks['med_epicon'] = {}
        self.landmarks["med_epicon"]["coords"] = med_epi[j, :]  # med_epi[j,:]
        self.landmarks["med_epicon"]["ID"] = self.vertex_lists['med_con'][j]

        if self.side == 'left':
            i = lat_epi.argmin(axis=0)[2]
            j = med_epi.argmax(axis=0)[2]
        else:
            i = lat_epi.argmax(axis=0)[2]
            j = med_epi.argmin(axis=0)[2]
        self.landmarks['lat_epicon_2'] = {}
        self.landmarks["lat_epicon_2"]["coords"] = lat_epi[i, :]  # lat_epi[i,:]
        self.landmarks["lat_epicon_2"]["ID"] = self.vertex_lists['lat_con'][i]
        self.landmarks['med_epicon_2'] = {}
        self.landmarks["med_epicon_2"]["coords"] = med_epi[j, :]  # med_epi[j,:]
        self.landmarks["med_epicon_2"]["ID"] = self.vertex_lists['med_con'][j]

        # re-align femur based on new epicondyles
        self.acs.update(*model_alignment.createFemurACSISB(self.landmarks['FHC']['coords'],
                                                           self.landmarks['med_epicon']['coords'],
                                                           self.landmarks['lat_epicon']['coords'],
                                                           side=self.side
                                                           ))
        u = np.array([self.acs.o, self.acs.o + self.acs.x, self.acs.o + self.acs.y, self.acs.o + self.acs.z])
        ut = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1]])

        tm = transform3D.directAffine(u, ut)
        self.vertices = transform3D.transformAffine(self.vertices, tm)

        # transform landmarks to new CS
        self.transform_landmarks(tm)

        # posterior condyles
        self.landmarks["lat_con_post"] = {}
        self.landmarks["med_con_post"] = {}
        lat_con = self.vertices[self.vertex_lists["lat_con"]]
        med_con = self.vertices[self.vertex_lists["med_con"]]
        i = lat_con.argmin(axis=0)[0]
        j = med_con.argmin(axis=0)[0]
        lat_con = lat_con[i, :]

        self.landmarks["lat_con_post"]["coords"] = lat_con
        self.landmarks["lat_con_post"]["ID"] = self.vertex_lists["lat_con"][i]
        med_con = med_con[j, :]
        self.landmarks["med_con_post"]["coords"] = med_con
        self.landmarks["med_con_post"]["ID"] = self.vertex_lists["med_con"][j]

        # epicondylar width
        squared_dist = np.sum((self.landmarks["med_epicon"]["coords"] - self.landmarks["lat_epicon"]["coords"]) ** 2,
                              axis=0)
        width = np.sqrt(squared_dist)

        self.landmarks["lengths"]["epicon_width"] = width

        # Find the epicondylar midpoint
        mid_point = (self.landmarks["med_epicon"]["coords"] + self.landmarks["lat_epicon"]["coords"]) / 2

        self.landmarks["epicon_mid"] = {}
        self.landmarks["epicon_mid"]["coords"] = mid_point

        # greater trochanter
        '''
        Define the greater trochanter of the femur by most dorsal (-x) and most 
        proximal (+y) by calculating distance from origin (epicon_mid). Find the most lateral 
        point of the greater trochanter (furthest point in both -z and x from origin)
        '''
        prox_lat_vertices = self.vertices[self.vertex_lists["great_troch"]]

        # greater trochanter proximal
        dist = 0
        for n in range(len(prox_lat_vertices)):
            y = prox_lat_vertices[n, 1] - self.landmarks["med_epicon"]["coords"][1]
            z = prox_lat_vertices[n, 2] - self.landmarks["med_epicon"]["coords"][2]
            length = np.sqrt(y ** 2 + z ** 2)
            if length > dist:
                dist = length
                i = n
        self.landmarks["great_trochant"] = {}
        self.landmarks["great_trochant"]["ID"] = self.vertex_lists['great_troch'][i]
        self.landmarks["great_trochant"]["coords"] = prox_lat_vertices[i, :]

        # greater trochanter lateral
        dist = 0
        for n in range(len(prox_lat_vertices)):
            x = prox_lat_vertices[n, 0] - self.landmarks["med_epicon"]["coords"][0]
            z = prox_lat_vertices[n, 2] - self.landmarks["med_epicon"]["coords"][2]
            length = np.sqrt(0.5 * x ** 2 + z ** 2)
            if length > dist:
                dist = length
                i = n
        self.landmarks["lat_great_trochant"] = {}
        self.landmarks["lat_great_trochant"]["ID"] = self.vertex_lists['great_troch'][i]
        self.landmarks["lat_great_trochant"]["coords"] = prox_lat_vertices[i, :]

        # bot condyles
        '''
        find the most distal points of the medial and lateral condyles of the femur
        '''
        lat_con = self.vertices[self.vertex_lists["lat_con"]]
        med_con = self.vertices[self.vertex_lists["med_con"]]
        i = lat_con.argmin(axis=0)[1]
        j = med_con.argmin(axis=0)[1]
        lat_bot_con = lat_con[i, :]
        med_bot_con = med_con[j, :]
        # lateral condyle
        self.landmarks["lat_bot_con"] = {}
        self.landmarks["lat_bot_con"]["ID"] = self.vertex_lists["lat_con"][i]
        self.landmarks["lat_bot_con"]["coords"] = lat_bot_con

        # medial condyle
        self.landmarks["med_bot_con"] = {}
        self.landmarks["med_bot_con"]["ID"] = self.vertex_lists['med_con'][j]
        self.landmarks["med_bot_con"]["coords"] = med_bot_con

        return tm

    def measurements(self):
        # femoral length
        '''
        calculate the length of the femur from the greater trochanter landmark
        to the lateral epicondyle landmark
        '''
        squared_dist = np.sum(
            (self.landmarks["lat_great_trochant"]["coords"] - self.landmarks["lat_epicon"]["coords"]) ** 2, axis=0)
        length = np.sqrt(squared_dist)

        self.landmarks["lengths"]["femoral_length"] = length

        # anteversion angle
        '''
        Calculation of the anteversion in 2d looking down the y axis. Uses the
        landmarks of the posterior condyles, head centre, and neck centre.
        '''
        lcp = self.landmarks['lat_con_post']['coords']
        mcp = self.landmarks['med_con_post']['coords']
        hjc = self.landmarks['FHC']['coords']
        nc = self.landmarks['neck_centre']['coords']
        if self.side == 'left':
            post_vec = _norm(mcp - lcp)
        else:
            post_vec = _norm(lcp - mcp)

        if post_vec[2] < 0:
            post_vec = post_vec * -1
        neck_vec = _norm(hjc - nc)
        if neck_vec[2] < 0:
            neck_vec = neck_vec * -1
        # calculate anteversion angle
        angle = np.rad2deg(np.arctan2(post_vec[2] * neck_vec[0] - post_vec[0] * neck_vec[2],
                                      post_vec[0] * neck_vec[0] + post_vec[2] * neck_vec[2]))
        self.landmarks["angles"] = {}
        self.landmarks["angles"]["AA"] = angle

        # neck shaft angle
        '''
        calculation of neck shaft angle in 2d looking down the x axis. Uses the
        landmarks of the head centre, neck centre, proximal shaft centre, and 
        condylar midpoint.
        '''
        con_mid = self.landmarks['epicon_mid']['coords']
        shaft_mid = self.landmarks['shaft_centre']['coords']
        hjc = self.landmarks['FHC']['coords']
        nc = self.landmarks['neck_centre']['coords']
        shaft_vec = _norm(con_mid - shaft_mid)
        if shaft_vec[1] > 0:
            shaft_vec = shaft_vec * -1
        neck_vec = _norm(hjc - nc)
        if neck_vec[1] < 0:
            neck_vec = neck_vec * -1
        angle = np.rad2deg(np.arctan2(shaft_vec[2] * neck_vec[1] - shaft_vec[1] * neck_vec[2],
                                      shaft_vec[1] * neck_vec[1] + shaft_vec[2] * neck_vec[2]))
        self.landmarks["angles"]["NSA"] = angle

        # mLDFA
        '''
        calculation of mechanical lateral distal femoral angle as the angle between
        the y axis and the bottom axis of the femur.
        '''
        lbc = self.landmarks['lat_bot_con']['coords']
        mbc = self.landmarks['med_bot_con']['coords']
        if self.side == 'left':
            bot_vec = _norm(lbc - mbc)
        else:
            bot_vec = _norm(mbc - lbc)
        #if bot_vec[2] > 0:
        #    bot_vec = bot_vec * -1
        y_vec = _norm(hjc - con_mid)
        angle = np.rad2deg(np.arctan2(y_vec[2] * bot_vec[1] - y_vec[1] * bot_vec[2], y_vec[1] * bot_vec[1] + y_vec[2] * bot_vec[2]))
        self.landmarks['angles']['mLDFA'] = angle

        # bicondylar angle
        '''
        calculation of the biconylar angle as the angle between the shaft axis of 
        the femur and the vertical axis.
        '''
        lbc = self.landmarks['lat_bot_con']['coords']
        mbc = self.landmarks['med_bot_con']['coords']
        con_mid = self.landmarks['epicon_mid']['coords']
        shaft_mid = self.landmarks['shaft_centre']['coords']
        if self.side == 'left':
            bot_vec = (lbc - mbc)
        else:
            bot_vec = (mbc - lbc)
        up_vec = ([0, bot_vec[2], -bot_vec[1]])
        shaft_vec = _norm(con_mid - shaft_mid)
        angle = np.rad2deg(np.arctan2(up_vec[2] * shaft_vec[1] - up_vec[1] * shaft_vec[2],
                                      up_vec[1] * shaft_vec[1] + up_vec[2] * shaft_vec[2]))
        self.landmarks['angles']['BA'] = angle

        return

    def init_align(self, fem_landmarks):
        # Define eigenvectors
        com = np.mean(self.vertices, axis=0)
        cov_mat = np.cov((self.vertices - com).T)
        u, s, v = np.linalg.svd(cov_mat)
        com = np.mean(self.vertices, axis=0)

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
        self.vertices = Transform.by_tm(self.vertices, tm)

        # Align first eigenvector with y axis instead of x
        self.vertices, tm_temp = Transform.by_degrees(self.vertices, "z", 90)
        tm = tm_temp.dot(tm)

        # align z axis
        lat_epi = np.asarray(self.vertices)[fem_landmarks[0], :]
        med_epi = np.asarray(self.vertices)[fem_landmarks[1], :]
        if self.side == 'left':
            ML_vector = med_epi - lat_epi
        else:
            ML_vector = lat_epi - med_epi
        vec_xz = np.array([ML_vector[0], 0, ML_vector[2]])
        vec_xz = vec_xz / norm(vec_xz)
        gamma = np.arccos(np.dot(vec_xz, [0, 0, 1]))
        if vec_xz[0] > 0:
            gamma = -gamma
        self.vertices, tm_temp = Transform.by_radians(self.vertices, "y", gamma)
        tm = tm_temp.dot(tm)
        return tm

    def init_align_ISB(self, fem_landmarks):
        '''
        aligns the femur to ISB coordinate system based on manually selected landmarks and femoral head nodes
        '''
        tm_1 = self.init_align(fem_landmarks)
        lat_epicon = np.asarray(self.vertices)[fem_landmarks[0], :]
        med_epicon = np.asarray(self.vertices)[fem_landmarks[1], :]
        femur_head = self.vertices[self.vertex_lists["fem_head"]]
        FHC, rad = gp.fitSphereAnalytic(femur_head)

        self.acs.update(*model_alignment.createFemurACSISB(FHC,
                                                           med_epicon,
                                                           lat_epicon,
                                                           side=self.side
                                                           ))
        u = np.array([self.acs.o, self.acs.o + self.acs.x, self.acs.o + self.acs.y, self.acs.o + self.acs.z])
        ut = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1]])

        tm_2 = transform3D.directAffine(u, ut)
        self.vertices = transform3D.transformAffine(self.vertices, tm_2)

        tm = tm_2.dot(tm_1)
        return tm

    def align_ISB(self, Pelvis):
        '''
        aligns the femur to the ISB coordinate system, consistent with the pelvis
        '''
        self.acs.update(*model_alignment.createFemurACSISB(self.landmarks['FHC']['coords'],
                                                           self.landmarks['med_epicon']['coords'],
                                                           self.landmarks['lat_epicon']['coords'],
                                                           side=self.side
                                                           ))
        # translate axes to femur system
        if self.side == 'left':
            op = Pelvis.landmarks['LHJC']['coords']
        else:
            op = Pelvis.landmarks['RHJC']['coords']
        cs_targ = np.array([op,
                            op + Pelvis.acs.x,
                            op + Pelvis.acs.y,
                            op + Pelvis.acs.z])

        of = self.landmarks['FHC']['coords']
        cs_source = np.array([of,
                              of + self.acs.x,
                              of + self.acs.y,
                              of + self.acs.z])
        tm = transform3D.directAffine(cs_source, cs_targ)
        self.vertices = transform3D.transformAffine(self.vertices, tm)
        return tm

    def update_acs(self):
        '''
        updates femoral anatomical coordinate system
        '''
        self.acs.update(*model_alignment.createFemurACSISB(self.landmarks['FHC']['coords'],
                                                           self.landmarks['med_epicon']['coords'],
                                                           self.landmarks['lat_epicon']['coords'],
                                                           side=self.side
                                                           ))
        return

    def back_transform_landmarks(self, tm):
        '''
        back transforms landmarks to the original coordinate system
        '''
        tm_inv = np.linalg.inv(tm)
        for ld in self.landmarks:
            if ld == 'lengths' or ld == 'angles':
                continue
            self.original_landmarks[ld] = {}
            self.original_landmarks[ld]["coords"] = list(
                transform3D.transformAffine(np.asarray([self.landmarks[ld]["coords"]]), tm_inv).flatten())
            try:
                self.original_landmarks[ld]["ID"] = self.landmarks[ld]['ID']
            except KeyError:
                continue
        return

    def transform_landmarks(self, tm):
        '''
        transforms landmarks to the updated coordinate system
        '''
        for ld in self.landmarks:
            if ld == 'lengths' or ld == 'angles':
                continue
            self.landmarks[ld]["coords"] = transform3D.transformAffine(np.asarray([self.landmarks[ld]["coords"]]), tm)[
                0]
        return


class TibfibPLY:
    def __init__(self, file_name, side):
        self.file_name = file_name.split(".")[0]
        self.vertices = None
        self.elements = None
        self.vertex_lists = {}
        self.landmarks = {}
        self.original_landmarks = {}
        self.acs = ACSCartesian((0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
        self.side = side

    def load_vertices(self, file_path):
        mesh = vtktools.loadpoly(file_path)
        self.vertices = mesh.v
        self.faces = mesh.f
        return

    def load_landmark_lists(self, lists_path):
        self.vertex_lists["lat_con_tib"] = node("lat_con_tib", lists_path)
        self.vertex_lists["med_con_tib"] = node("med_con_tib", lists_path)
        self.vertex_lists["pos_lat_con_tib"] = node("pos_lat_con_tib", lists_path)
        self.vertex_lists["pos_med_con_tib"] = node("pos_med_con_tib", lists_path)
        self.vertex_lists["med_mal"] = node("med_mal", lists_path)
        self.vertex_lists["lat_mal"] = node("lat_mal", lists_path)
        self.vertex_lists['tib_plateau_lat'] = node('tib_plateau_lat', lists_path)
        self.vertex_lists['tib_plateau_med'] = node('tib_plateau_med', lists_path)
        self.vertex_lists['fem_knee_surface'] = node('fem_knee_surface', lists_path)
        self.vertex_lists['tib_knee_surface'] = node('tibfib_knee_surface', lists_path)

        return

    def save_ply(self, save_path):
        save_name = os.path.join(save_path, self.file_name + '.ply')
        mesh = vtktools.Writer(filename=save_name, v=self.vertices, f=self.faces)
        mesh.write()
        return

    def find_landmarks(self):
        # malleoli
        '''
        find the most distal point of the medial and lateral tibfib
        '''
        self.landmarks["lateral_malleolus"] = {}
        self.landmarks["medial_malleolus"] = {}
        self.landmarks["int_mal"] = {}
        i = self.vertices[self.vertex_lists["lat_mal"]].argmin(axis=0)[1]
        j = self.vertices[self.vertex_lists["med_mal"]].argmin(axis=0)[1]
        self.landmarks["lateral_malleolus"]["coords"] = self.vertices[self.vertex_lists["lat_mal"][i]]
        self.landmarks["lateral_malleolus"]["ID"] = self.vertex_lists["lat_mal"][i]
        self.landmarks["medial_malleolus"]["coords"] = self.vertices[self.vertex_lists["med_mal"][j]]
        self.landmarks["medial_malleolus"]["ID"] = self.vertex_lists["med_mal"][j]

        int_mal = (self.landmarks["lateral_malleolus"]["coords"] + self.landmarks["medial_malleolus"]["coords"]) / 2
        self.landmarks["int_mal"]["coords"] = int_mal

        # condyles
        '''
        aligns z axis to condyles and calculates most medial and lateral point and posterior condyle
        '''

        self.landmarks["lateral_condyle"] = {}
        self.landmarks["medial_condyle"] = {}
        self.landmarks["PLC"] = {}
        self.landmarks["PMC"] = {}
        if self.side == 'left':
            i = self.vertices[self.vertex_lists["lat_con_tib"]].argmin(axis=0)[2]
            j = self.vertices[self.vertex_lists["med_con_tib"]].argmax(axis=0)[2]
        else:
            i = self.vertices[self.vertex_lists["lat_con_tib"]].argmax(axis=0)[2]
            j = self.vertices[self.vertex_lists["med_con_tib"]].argmin(axis=0)[2]

        lat_con = self.vertices[self.vertex_lists["lat_con_tib"]][i, :]
        med_con = self.vertices[self.vertex_lists["med_con_tib"]][j, :]

        self.landmarks["lateral_condyle"]["coords"] = lat_con
        self.landmarks['lateral_condyle']['ID'] = self.vertex_lists['lat_con_tib'][i]
        self.landmarks["medial_condyle"]["coords"] = med_con
        self.landmarks['medial_condyle']['ID'] = self.vertex_lists['med_con_tib'][j]

        # posterior condyles
        i = self.vertices[self.vertex_lists["pos_lat_con_tib"]].argmin(axis=0)[0]
        j = self.vertices[self.vertex_lists["pos_med_con_tib"]].argmin(axis=0)[0]

        PLC = self.vertices[self.vertex_lists["pos_lat_con_tib"]][i, :]
        PMC = self.vertices[self.vertex_lists["pos_med_con_tib"]][j, :]

        self.landmarks["PLC"]["coords"] = PLC
        self.landmarks['PLC']['ID'] = self.vertex_lists["pos_lat_con_tib"][i]
        self.landmarks["PMC"]["coords"] = PMC
        self.landmarks['PMC']['ID'] = self.vertex_lists["pos_med_con_tib"][j]

        # Find the condylar midpoint
        mid_point = (self.landmarks["medial_condyle"]["coords"] + self.landmarks["lateral_condyle"]["coords"]) / 2

        self.landmarks["con_mid"] = {}
        self.landmarks["con_mid"]["coords"] = mid_point

        tp_lat = np.mean(self.vertices[self.vertex_lists['tib_plateau_lat']], axis=0)
        tp_med = np.mean(self.vertices[self.vertex_lists['tib_plateau_med']], axis=0)

        self.landmarks['tib_plateau_lat'] = {}
        self.landmarks['tib_plateau_lat']['coords'] = tp_lat
        self.landmarks['tib_plateau_med'] = {}
        self.landmarks['tib_plateau_med']['coords'] = tp_med

        return

    def measurements(self):
        self.landmarks["lengths"] = {}
        squared_dist = np.sum(
            (self.landmarks["medial_malleolus"]["coords"] - self.landmarks["lateral_malleolus"]["coords"]) ** 2, axis=0)
        width = np.sqrt(squared_dist)

        self.landmarks["lengths"]["malleolar_width"] = width

        squared_dist = np.sum(
            (self.landmarks["medial_condyle"]["coords"] - self.landmarks["lateral_condyle"]["coords"]) ** 2, axis=0)
        width = np.sqrt(squared_dist)

        self.landmarks["lengths"]["condylar_width"] = width

        squared_dist = np.sum(
            (self.landmarks["lateral_condyle"]["coords"] - self.landmarks["lateral_malleolus"]["coords"]) ** 2, axis=0)
        length = np.sqrt(squared_dist)

        self.landmarks["lengths"]["tibial_length"] = length
        # tibial torsion
        '''
        Find the tibial torsion as the angle between the posterior tibial condyles
        and the malleoli. In 2d looking down the y axis
        '''
        self.landmarks["angles"] = {}
        pmc = self.landmarks['PMC']['coords']
        plc = self.landmarks['PLC']['coords']
        mm = self.landmarks['medial_malleolus']['coords']
        lm = self.landmarks['lateral_malleolus']['coords']
        prox_vec = _norm(pmc - plc)
        if prox_vec[2] < 0:
            prox_vec = prox_vec * -1
        dist_vec = _norm(mm - lm)
        if dist_vec[2] < 0:
            dist_vec = dist_vec * -1
        # calculate tibial torsion
        angle = np.rad2deg(np.arctan2(prox_vec[2] * dist_vec[0] - prox_vec[0] * dist_vec[2],
                                      prox_vec[0] * dist_vec[0] + prox_vec[2] * dist_vec[2]))
        self.landmarks["angles"]["TT"] = angle

        # mMPTA
        mc = self.landmarks["tib_plateau_med"]["coords"]
        lc = self.landmarks["tib_plateau_lat"]["coords"]
        ic = self.landmarks['con_mid']['coords']
        im = self.landmarks['int_mal']['coords']
        if self.side == 'left':
            vec = _norm(mc - lc)
        else:
            vec = _norm(lc - mc)
        if vec[2] < 0:
            vec = vec * -1

        vec_2 = _norm(im - ic)
        # calculate mMPTA
        angle = np.rad2deg(np.arctan2(vec_2[2] * vec[1] - vec_2[1] * vec[2], vec_2[1] * vec[1] + vec_2[2] * vec[2]))
        self.landmarks["angles"]["mMPTA"] = angle
        return

    def init_align(self, tib_landmarks):
        # Define eigenvectors
        com = np.mean(self.vertices, axis=0)
        cov_mat = np.cov((self.vertices - com).T)
        u, s, v = np.linalg.svd(cov_mat)
        com = np.mean(self.vertices, axis=0)

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
        self.vertices = Transform.by_tm(self.vertices, tm)

        # Align first eigenvector with y axis instead of x (main eignenvector)
        self.vertices, tm_temp = Transform.by_degrees(self.vertices, "z", 90)
        tm = tm_temp.dot(tm)

        # align z axis
        lat_mal = np.asarray(self.vertices)[tib_landmarks[0], :]
        med_mal = np.asarray(self.vertices)[tib_landmarks[1], :]

        if self.side == 'left':
            ML_vector = med_mal - lat_mal
        else:
            ML_vector = lat_mal - med_mal

        # Rotation about y-axis that aligns mediolateral axis with z-axis
        vec_xz = np.array([ML_vector[0], 0, ML_vector[2]])
        vec_xz = vec_xz / norm(vec_xz)
        gamma = np.arccos(np.dot(vec_xz, [0, 0, 1]))
        if vec_xz[0] > 0:
            gamma = -gamma

        # Apply rotation about the long axis(y)
        self.vertices, tm_temp = Transform.by_radians(self.vertices, "y", gamma)
        tm = tm_temp.dot(tm)

        return tm
    def init_align_ISB(self, tib_landmarks):
        '''
        aligns the tibfib to the ISB coordinate system using manually identified landmarks
        '''

        tm_1 = self.init_align(tib_landmarks)

        lat_mal = np.asarray(self.vertices)[tib_landmarks[0], :]
        med_mal = np.asarray(self.vertices)[tib_landmarks[1], :]
        lat_con = np.asarray(self.vertices)[tib_landmarks[2], :]
        med_con = np.asarray(self.vertices)[tib_landmarks[3], :]

        # self.acs.update(*model_alignment.createTibiaFibulaACSISB_2(med_mal,
        #                                                            lat_mal,
        #                                                            med_con,
        #                                                            lat_con,
        #                                                            side=self.side
        #                                                            ))

        self.acs.update(*model_alignment.createTibiaFibulaACSISB(med_mal,
                                                                   lat_mal,
                                                                   med_con,
                                                                   lat_con,
                                                                   side=self.side
                                                                   ))

        u = np.array([self.acs.o, self.acs.o + self.acs.x, self.acs.o +self.acs.y, self.acs.o +self.acs.z])
        ut = np.array([[0, 0, 0],
                       [1, 0, 0],
                       [0, 1, 0],
                       [0, 0, 1]])

        tm_2 = transform3D.directAffine(u, ut)
        #tm = transform3D.calcRigidAffineMatrix(u)

        #t = calcTranslationBetweenCoordSystems(ut, u)
        #tm = transform3D.calcRigidAffineMatrix(t)
        #self.vertices, tm = transform3D.transformRigid3D(self.vertices, t)

        self.vertices = transform3D.transformAffine(self.vertices, tm_2)

        tm = tm_2.dot(tm_1)

        return tm

    def align_ISB(self, Femur):
        '''
        aligns the tibfib to the ISB coordinate system consistent with the femoral coordinate system
        :param Femur:
        :return:
        '''
        # self.acs.update(*model_alignment.createTibiaFibulaACSISB_2(self.landmarks['medial_malleolus']['coords'],
        #                                                            self.landmarks['lateral_malleolus']['coords'],
        #                                                            self.landmarks['medial_condyle']['coords'],
        #                                                            self.landmarks['lateral_condyle']['coords'],
        #                                                            side=self.side
        #                                                            ))

        self.acs.update(*model_alignment.createTibiaFibulaACSISB(self.landmarks['medial_malleolus']['coords'],
                                                                   self.landmarks['lateral_malleolus']['coords'],
                                                                   self.landmarks['medial_condyle']['coords'],
                                                                   self.landmarks['lateral_condyle']['coords'],
                                                                   side=self.side
                                                                   ))       

        of = Femur.acs.o
        cs_targ = np.array([of,
                            of + Femur.acs.x,
                            of + Femur.acs.y,
                            of + Femur.acs.z,
                            ])

        ot = 0.5 * (self.landmarks['medial_condyle']['coords'] +
                    self.landmarks['lateral_condyle']['coords'])

        cs_source = np.array([ot,
                              ot + self.acs.x,
                              ot + self.acs.y,
                              ot + self.acs.z,
                              ])

        tm = transform3D.directAffine(cs_source, cs_targ)
        self.vertices = transform3D.transformAffine(self.vertices, tm)

        return tm

    def reset_knee_gap(self, Femur):
        # reset knee gap
        current_knee_gap = self.calc_knee_gap(Femur)
        knee_shift = current_knee_gap - 5
        shift_t = knee_shift * self.acs.y
        self.vertices = self.vertices + shift_t
        tm_temp = transform3D.calcAffineMatrix(trans=shift_t)
        return tm_temp

    def reset_knee_gap_2(self, Femur):

        # build ckdtree of femoral condyle points
        femur_tree = cKDTree(np.asarray(Femur.vertices[Femur.vertex_lists['fem_knee_surface']]))

        # calculate and apply varus-valgus angle (about floating-x)
        varus_angle = self.calc_varus_angle(femur_tree, Femur)

        floating_x = math.norm(np.cross(self.acs.y, Femur.acs.z))
        self.vertices, tm_temp = transform3D.transformRotateAboutAxis(self.vertices,
                                                                      -varus_angle,
                                                                      Femur.acs.o,
                                                                      Femur.acs.o + floating_x, retmat=True)
        tm_2 = self.reset_knee_gap(Femur)
        return tm_2.dot(tm_temp)

    def calc_knee_gap(self, Femur):
        """Calculate the closest distance between points on the femoral
        knee articular surface and the tibial knee articular surface.
        """
        # evaluate points
        femur_points = np.asarray(Femur.vertices[Femur.vertex_lists["fem_knee_surface"]])
        tibia_points = np.asarray(self.vertices[self.vertex_lists['tib_knee_surface']])

        # find closest neighbours
        femur_tree = cKDTree(femur_points)
        dist, femur_points_i = femur_tree.query(tibia_points)

        # calc distance in the tibial Y for each pair
        tibia_y_dist = np.dot(femur_points[femur_points_i] - tibia_points, self.acs.y)
        return tibia_y_dist.min()

    def calc_varus_angle(self, femur_points_tree, Femur):
        # evaluate tibial plateau points
        tp_lat = self.landmarks['tib_plateau_lat']['coords']
        tp_med = self.landmarks['tib_plateau_med']['coords']

        # tibial plateau vector
        if self.side == 'left':
            tp_v = (tp_med - tp_lat)
        else:
            tp_v = (tp_lat - tp_med)

        fc_med_i = Femur.landmarks['med_bot_con']['coords']
        fc_lat_i = Femur.landmarks['lat_bot_con']['coords']

        if self.side == 'left':
            fc_v = fc_med_i - fc_lat_i
        else:
            fc_v = fc_lat_i - fc_med_i
        # project vectors to Z-Y plane
        z = Femur.acs.z
        y = self.acs.y
        tp_v_zy = np.array([np.dot(tp_v, z), np.dot(tp_v, y)])
        fc_v_zy = np.array([np.dot(fc_v, z), np.dot(fc_v, y)])

        # calc angle in Z-Y plane
        tp_fc_angle = abs(math.angle(tp_v_zy, fc_v_zy))

        # if lat condyle is higher than med condyle, negative rotation is needed
        if self.side == 'left':
            if abs(fc_lat_i[1] - tp_lat[1]) > abs(fc_med_i[1] - tp_med[1]):
                tp_fc_angle *= -1.0
        else:
            if abs(fc_lat_i[1] - tp_lat[1]) < abs(fc_med_i[1] - tp_med[1]):
                tp_fc_angle *= -1.0

        return tp_fc_angle

    def update_acs(self):
        # self.acs.update(*model_alignment.createTibiaFibulaACSISB_2(self.landmarks['medial_malleolus']['coords'],
        #                                                          self.landmarks['lateral_malleolus']['coords'],
        #                                                          self.landmarks['medial_condyle']['coords'],
        #                                                          self.landmarks['lateral_condyle']['coords'],
        #                                                          side=self.side
        #                                                          ))

        self.acs.update(*model_alignment.createTibiaFibulaACSISB(self.landmarks['medial_malleolus']['coords'],
                                                                 self.landmarks['lateral_malleolus']['coords'],
                                                                 self.landmarks['medial_condyle']['coords'],
                                                                 self.landmarks['lateral_condyle']['coords'],
                                                                 side=self.side
                                                                 ))

    def back_transform_landmarks(self, tm):
        tm_inv = np.linalg.inv(tm)
        for ld in self.landmarks:
            if ld == 'lengths' or ld == 'angles':
                continue
            self.original_landmarks[ld] = {}
            self.original_landmarks[ld]["coords"] = list(
                transform3D.transformAffine(np.asarray([self.landmarks[ld]["coords"]]), tm_inv).flatten())
            try:
                self.original_landmarks[ld]["ID"] = self.landmarks[ld]['ID']
            except KeyError:
                continue

    def transform_landmarks(self, tm):
        '''
        transforms landmarks to the updated coordinate system
        '''
        for ld in self.landmarks:
            if ld == 'lengths' or ld == 'angles':
                continue
            self.landmarks[ld]["coords"] = transform3D.transformAffine(np.asarray([self.landmarks[ld]["coords"]]), tm)[
                0]
        return


def node(file_name, file_path):
    list_path = os.path.join(file_path, file_name + ".txt")
    with open(list_path, 'r') as myfile:
        list_file = myfile.read()

    list_file = list_file.split("\n")

    node_id = np.array([]).astype(int)
    for i in list_file:
        node = i.split("..")
        if len(node) == 2:
            missing = np.arange(int(node[0]), int(node[1]) + 1).astype(int)
            node_id = np.append(node_id, missing)
        else:
            node_id = np.append(node_id, int(node[0]))

    return node_id


def calcTranslationBetweenCoordSystems(u, u1):
    '''
    calculates the rigid translation matrix between two sets of coordinate systems as t = (tx,ty,tz,rx,ry,rz)
    :param u: target coordinate system
    where u[0,:] = origin, u[1,:] = x_vector, u[2,:] = y_vector, u[3,:] = z_vector
    e.g. u = [[0,0,0],[1,0,0],[0,1,0],[0,0,1]]
    :param u1: coordinate system to be transformed
    :return t: transformation matrix between coordinate systems
    '''

    t = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    t[0] = u[0, 0] - u1[0, 0]
    t[1] = u[0, 1] - u1[0, 1]
    t[2] = u[0, 2] - u1[0, 2]
    t[3] = np.degrees(np.arccos(np.dot(u[1, :], u1[1, :])))
    t[4] = np.degrees(np.arccos(np.dot(u[2, :], u1[2, :])))
    t[5] = np.degrees(np.arccos(np.dot(u[3, :], u1[3, :])))

    print(t)
    return t

def length(v):
  return np.sqrt(np.dot(v, v))