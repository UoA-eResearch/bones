'''
Code to merge pelvis, femur and tibfib meshes for creation of an articulated shape model. Perform before PCA.
author: Laura Carman
3/5/22
'''
import os
from gias3.mesh import vtktools
from gias3.mesh import simplemesh
import numpy as np

os.chdir('C:/Users/lcar475/Documents/Combined_shape_model/shape_model_hemi_left_mirrored/aligned_meshes_CS_3_hemi_mirrored')
root_dir = os.getcwd()
side = 'right'

#set directories for bones to be merged
pelvis_fold = root_dir #os.path.join(root_dir,'pelvis','right_hemi','aligned_meshes')
femur_fold = root_dir #os.path.join(root_dir,'femur', 'right', 'aligned_meshes')
tibfib_fold = root_dir #os.path.join(root_dir,'tibfib', 'right', 'aligned_meshes')

#set directory where you want the merged meshes to go
aligned_fold = os.path.join(root_dir)
if not os.path.exists(aligned_fold):
    os.makedirs(aligned_fold)

#find the bones to be merged
pelvis_l = [f for f in os.listdir(pelvis_fold) if "pelvis" in f and "Left" in f]
pelvis_r = [f for f in os.listdir(pelvis_fold) if "pelvis" in f and "Right" in f]
femur_l = [f for f in os.listdir(femur_fold) if "femur" in f and "Left" in f]
femur_r = [f for f in os.listdir(femur_fold) if "femur" in f and "Right" in f]
tibfib_l = [f for f in os.listdir(tibfib_fold) if "tibfib" in f and "Left" in f]
tibfib_r = [f for f in os.listdir(tibfib_fold) if "tibfib" in f and "Right" in f]

#loop through first bone and make sure the case is present in the rest of the bones
if side == 'left' or side == 'both':
    for pel in sorted(pelvis_l):
        case = pel[:-15] #change depending on file names
        #for fem in sorted(femur)
        if case+'Left_femur.ply' in sorted(femur_l):
            if case+'Left_tibfib.ply' in sorted(tibfib_l):
                #load first mesh
                pelvis_path_l = os.path.join(pelvis_fold, pel)
                pel_mesh_l = vtktools.loadpoly(pelvis_path_l)

                #load second mesh
                femur_path_l = os.path.join(femur_fold, case+'Left_femur.ply')
                fem_mesh_l = vtktools.loadpoly(femur_path_l)

                #load third mesh
                tibfib_path_l = os.path.join(tibfib_fold, case+'Left_tibfib.ply')
                tib_mesh_l = vtktools.loadpoly(tibfib_path_l)

                #create a new merged mesh object
                merged_mesh = simplemesh.SimpleMesh()
                #merge vertices of all meshes
                merged_mesh_vert = np.concatenate((pel_mesh_l.v,fem_mesh_l.v,tib_mesh_l.v),axis=0)
                #merge faces of all meshes, need to add number of vertices for all previous meshes so the face numbers are correct
                a = len(pel_mesh_l.v)
                b = len(fem_mesh_l.v)
                c = len(tib_mesh_l.v)

                merged_mesh_faces = np.concatenate((pel_mesh_l.f,fem_mesh_l.f+a,tib_mesh_l.f+a+b),axis=0)
                #assign vertices and faces to new mesh
                merged_mesh.v = merged_mesh_vert
                merged_mesh.f = merged_mesh_faces
                #save merged mesh
                save_path = os.path.join(aligned_fold, case+'left_side_hemi')
                vtktools.savepoly(merged_mesh,save_path+'.ply')
                print("finished case left:", case)

if side == 'right' or side == 'both':
    #loop through first bone and make sure the case is present in the rest of the bones
    for pel in sorted(pelvis_r):
        case = pel[:-25] #change depending on file names
        #for fem in sorted(femur)
        if case+'Right_femur_mirrored.ply' in sorted(femur_r):
            if case+'Right_tibfib_mirrored.ply' in sorted(tibfib_r):
                #load first mesh
                pelvis_path_r = os.path.join(pelvis_fold, pel)
                pel_mesh_r = vtktools.loadpoly(pelvis_path_r)

                #load second mesh
                femur_path_r = os.path.join(femur_fold, case+'Right_femur_mirrored.ply')
                fem_mesh_r = vtktools.loadpoly(femur_path_r)

                #load third mesh
                tibfib_path_r = os.path.join(tibfib_fold, case+'Right_tibfib_mirrored.ply')
                tib_mesh_r = vtktools.loadpoly(tibfib_path_r)

                #create a new merged mesh object
                merged_mesh = simplemesh.SimpleMesh()
                #merge vertices of all meshes
                merged_mesh_vert = np.concatenate((pel_mesh_r.v,fem_mesh_r.v,tib_mesh_r.v),axis=0)
                #merge faces of all meshes, need to add number of vertices for all previous meshes so the face numbers are correct
                a = len(pel_mesh_r.v)
                b = len(fem_mesh_r.v)
                c = len(tib_mesh_r.v)

                merged_mesh_faces = np.concatenate((pel_mesh_r.f,fem_mesh_r.f+a,tib_mesh_r.f+a+b),axis=0)
                #assign vertices and faces to new mesh
                merged_mesh.v = merged_mesh_vert
                merged_mesh.f = merged_mesh_faces
                #save merged mesh
                save_path = os.path.join(aligned_fold, case+'right_side_mirrored_hemi')
                vtktools.savepoly(merged_mesh,save_path+'.ply')
                print("finished case right:", case)
