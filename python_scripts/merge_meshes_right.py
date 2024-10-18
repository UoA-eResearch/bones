'''
Code to merge pelvis, femur and tibfib meshes for creation of an articulated shape model. Perform before PCA.
author: Laura Carman
3/5/22
'''
import os
from gias3.mesh import vtktools
from gias3.mesh import simplemesh
import numpy as np

os.chdir('C:/Users/lcar475/Documents/Combined_shape_model')
root_dir = os.getcwd()

#set directories for bones to be merged
pelvis_fold = os.path.join(root_dir,'pelvis','aligned_meshes')
femur_fold_r = os.path.join(root_dir,'femur', 'right', 'aligned_meshes')
tibfib_fold_r = os.path.join(root_dir,'tibfib', 'right', 'aligned_meshes')

#pelvis_fold = os.path.join(root_dir, 'pelvis-sacrum', 'aligned_meshes')
#femur_fold = os.path.join(root_dir, 'femur', 'aligned_meshes')
#tibfib_fold = os.path.join(root_dir, 'tibfib', 'aligned_meshes')

#set directory where you want the merged meshes to go
aligned_fold = os.path.join(root_dir,'aligned_meshes_right')
if not os.path.exists(aligned_fold):
    os.makedirs(aligned_fold)

#find the bones to be merged
pelvis = [f for f in os.listdir(pelvis_fold) if "Pelvis" in f]
femur_r = [f for f in os.listdir(femur_fold_r) if "femur" in f]
tibfib_r = [f for f in os.listdir(tibfib_fold_r) if "tibfib" in f]

#loop through first bone and make sure the case is present in the rest of the bones
for pel in sorted(pelvis):
    case = pel[:-26] #change depending on file names
    #for fem in sorted(femur)
    if case+'Right_femur_rbfreg_rigidreg.ply' in sorted(femur_r):
                if case + 'Right_tibfib_rbfreg_rigidreg.ply' in sorted(tibfib_r):
                    #load first mesh
                    pelvis_path = os.path.join(pelvis_fold, pel)
                    pel_mesh = vtktools.loadpoly(pelvis_path)

                    # load second mesh
                    femur_path_r = os.path.join(femur_fold_r, case + 'Right_femur_rbfreg_rigidreg.ply')
                    fem_mesh_r = vtktools.loadpoly(femur_path_r)

                    # load third mesh
                    tibfib_path_r = os.path.join(tibfib_fold_r, case + 'Right_tibfib_rbfreg_rigidreg.ply')
                    tib_mesh_r = vtktools.loadpoly(tibfib_path_r)

                    #create a new merged mesh object
                    merged_mesh = simplemesh.SimpleMesh()
                    #merge vertices of all meshes
                    merged_mesh_vert = np.concatenate((pel_mesh.v,fem_mesh_r.v,tib_mesh_r.v),axis=0)
                    #merge faces of all meshes, need to add number of vertices for all previous meshes so the face numbers are correct
                    a = len(pel_mesh.v)
                    d = len(fem_mesh_r.v)
                    e = len(tib_mesh_r.v)

                    merged_mesh_faces = np.concatenate((pel_mesh.f,fem_mesh_r.f+a,tib_mesh_r.f+a+d),axis=0)
                    #assign vertices and faces to new mesh
                    merged_mesh.v = merged_mesh_vert
                    merged_mesh.f = merged_mesh_faces
                    #save merged mesh
                    save_path = os.path.join(aligned_fold, case+'right_side')
                    vtktools.savepoly(merged_mesh,save_path+'.ply')
                    print("finished case:", case)


print(a,d,e)