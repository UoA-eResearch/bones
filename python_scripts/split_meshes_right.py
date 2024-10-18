'''
Code to split pelvis, femur, and tibfib meshes from a merged mesh created for the articulated shape model
author: laura carman
3/5/22
'''
import os
from gias3.mesh import vtktools
from gias3.mesh import simplemesh
import numpy as np

os.chdir("C:/Users/lcar475/Documents/Combined_shape_model/shape_model_right_side")
#os.chdir("/Volumes/Laura2TB/ASM")
root_dir = os.getcwd()

aligned_fold = os.path.join(root_dir) #where are the meshes to be split
split_fold = os.path.join(root_dir, "split") #where will the split meshes go
if not os.path.exists(split_fold):
    os.makedirs(split_fold)

bones = [f for f in os.listdir(aligned_fold) if f.endswith('.ply') and '._' not in f]
# mesh_path = os.path.join(aligned_fold, 'ASM_Left_side_mean.ply')

for bone in sorted(bones):
    all_bones = os.path.join(aligned_fold, bone)
    all_bone_mesh = vtktools.loadpoly(all_bones)

    pel = simplemesh.SimpleMesh()
    fem_r = simplemesh.SimpleMesh()
    tib_r = simplemesh.SimpleMesh()

    pel.v = all_bone_mesh.v[0:9293,:]
    fem_r.v = all_bone_mesh.v[9293:14409, :]
    tib_r.v = all_bone_mesh.v[14409:21133, :]

    pel.f = all_bone_mesh.f[0:18586,:]
    fem_r.f = all_bone_mesh.f[18586:28814,:]-9293
    tib_r.f = all_bone_mesh.f[28814:42258, :]-14409

    save_path_pelvis = os.path.join(split_fold, bone[:-4]+'_Pelvis.ply')
    vtktools.savepoly(pel, save_path_pelvis)

    save_path_fem_r = os.path.join(split_fold, bone[:-4] + '_Right_femur.ply')
    vtktools.savepoly(fem_r, save_path_fem_r)

    save_path_tib_r = os.path.join(split_fold, bone[:-4] + '_Right_tibfib.ply')
    vtktools.savepoly(tib_r, save_path_tib_r)

    print('finished case:', bone)


