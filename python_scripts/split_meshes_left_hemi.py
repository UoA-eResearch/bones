'''
Code to split pelvis, femur, and tibfib meshes from a merged mesh created for the articulated shape model
author: laura carman
3/5/22
'''
import os
from gias3.mesh import vtktools
from gias3.mesh import simplemesh
import numpy as np

os.chdir("C:/Users/lcar475/Documents/Combined_shape_model/shape_model_hemi_left_mirrored")
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
    fem_l = simplemesh.SimpleMesh()
    tib_l = simplemesh.SimpleMesh()

    pel.v = all_bone_mesh.v[0:4683,:]
    fem_l.v = all_bone_mesh.v[4683:9799, :]
    tib_l.v = all_bone_mesh.v[9799:, :]

    pel.f = all_bone_mesh.f[0:9366,:]
    fem_l.f = all_bone_mesh.f[9366:19594,:]-4683
    tib_l.f = all_bone_mesh.f[19594:, :]-9799

    save_path_pelvis = os.path.join(split_fold, bone[:-4]+'_Left_pelvis.ply')
    vtktools.savepoly(pel, save_path_pelvis)

    save_path_fem_l = os.path.join(split_fold, bone[:-4] + '_Left_femur.ply')
    vtktools.savepoly(fem_l, save_path_fem_l)

    save_path_tib_l = os.path.join(split_fold, bone[:-4] + '_Left_tibfib.ply')
    vtktools.savepoly(tib_l, save_path_tib_l)


    print('finished case:', bone)


