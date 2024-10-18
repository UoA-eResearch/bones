'''
Code to split pelvis, femur, and tibfib meshes from a merged mesh created for the articulated shape model
author: laura carman
3/5/22
'''
import os
from gias3.mesh import vtktools
from gias3.mesh import simplemesh
import numpy as np

os.chdir("C:/Users/lcar475/Documents/Combined_shape_model/pelvis/fitted_meshes")
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

    pel_l = simplemesh.SimpleMesh()
    pel_r = simplemesh.SimpleMesh()

    pel_l.v = all_bone_mesh.v[4610:,:]
    pel_r.v = all_bone_mesh.v[:4610,:]

    pel_l.f = all_bone_mesh.f[9220:,:]-4610
    pel_r.f = all_bone_mesh.f[:9220, :]

    save_path_pelvis_l = os.path.join(split_fold, bone[:-18]+'_Left_pelvis.ply')
    vtktools.savepoly(pel_l, save_path_pelvis_l)

    save_path_pel_r = os.path.join(split_fold, bone[:-18] + '_Right_pelvis.ply')
    vtktools.savepoly(pel_r, save_path_pel_r)

    print('finished case:', bone)


