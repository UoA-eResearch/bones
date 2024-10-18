'''
Code to split pelvis, femur, and tibfib meshes from a merged mesh created for the articulated shape model
author: laura carman
3/5/22
'''
import os
from gias3.mesh import vtktools
from gias3.mesh import simplemesh
import numpy as np

os.chdir("C:/Users/lcar475/Documents/Combined_shape_model/tibfib/both")
#os.chdir("/Volumes/Laura2TB/ASM")
root_dir = os.getcwd()

aligned_fold = os.path.join(root_dir, 'aligned_meshes') #where are the meshes to be split
split_fold = os.path.join(root_dir, "split") #where will the split meshes go
if not os.path.exists(split_fold):
    os.makedirs(split_fold)

bones = [f for f in os.listdir(aligned_fold) if f.endswith('.ply') and '._' not in f]
# mesh_path = os.path.join(aligned_fold, 'ASM_Left_side_mean.ply')

for bone in sorted(bones):
    all_bones = os.path.join(aligned_fold, bone)
    all_bone_mesh = vtktools.loadpoly(all_bones)

    tib = simplemesh.SimpleMesh()
    fib = simplemesh.SimpleMesh()

    tib.v = all_bone_mesh.v[2078:,:]
    fib.v = all_bone_mesh.v[:2078,:]

    tib.f = all_bone_mesh.f[4152:,:]-2078
    fib.f = all_bone_mesh.f[:4152, :]

    save_path_tib = os.path.join(split_fold, 'tibia', bone[:-23]+'ia_rbfreg_rigidreg.ply')
    vtktools.savepoly(tib, save_path_tib)

    save_path_fib = os.path.join(split_fold, 'fibula', bone[:-26] + 'fibula_rbfreg_rigidreg.ply')
    vtktools.savepoly(fib, save_path_fib)

    print('finished case:', bone)


