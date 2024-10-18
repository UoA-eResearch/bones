"""
Code to align meshes to a template mesh using gias2 for the femur, pelvis, and tibfib
author: Laura Carman
1/5/22
"""
import os
from subprocess import call

os.chdir("C:/Users/lcar475/Documents/Combined_shape_model/pelvis/right_hemi")
root_dir = os.getcwd()
template_mesh = "Left_pelvis_mean"

fitted_path = os.path.join(root_dir, "fitted_meshes")
aligned_path = os.path.join(root_dir, "aligned_meshes")
if not os.path.exists(aligned_path):
    os.makedirs(aligned_path)

file_path = os.path.join(root_dir, "rigidreg_list.txt")
if os.path.exists(file_path):
    os.remove(file_path)

file_align = open("rigidreg_list.txt", "w+")
file_align.write("fitted_meshes/" + template_mesh + ".ply\n")
bone_names_align = [f for f in os.listdir(fitted_path) if f.endswith('.ply')]
for case in sorted(bone_names_align):
    if not case == template_mesh + ".ply":
        file_align.write("fitted_meshes/" + str(case) + "\n")
file_align.close()

batch_file = os.path.join(root_dir, "rigidreg_list.txt")
call(["gias-rigidreg", "corr_r", "-b", batch_file, "-d", aligned_path, "--outext", ".ply"])
print("- Completed alignment")
