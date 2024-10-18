"""
Code to fit meshes to a template mesh using gias3 for the femur, pelvis, and tibfib
author: Laura Carman
1/5/22
"""
import os
from subprocess import call

os.chdir("C:/Users/lcar475/Documents/Combined_shape_model/pelvis/both_hemi")
root_dir = os.getcwd()

template_mesh_path = "Left/19_3185_Left_pelvis"

segmented_path = os.path.join(root_dir, "split", "right_mirrored")
fitted_path = os.path.join(root_dir, "Right_mirrored")
#aligned_path = os.path.join(root_dir, "aligned_meshes")
if not os.path.exists(fitted_path):
    os.makedirs(fitted_path)

file_path = os.path.join(root_dir, "rbfreg_list.txt")
if os.path.exists(file_path):
    os.remove(file_path)

file_align = open("rbfreg_list.txt", "w+")
file_align.write(template_mesh_path+".ply\n")
bone_names_align = [f for f in os.listdir(segmented_path) if f.endswith('.ply')]
for case in sorted(bone_names_align):
    if "mean" not in case and "19_3185" not in case:
        file_align.write("split/right_mirrored/" + str(case) + "\n")
file_align.close()

batch_file = os.path.join(root_dir, "rbfreg_list.txt")
call(["gias-rbfreg", "-b", batch_file, "-d", fitted_path, "--outext", ".ply"])
print("- Completed alignment")
