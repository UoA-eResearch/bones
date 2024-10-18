'''
Code to run Principal Component Analysis using gias2
author: laura carman
3/5/22
'''

import os
from subprocess import call

os.chdir("C:/Users/lcar475/Documents/Combined_shape_model/tibia")
os.chdir('E:/SSM_final/fibula/whole_dataset')
root_dir = os.getcwd()
template_mesh = "01_07_1175_Left_fibula_rbfreg_rigidreg"
create_list = True

aligned_path = os.path.join(root_dir, "aligned_meshes")

file_path = os.path.join(root_dir, "pca_list_2.txt")

if create_list == True:
    if os.path.exists(file_path):
        os.remove(file_path)

    file_pca = open(file_path, "w+")
    file_pca.write("aligned_meshes/" + template_mesh + ".ply\n")
    bone_names_pca = [f for f in os.listdir(aligned_path) if f.endswith('.ply') and '._' not in f]
    for case in sorted(bone_names_pca):
        if not case == template_mesh + ".ply":
            file_pca.write("aligned_meshes/" + str(case) + "\n")
    file_pca.close()

file1 = open(file_path, 'r')
Lines = file1.readlines()
n = len(Lines)
print("number of bones:", n)

pca_folder = os.path.join(root_dir, "shape_model_2")
pca_file = os.path.join(root_dir, "shape_model_2", "shape_model_2")
if not os.path.exists(pca_folder):
    os.makedirs(pca_folder)
call(["gias-trainpcashapemodel", file_path, "-n", "658", "-r", "0", "1", "2", "3", "4", "-o", pca_file])
print("- Completed pca")