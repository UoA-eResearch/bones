# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 13:12:53 2020

@author: lcar475
"""

# -*- coding: utf-8 -*-
"""
find prediction error between predicted meshes generated from PCA weights and the original segmented meshes. 
Uses cloudcompare to compute distances
to change:
- ccpath
- pred path
- aligned path
- data path
- cc_aligned_path
- might need to change row 60: case_name = case[:-25] depending on how your files are named
"""
#%%#########################################################################%%#
import os
#import open3d as o3d
from packages.export import Export
from subprocess import call
import numpy as np
from gias3.mesh import simplemesh, vtktools

#%%#########################################################################%%#
# Directory parameters
#os.chdir("/Users/lauracarman/Desktop/ASM")
os.chdir('/Users/lauracarman/Desktop')
root_dir = os.getcwd()
#path to cloud compare on your computer
ccpath = '/Applications/CloudCompare.app/Contents/MacOS/CloudCompare'
#ground truth meshes path
aligned_path = os.path.join(root_dir, "ASM", "aligned_meshes")
#predicted or fitted meshes path
pred_path = os.path.join(root_dir, "articulated-ssm", "pred_meshes", 'LOO_4', "split")
#where the error results will be saved
data_path = os.path.join(root_dir, "articulated-ssm", "pred_meshes", 'LOO_4', 'data', "LOO_rms")
#where the cloud compare aligned meshes will be saved
cc_aligned_path = os.path.join(root_dir, "articulated-ssm", "pred_meshes", 'LOO_4', "cc_aligned")
if not os.path.exists(cc_aligned_path):
    os.makedirs(cc_aligned_path)
#%%#########################################################################%%#
# Prepare directories
bone = "Pelvis"
if not os.path.exists(data_path):
    os.makedirs(data_path)

rms_file = os.path.join(data_path, bone + "_prediction_error_cc" + ".txt")
if os.path.exists(rms_file):
    os.remove(rms_file)
rms = open(rms_file, "w")
rms.write("Case           \tRMSE            \tMAE            \tSD            \tMin            \tMax\n")

# Setup ascii and stl folder paths
# Create list of meshes
ply_list = sorted([f for f in os.listdir(pred_path) if "Pelvis" in f and '._' not in f and f.endswith('.ply')])
for i in range(0,len(ply_list)):
    ply_file = ply_list[i]
    case = ply_file[:-25]
    pred_mesh = os.path.join(pred_path, ply_file)
    true_mesh = os.path.join(aligned_path, case+bone+'_rbfreg.ply')
    # Convert compare mesh to ascii

    #ply_path = os.path.join(aligned_path, true_mesh)
    print(true_mesh)
    try:
        call([ccpath, '-SILENT', '-M_EXPORT_FMT', 'PLY', '-o', pred_mesh, '-o', true_mesh, '-ICP'])
    except:
        print("!! There was a problem calling CloudCompare !!")
    # Calculate cloud to data distances
    ply_registered = \
    [f for f in os.listdir(pred_path) if f.startswith(ply_file[:-4] + "_REGISTERED") and f.endswith('.ply')][0]
    ply_path_registered = os.path.join(pred_path, ply_registered)
    mesh = vtktools.loadpoly(ply_path_registered)
    #mesh = o3d.io.read_triangle_mesh(ply_path_registered)
    Export.to_ascii(mesh, pred_mesh[:-4])

    # Setup paths for CloudCompare
    asc_path_registered = pred_mesh[:-4] + ".asc"
    print(asc_path_registered)
    tm_registered = [f for f in os.listdir(pred_path) if f.startswith(ply_file[:-4] + "_REGISTRATION") and f.endswith('.txt')][0]
    tm_path_registered = os.path.join(pred_path,tm_registered)

    try:
        call([ccpath, '-SILENT', '-C_EXPORT_FMT', 'ASC', '-o', asc_path_registered, '-o', true_mesh, '-C2M_Dist'])
    except:
        print("!! There was a problem calling CloudCompare !!")

    # Determine mean and max rms distances
    cm2dist_file = [f for f in os.listdir(pred_path) if f.startswith(ply_file[:-4] + "_C2M")][0]
    cm2dist_path = os.path.join(pred_path, cm2dist_file)
    print(cm2dist_path)

    vertices = np.loadtxt(cm2dist_path)

    dists = np.absolute(vertices[:, 3])
    dists_mae = np.mean(dists)
    dists_std = np.std(dists)
    dists_rms = np.sqrt(np.mean(np.square(dists)))
    dists_min = np.amin(dists)
    dists_max = np.amax(dists)

    # Remove text files
    os.remove(cm2dist_path)
    os.remove(asc_path_registered)
    #os.remove(asc_path_registered)
    os.remove(tm_path_registered)
    os.rename(ply_path_registered, os.path.join(cc_aligned_path, ply_file))

    # Write to file
    rms.write("%s\t%5.15f\t%5.15f\t%5.15f\t%5.15f\t%5.15f\n" % (case, dists_rms, dists_mae, dists_std, dists_min, dists_max))
    print(ply_file)
rms.close()