# -*- coding: utf-8 -*-
"""
These objects and functions relate to fitting meshes to a template mesh using
the gias2 library (Ju Zhang).

Created by Desney Greybe, 2018.

"""
import os
import numpy as np
import open3d as o3d

from subprocess import call

from packages.export import Export
from packages.stl_mesh import STL

#ccpath = 'F:/Programs/Laura/CloudCompare/CloudCompare.exe'
#ccpath = 'C:/Program Files/CloudCompare/CloudCompare.exe'
ccpath = '/Applications/CloudCompare.app/Contents/MacOS/CloudCompare.exe'

# %%#########################################################################%%#
class FitMeshes(object):
    ## Initiate class #########################################################

    def __init__(self, mesh_path, data_path, bone):
        self.fitted_path = mesh_path
        self.data_path = data_path
        self.bone = bone
        self.fit_count = 0

    ## Calculate the rms fitting error ########################################
    def calculate_fit_error(self, mesh_fold, compare_fold, suffix, rms_path, n):
        print("--- Calculating rms fitting error")
        # Create rms error file
        #rms_path = os.path.join(self.data_path, "fitting_error")
        if not os.path.exists(rms_path):
            os.makedirs(rms_path)

        rms_file = os.path.join(rms_path, self.bone + "_rms_fit_" + str(self.fit_count) + ".txt")
        rms = open(rms_file, "w")
        rms.write("Data           \tRMSE            \tMean            \tSD            \tMax\n")

        # Setup ascii and stl folder paths
        compare_path = os.path.join(compare_fold)
        mesh_path = os.path.join(mesh_fold)
        reference_path = os.path.join(self.fitted_path)

        # Create list of meshes
        ply_list = [f for f in os.listdir(compare_path) if f.endswith('.ply')]

        for ply_file in sorted(ply_list):
            case_path = os.path.join(compare_path, ply_file)
            # Convert compare mesh to ascii
            mesh = o3d.io.read_triangle_mesh(case_path)
            Export.to_ascii(mesh, compare_path)
            
            ref_file_name = stl_mesh.file_name[:(-n+4)]+"_rbfreg_rigidreg"

            # Setup paths for CloudCompare
            asc_path = os.path.join(compare_path, stl_mesh.file_name + ".asc")
            print(asc_path)
            stl_path = os.path.join(mesh_path, ref_file_name + ".stl")
            print(stl_path)

            # Calculate cloud to data distances
            try:
                call([ccpath, '-SILENT', '-C_EXPORT_FMT', 'ASC', '-o', stl_path, '-o', asc_path, '-C2M_Dist'])
            except:
                print("!! There was a problem calling CloudCompare !!")

            # Determine mean and max rms distances

            cm2dist_file = [f for f in os.listdir(compare_path) if f.startswith(stl_mesh.file_name + "_C2M")][0]
            cm2dist_path = os.path.join(compare_path, cm2dist_file)

            vertices = np.loadtxt(cm2dist_path)

            dists = np.absolute(vertices[:, 3])
            dists_mean = np.mean(dists)
            dists_std = np.std(dists)
            dists_rms = np.sqrt(np.mean(np.square(dists)))
            dists_max = np.amax(dists)

            # Remove text files
            os.remove(cm2dist_path)
            os.remove(asc_path)

            # Write to file
            rms.write("%s\t%5.15f\t%5.15f\t%5.15f\t%5.15f\n" % (stl_mesh.file_name, dists_rms, dists_mean, dists_std, dists_max))

        rms.close()

        return

        # %%#########################################################################%%#
