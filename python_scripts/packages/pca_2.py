# -*- coding: utf-8 -*-
"""
Created by Desney Greybe, 2018.

"""
import os
import time
import numpy as np

from subprocess import call

#from packages.transform import Transform


# %%#########################################################################%%#
class CalculatePCA(object):
    ## Initiate class #########################################################
    def __init__(self, batch_path, mesh_path, input_path, output_path, bone):
        self.batch_path = batch_path
        self.mesh_path = mesh_path
        self.input_path = input_path
        self.output_path = output_path
        self.bone = bone
        self.batch_count = 0

    ## Create batch file for fitting ##########################################
    def create_batch_file(self):
        try:
            len(self.input_path)
        except:
            print("!! Specify path to meshes !!")
            return

        self.batch_count = self.batch_count + 1
        print("-- Creating batch file " + str(self.batch_count))

        # Start writing batch file
        self.batch_file = os.path.join(self.batch_path, self.bone + "_pca_" + str(self.batch_count) + ".txt")
        if not os.path.exists(self.batch_path):
            os.makedirs(self.batch_path)
        batch_file = open(self.batch_file, "w")

        # Create list of meshes to fit in each bone folder
        stl_list = [f for f in os.listdir(self.input_path) if f.endswith('.stl')]
        for stl_file in stl_list:
            if stl_file.startswith("mean_"):
                continue
            stl_path = os.path.join(self.input_path, stl_file)
            batch_file.write(stl_path + "\n")

        batch_file.close()

        return

    ## Perform rigid alignment of meshes ######################################
    def align_meshes(self):
        print("-- Aligning meshes")
        # Create output directory
        align_path = os.path.join(self.mesh_path, "3_aligned_pca", self.bone)
        if not os.path.exists(align_path):
            os.makedirs(align_path)
        print(self.batch_file)
        # Call the aligning function
        call(["gias-rigidreg", "corr_r", "-b", self.batch_file, "-d", align_path, "--outext", ".stl"])

        # Rename output files without the "_rigidreg" suffix
        stl_list = os.listdir(align_path)
        for stl_file in stl_list:
            if stl_file.endswith("_rigidreg.stl"):
                name = stl_file.split("_rigidreg.stl")[0]
                old_name = os.path.join(align_path, stl_file)
                new_name = os.path.join(align_path, name + ".stl")
                if os.path.exists(new_name):
                    os.remove(new_name)
                os.rename(old_name, new_name)

        # Update input path
        self.input_path = align_path

        return

    ## Calculate principle components #########################################
    def calculate_pca(self, visualise, pcs):
        print("-- Calculating principle components")
        # Create output directory
        self.pca_path = os.path.join(self.output_path)
        #if not os.path.exists(self.pca_path):
            #os.makedirs(self.pca_path)
        self.pca_file = os.path.join(self.pca_path, self.bone + "_pca")

        print(self.batch_file)
        print(self.pca_file)

        # Call the pca function
        if visualise == "yes":
            print( "yes")
            call(["gias-trainpcashapemodel", self.batch_file, "-n", str(pcs), "-r", "0", "1", "2", "-o", self.pca_file,
                  "--plot_pcs", "9", "-v"])
        elif visualise == "no":
            print( "no")
            call(["gias-trainpcashapemodel", self.batch_file, "-n", str(pcs), "-r", "0", "1", "2", "-o", self.pca_file,
                  "--plot_pcs", "9"])
        else:
            print("!!! Visualisation not specified correctly !!!")

        return

    # %%#########################################################################%%#


class ShapeModel():
    ## Initiate class #########################################################
    def __init__(self, input_path, bone):
        self.input_path = input_path
        self.bone = bone

    ## Load shape model #######################################################
    def load_shape_model(self):
        print( "-- Loading shape model")
        # Set pca file path
        self.pca_file = os.path.join(self.input_path, self.bone + ".pc.npz")

        # Load shape model
        shape = np.load(self.pca_file)

        self.weights = shape['projectedWeights']
        self.vertices = shape['mean']
        self.modes = shape['modes']
        self.variance = shape['weights']
        self.stdev = np.std(shape['projectedWeights'], axis=1)

        return

    ## Reorganise shape model vertices and modes to match mean ################
    def correct_vertex_order(self):
        # Load mean mesh vertices
        try:
            len(self.mean_vertices)
        except:
            mean_mesh = STL(self.bone + "_mean.stl")
            mean_mesh.load_stl(self.input_path)
            self.mean_vertices = mean_mesh.unique_vertices()

        # Reorganise shape model vertices and modes
        modes = np.shape(self.modes)[1]
        vertices = len(self.mean_vertices)

        pca_vertices = np.reshape(self.vertices, [vertices, 3])
        pca_modes = np.reshape(self.modes, [vertices, 3 * modes])


        new_pca_vertices = np.zeros([vertices, 3])
        new_pca_modes = np.zeros([vertices, 3 * modes])
        for c in range(vertices):
            i = np.where(np.all(np.isclose(pca_vertices, self.mean_vertices[c]), axis=1))

            new_pca_vertices[c, :] = pca_vertices[i, :]
            new_pca_modes[c, :] = pca_modes[i, :]

        # Replace shape model vertices and modes
        self.vertices = np.reshape(new_pca_vertices, [vertices * 3])
        self.modes = np.reshape(new_pca_modes, [vertices * 3, modes])

        # Flag
        self.ordered = 1

        return

        ## Load mean mesh, obtain connectivity ####################################

    def load_mean_mesh(self):
        print( "-- Loading mean mesh")
        # Load vertices and elements of mean mesh
        mean_mesh = STL(self.bone + "_mean.stl")
        mean_mesh.load_stl(self.input_path)
        try:
            len(self.mean_vertices)
        except:
            self.mean_vertices = mean_mesh.unique_vertices()
        self.mean_elements = mean_mesh.connectivity()

        return

    ## Calculate mode weights from standard deviations ########################
    def calculate_weights_from_stdev(self, stdevs):
        weights = self.stdev * stdevs

        return weights

    ## Calculate mesh vertices from mode weights ##############################
    def calculate_vertices_from_weights(self, weights):
        # Pad weights to match number of modes
        full_weights = np.zeros([np.shape(self.modes)[1]])
        full_weights[:len(weights)] = weights

        # Calculate new mesh vertices
        modes = self.modes * full_weights
        vertices = self.vertices + np.sum(modes, axis=1)
        n=len(vertices)
        x = int(n/3)
        vertices = np.reshape(vertices, (x, 3))

        return vertices

    ## Export mesh from weights ###############################################
    def create_mesh_from_weights(self, mesh_path, mesh_name, weights, reflect=None):
        try:
            len(self.mean_vertices)
            len(self.mean_elements)
        except:
            self.load_mean_mesh()
        try:
            len(self.ordered)
        except:
            self.correct_vertex_order()

        # Determine new mesh vertices
        new_mesh = STL(mesh_name)
        new_mesh.vertices = self.calculate_vertices_from_weights(weights)
        new_mesh.elements = self.mean_elements


        ## Reflect left mesh (optional)
        # if reflect.lower() == "yes":
        # if "L_" in mesh_name and "L_" not in self.bone:
        # new_mesh.vertices = Transform.by_reflection(new_mesh.vertices, "yz")[0]

        # Save stl mesh
        new_mesh.save_stl(mesh_path)

        return


# %%#########################################################################%%#
class ExportPCA:
    ## Reorganise PCA weights for export ###################################%%#    
    @staticmethod
    def weights(data):
        print( "-- Reorganising and exporting PCA weights")
        # Get list of bones and subjects
        with open(data.batch_file) as f:
            batch = f.readlines()

        cases = []
        case_ids = []
        bones = []
        for i in range(len(batch)):
            case = batch[i].split("\\")[-1].split(".")[0]
            case_id = "_".join(case.split("_")[:-2])
            bone = "_".join(case.split("_")[-2:])

            cases.append(case)
            if case_id not in case_ids:
                case_ids.append(case_id)
            if bone not in bones:
                bones.append(bone)

        # Create new data array
        n_case_ids = len(case_ids)
        n_modes = np.shape(data.weights)[0]
        n_weights = np.shape(data.weights)[1]
        n_bones = len(bones)
        weights = np.zeros([n_case_ids, n_modes * n_bones])

        for d in range(n_weights):
            case = cases[d]
            case_id = "_".join(cases[d].split("_")[:-2])
            bone = "_".join(cases[d].split("_")[-2:])
            for c in range(n_case_ids):
                if case_ids[c] == case_id:
                    row = c
                    break
            for b in range(n_bones):
                if bones[b] == bone:
                    break
            for m in range(n_modes):
                column = m * n_bones + b
                weights[row, column] = data.weights[m, d]

        # Calculate percentage variance
        total_v = np.sum(data.variance)
        data.variance_percent = []
        for v in data.variance:
            data.variance_percent.append(v / total_v * 100)

        # Write out weights
        weights_path = os.path.join(data.input_path, data.bone + "_weights.txt")
        weights_file = open(weights_path, 'w')

        weights_file.write("Mode:")
        for m in range(n_modes):
            weights_file.write("\tMode %s\t\t" % (str(m + 1)))

        weights_file.write("\nMode%:")
        total_v = np.sum(data.variance)
        for m in range(n_modes):
            weights_file.write("\t%2.10f\t\t" % (data.variance[m] / total_v * 100))

        weights_file.write("\nTotal%:")
        sum_v = 0
        for m in range(n_modes):
            sum_v = sum_v + data.variance[m]
            weights_file.write("\t%2.10f\t\t" % (sum_v / total_v * 100))

        weights_file.write("\nBone:")
        for m in range(n_modes):
            for b in bones:
                weights_file.write("\t%s" % (b))
            weights_file.write("\t")

        for s in range(n_case_ids):
            weights_file.write("\n%s\t" % (case_ids[s]))
            c = 0
            for m in range(n_modes * n_bones):
                weights_file.write("%8.15f\t" % weights[s, m])
                c = c + 1
                if c % 2 == 0:
                    weights_file.write("\t")

        weights_file.close()

        return

    ## Export meshes from -3 to +3 SDs for each PC (0.2 SD steps) ##########%%#    
    @staticmethod
    def all_meshes(data, output_path, pcs):
        start = -2.0
        step_size = 0.1
        stop = 2.0 + step_size

        for pc in range(pcs):
            print( "--- Starting PC " + str(pc + 1))
            start_time = time.time()
            mesh_path = os.path.join(output_path, data.bone, data.bone + "_pc_" + str(pc + 1))
            if not os.path.exists(mesh_path):
                os.makedirs(mesh_path)

            steps = np.arange(start, stop, step_size)
            for s in range(len(steps)):
                modes = np.shape(data.modes)[1]
                stdev = np.zeros(modes)
                stdev[pc] = steps[s]

                mesh_name = data.bone + "_pc_" + str(pc + 1) + "_s" + "{0:0>2}".format(s)

                weights = data.calculate_weights_from_stdev(np.transpose(stdev))

                data.create_mesh_from_weights(mesh_path, mesh_name, weights)

                print( "---- Completed SD " + "{0:+}".format(round(steps[s], 1)))

            print( "--- Completed PC " + str(pc + 1) + " in " + "{:.2f}".format(time.time() - start_time) + "s")

        return
    ## export meshes from -2 to +2 SD using the specified number of PCs
    def all_pcs(data, output_path, pcs):
        start = -2.0
        step_size = 0.1
        stop = 2.0 + step_size
        mesh_path = os.path.join(output_path, "all_pcs")
        if not os.path.exists(mesh_path):
            os.makedirs(mesh_path)
        modes = np.shape(data.modes)[1]         
        steps = np.arange(start, stop, step_size)
        for s in range(len(steps)):
            stdev = np.zeros(modes)
            stdev[:pcs] = steps[s]          
            mesh_name = data.bone + "_s" + "{0:0>2}".format(s)
            print( "---- Completed SD " + "{0:+}".format(round(steps[s], 1)))
            weights = data.calculate_weights_from_stdev(np.transpose(stdev))
            data.create_mesh_from_weights(mesh_path, mesh_name, weights)

        return

# %%#########################################################################%%#
