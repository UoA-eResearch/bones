# -*- coding: utf-8 -*-
"""
Created on Thu Jul  2 09:19:14 2020

@author: laura
"""
import numpy as np
import os
from gias3.mesh import vtktools

def by_tm(data, tm):
    if len(np.shape(data)) == 1:
        data = np.append(data, 1)
        data = np.transpose(data)
        data = np.matmul(tm, data)
        data = np.transpose(data)
        data = data[:3]
    else:
        ones = np.ones((np.shape(data)[0],1))
        data = np.append(data, ones, axis=1)
        data = np.transpose(data)
        data = np.matmul(tm, data)
        data = np.transpose(data)    
        data = data[:,0:3]
    return data

def by_reflection(data, plane):        
    # Define transformation matrix        
    if plane == "xy":    # reflect across z-axis
        tm = np.array([[1,0,0,0],[0,1,0,0],[0,0,-1,0],[0,0,0,1]])
    if plane == "xz":    # reflect across y-axis
        tm = np.array([[1,0,0,0],[0,-1,0,0],[0,0,1,0],[0,0,0,1]])
    if plane == "yz":    # reflect across x-axis
        tm = np.array([[-1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])

    # Apply transformation
    data = by_tm(data, tm)

    return data, tm

def reflect(self, plane):
    self.v, tm = by_reflection(self.v, plane)
    f1 = np.array([self.f[:,0]])
    self.f[:,0] = self.f[:,2]
    self.f[:,2] = f1
    return tm
 
os.chdir('C:/Users/lcar475/Documents/Combined_shape_model/shape_model_hemi_left_mirrored/split')
os.chdir('C:/Users/lcar475/Downloads/Sara_pred_meshes_asm/Predicted_meshes_asm/PLB_03/split')
root_dir=os.getcwd()
mirror_path=os.path.join(root_dir)
mirrored_path='C:/Users/lcar475/Documents/Combined_shape_model/shape_model_hemi_left_mirrored/split_mirrored'
mirrored_path = 'C:/Users/lcar475/Downloads/Sara_pred_meshes_asm/Predicted_meshes_asm/PLB_03/split'
if not os.path.exists(mirrored_path):
    os.makedirs(mirrored_path)

bone_names = [f for f in os.listdir(mirror_path) if f.endswith('.ply') and 'Left' in f]

for bone in bone_names:
        # Step through meshes
        bone_path = os.path.join(mirror_path, bone)
        bone_mesh = vtktools.loadpoly(bone_path)

        reflect(bone_mesh,"xy")
        save_path_bone = os.path.join(mirrored_path, bone[:-4]+'_mirrored.ply')
        vtktools.savepoly(bone_mesh, save_path_bone)

        print(bone)