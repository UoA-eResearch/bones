# -*- coding: utf-8 -*-
"""
Created by Desney Greybe, 2018.

"""
import os
import numpy as np

#%%#########################################################################%%#
class Export:
    ## Export mesh as pts and tri files #######################################
    @staticmethod
    def to_pts_tri(data, output_path):
        # Check data
        try:
            len(data.vertices)
        except:
            try:
                data.unique_vertices()
                data.connectivity()
            except:
                print("!! There is a problem with the input data !!")
                return            
            
        # Prepare file paths
        pts_path = os.path.join(output_path, data.file_name + ".pts")
        tri_path = os.path.join(output_path, data.file_name + ".tri")

        # Write out data files
        if not os.path.exists(output_path):
            os.makedirs(output_path)                      
        np.savetxt(pts_path, data.vertices, fmt = '%5.15f', delimiter=' ')
        np.savetxt(tri_path, data.elements, fmt = '%1.0f', delimiter=' ')
    
        return

    ## Save data to ascii file ################################################
    @staticmethod
    def to_ascii(data, output_path):
        # Check data
        try:
            len(data.vertices)
        except:
            try:
                data.unique_vertices()
            except:
                print("!! There is a problem with the input data !!")
                return  
            
        # Set file path
        ascii_out = output_path+".asc"

        np.savetxt(ascii_out, data.vertices, fmt = '%5.15f', delimiter='\t')
        
        return    
    
    ## Save data points to exdata file ########################################
    @staticmethod
    def to_exdata(data, output_path):
        # Check data
        try:
            len(data.vertices)
        except:
            try:
                data.unique_vertices()
            except:
                print("!! There is a problem with the input data !!")
                return   
            
        # Check for output folder
        if not os.path.exists(output_path):
            os.makedirs(output_path)            
        
        # Write to exdata file
        exdata_path = os.path.join(output_path, data.file_name + ".exdata")
        exdata_file = open(exdata_path, 'w')
        exdata_file.write(" Group name: " + data.file_name + "\n")
    
        exdata_file.write(" #Fields=1\n")
        exdata_file.write(" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
        exdata_file.write("   x.  Value index= 1, #Derivatives=0\n")
        exdata_file.write("   y.  Value index= 2, #Derivatives=0\n")
        exdata_file.write("   z.  Value index= 3, #Derivatives=0\n")
     
        for i in range(len(data.vertices)):
            exdata_file.write(" Node: %8.0f\n" % (i+1))
            exdata_file.write("\t%6.6f\n\t%6.6f\n\t%6.6f\n" % (data.vertices[i,0], data.vertices[i,1], data.vertices[i,2]))
    
        exdata_file.close();
        
        return
        
    ## Save data points to exnode file ########################################
    @staticmethod
    def to_exnode(data, output_path):
        # Check data
        try:
            len(data.vertices)
        except:
            try:
                data.unique_vertices()
            except:
                print("!! There is a problem with the input data !!")
                return   
            
        # Check for output folder
        if not os.path.exists(output_path):
            os.makedirs(output_path)            
        
        # Write to exdata file
        exdata_path = os.path.join(output_path, data.file_name + ".exnode")
        exdata_file = open(exdata_path, 'w')
        exdata_file.write(" Group name: " + data.file_name + "\n")
    
        exdata_file.write(" #Fields=1\n")
        exdata_file.write(" 1) coordinates, coordinate, rectangular cartesian, #Components=3\n")
        exdata_file.write("   x.  Value index= 1, #Derivatives=0\n")
        exdata_file.write("   y.  Value index= 2, #Derivatives=0\n")
        exdata_file.write("   z.  Value index= 3, #Derivatives=0\n")
     
        for i in range(len(data.vertices)):
            exdata_file.write(" Node: %8.0f\n" % (i+1))
            exdata_file.write("\t%6.6f\n\t%6.6f\n\t%6.6f\n" % (data.vertices[i,0], data.vertices[i,1], data.vertices[i,2]))
    
        exdata_file.close();
        
        return
       
    ## Create a view com file to load exdata files in cmgui ###################
    @staticmethod
    def viewcom_data(com_path, data_names):
        # Set data point colour
        colour = "bone"
        
        # Check whether data offset is necessary.
        data_files = len(data_names)
        if data_files > 1:
            offset = "100000"
        else:
            offset = "0"
        
        # Open com file 
        if not os.path.exists(com_path):
            os.makedirs(com_path)      
        com_path = os.path.join(com_path, "view_exdata.com")
        com_file = open(com_path, 'w');
        
        # Write subroutine
        com_file.write("$i = %s\n" % (offset))
        com_file.write("$int = %s\n\n" % (offset))
        
        com_file.write("###########################\n")
        com_file.write("sub setup_data\n{\n")
        com_file.write("${data} = \"$_[0]\";\n")
        com_file.write("${mat} = \"$_[1]\";\n\n")
        com_file.write("gfx read data ${data}.exdata;\n")
        com_file.write("$offset = $i;\n")
        com_file.write("gfx change_id group ${data} data_offset $offset;\n\n")
        com_file.write("gfx mod g_elem ${data} data_points glyph sphere general size \"1*1*1\" centre 0,0,0 font default select_on material ${mat} selected_material default_selected;\n")
        com_file.write("$i = $i + $int;\n}\n")
        com_file.write("###########################\n\n")
        
        # Write data sets
        for i in range(data_files):
            com_file.write("setup_data(\"%s\" , \"%s\")\n" % (data_names[i], colour))
        
        # Write the rest
        com_file.write("###########################\n\n")
        
        com_file.write("gfx cre win\n\n")
        
        com_file.write("gfx cre axes length 200\n")
        com_file.write("gfx draw axes\n\n")
        com_file.write("# gfx cre data_viewer\n")
        
        com_file.close();
    
        return

#%%#########################################################################%%#
class Write:
    ## Write rms distance files ###############################################
    @staticmethod
    def rms_distances(cases, data, percentiles, data_path):
        for key in data:
            # Open text file
            rms_path = os.path.join(data_path, key + ".txt")
            rms_file = open(rms_path, 'w')
            
            # Write headers
            rms_file.write("Case\tMean\tSD\tMin\tMax")
            for p in percentiles:
                rms_file.write("\tPerc " + str(p))
            rms_file.write("\n")
            
            # Write data
            for s in range(np.shape(data[key])[0]):
                rms_file.write("%s" % (cases[s]))
                for d in range(np.shape(data[key])[1]):
                    rms_file.write("\t%6.6f" % (data[key][s][d]))
                rms_file.write("\n")
                
            # Close text file
            rms_file.close()
            
        return
        
#%%#########################################################################%%#