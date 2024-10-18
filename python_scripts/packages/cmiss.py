# -*- coding: utf-8 -*-
"""
The cm class allows ipdata, exdata, ipnode, exnode, ipfiel and ipelem files to
be imported (read) and exported (written). The data in these files is stored as
follows:
ipdata/exdata:
    self.data["ids"] = list of data point numbers.
    self.data[field_name]["components"] = list of component names (for "coordinates" field, these are x, y, z).
    self.data[field_name]["values"] = array of values (for "coordinates" field, these are the x, y, z ceoordinates of data points).

    self.node["ids"] = list of node numbers.
    self.node[field_name]["components"] = list of component names (for "coordinates" field, these are x, y, z).
    self.node[field_name]["values"] = array of values (for "coordinates" field, these are the x, y, z ceoordinates of data points).
    self.node[field_name]["derivatives"][1] = dictionary of key:array pairs for each derivative. If a component does not have derivatives, array will contain NaNs.
    
    self.elem["elems"] = list of element numbers.
    self.elem["nodes"] = array of node numbers associated with each element.
    self.elem[field_name][" bases"] = list of basis functions associated with each component/coordinate.
    self.elem[field_name]["base_order"] = an ordered list of basis functions used in the mesh.
    self.elem[field_name]["components"] = list of component names (for "coordinates" field, these are x, y, z).
    self.elem[field_name]["interpolation"] = string defining the interpolation for that field (eg. "l.Lagrange*c.Hermite*l.Lagrange).
    self.elem[field_name]["scale_factors"] = the number of scale factors per component.

Note: There are assumptions with certain file types:
    - There are assumed to be no node versions in any of the files.
    - There are assumed to be no derivatives in exdata/ipdata files. If 
      derivatives exist, they are removed.
    - The first three values in an ipdata file are assumed to be coordinates. 
      Any remaining values (besides scale factors) are treated as components of
      a single additional field.
    - There is assumed to be only one componenet per field variable in an 
      ipfiel file. Also, the ipfiel functions have not been tested with derivatives.
    - Elements are assumed to have the same number of coordinates and the same
      basis functions.
    - Nodes must be defined before exporting an exelem file, as the nodal values
      are used to determine interpolation functions.
    - In exelem files, only l.Lagrange and c.Hermite interpolations are supported.
    
Created by Desney Greybe, 2018.

"""
import os
import numpy as np

from numpy.linalg import norm

from packages.transform import Transform

#%%#########################################################################%%#
class cm(object):
    '''This class creates an stl object'''
    ## Initiate class #########################################################
    def __init__(self, file_name, group_name = None):
        self.file_name = file_name.split(".")[0]
        self.group_name = group_name
        
        self.data = {}
        self.data["coordinates"] = {}
        self.node = {}
        self.node["coordinates"] = {}
        self.elem = {}
        self.elem["coordinates"] = {}
        
        self.precision = 15  # decimal places to print
    
    ## Function to sort derivative keys #######################################
    def sort_derivatives(self, keys):    
        sorted_keys = []
        for k in keys:
            sorted_keys.append("".join([s for s in k if s.isdigit()]))
        sorted_keys = sorted(sorted_keys)
    
        d = {}
        for k in sorted_keys:
            if len(k) not in d:
                d[len(k)] = []
            d[len(k)].append(k)
    
        derivatives = []    
        try:
            derivatives.append(d[1][0])    
            derivatives.append(d[1][1])
            derivatives.append(d[2][0])
            derivatives.append(d[1][2])
            derivatives.append(d[2][1])
            derivatives.append(d[2][2])
            derivatives.append(d[3][0])
        except:
            pass
        
        return derivatives    
    
    #%%#####################################################################%%#
    ## Read data files ########################################################
    def read_ipdata(self, file_path, file_name = None):
    # Note: There are assumed to be 3 coordinates in an ipdata file, with no 
    #       nodal derivatives. Any remaining values are treated as components of a 
    #       single additional field.            

        # Read ipdata file
        if file_name == None:
            file_name = self.file_name
            
        with open(os.path.join(file_path, file_name + ".ipdata"), 'r') as myfile:
            ipdata_file = myfile.read()
        ipdata_file = ipdata_file.split("\n")
        
        # Obtain group name
        if self.group_name == None or self.group_name == "":
            self.group_name = ipdata_file[0].strip()
        ipdata_file = ipdata_file[1:]
        
        # Organise variables
        self.data["coordinates"]
        for i in range(len(ipdata_file)):
            temp = ipdata_file[i].split()
            if len(temp) != 0:
                if i == 0:
                    self.data["ids"] = np.array(int(temp[0]))
                    self.data["coordinates"]["values"] = np.array(map(float, temp[1:4]))
                    if len(temp) > 7:
                        field_components = (len(temp)-1)/2 + 1
                        self.data["general"] = {}
                        self.data["general"]["values"] = np.array(map(float, temp[4:field_components]))
                else:
                    self.data["ids"] = np.vstack((self.data["ids"], int(temp[0])))
                    self.data["coordinates"]["values"] = np.vstack((self.data["coordinates"]["values"], map(float, temp[1:4])))
                    if len(temp) > 7:
                        field_components = (len(temp)-1)/2 + 1
                        self.data["general"]["values"] = np.vstack((self.data["general"]["values"], map(float, temp[4:field_components])))
        
        return

    def read_exdata(self, file_path, file_name = None):
    # Note: There are assumed to be no nodal derivatives in exdata files. If 
    #       derivativer are found, they are removed.
    
        # Read exdata file
        if file_name == None:
            file_name = self.file_name
            
        with open(os.path.join(file_path, file_name + ".exdata"), 'r') as myfile:
            exdata_file = myfile.read()
        exdata_file = exdata_file.split("\n")
        
        # Get group name
        if self.group_name == None or self.group_name == "":
            self.group_name = exdata_file[0].split(":")[1].strip()
        exdata_file = exdata_file[1:]

        # Read file header
        fields = {}
        fields["field_order"] = []
        for f in range(int(exdata_file[0].split("=")[1])):
            exdata_file = exdata_file[1:]
            
            # Read field name
            field = exdata_file[0].split(",")[0].split()[1]
            fields["field_order"].append(field)
            fields[field] = {}
            
            # Read field components
            components = int(exdata_file[0].split(",")[-1].split("=")[1].strip())
            fields[field]["components"] = []
            fields[field]["derivatives"] = {}
            for i in range(components):
                exdata_file = exdata_file[1:]
                component = exdata_file[0].split(".")[0].strip()
                fields[field]["components"].append(component)
                fields[field]["derivatives"][component] = []                        
                for t in exdata_file[0].split("#"):                    
                    if "Versions" in t.split("="):
                        if int(t.split("=")[1]) != 1:
                            print ("!!! Error reading exdata: Data have multiple versions !!!")
                            return
                    if "Derivatives" in t.split("="):
                        derivatives = int(t.split("=")[1].split(",")[0].split()[0])
                        if derivatives != 0:         
                            derivatives = t.split("=")[-1].split()[1].split(",")
                            for d in derivatives:
                                if len(d) > 0:
                                    fields[field]["derivatives"][component].append(d.split("/")[1].replace("(","").replace(")",""))
        exdata_file = exdata_file[1:]
        
        # Check for derivatives
        for field in fields["field_order"]:
            for component in fields[field]["components"]:                
                if len(fields[field]["derivatives"][component]) != 0:
                    print ("!!! Warning: exdata file contains derivatives, which have been removed !!!")
                    break
            
        # Split remaining file into list of data points
        exdata_file = "\n".join(exdata_file).split("Node:")
        if len(exdata_file[0]) == 1:
            exdata_file = exdata_file[1:]
        
        # Initiate dictionaries
        self.data["ids"] = np.zeros(len(exdata_file), dtype=int)
        for field in fields["field_order"]:
            self.data[field] = {}
            self.data[field]["values"] = np.zeros((len(exdata_file), len(fields[field]["components"])))
            self.data[field]["components"] = []
            for c in fields[field]["components"]:
                if c not in self.data[field]["components"]:
                    self.data[field]["components"].append(c)
        
        # Read data points
        for i in range(len(exdata_file)):
            temp = exdata_file[0].split()
            self.data["ids"][i] = int(temp[0])
            temp = temp[1:]
            for field in fields["field_order"]:
                for c in range(len(fields[field]["components"])):
                    self.data[field]["values"][i,c] = float(temp[0])
                    temp = temp[1:]
                    component = self.data[field]["components"][c]
                    for d in fields[field]["derivatives"][component]:
                        temp = temp[1:]
            exdata_file = exdata_file[1:]
        
        return

    ## Write data files #######################################################
    def export_ipdata(self, file_path, file_name = None, group_name = None, precision = None):
        
        if file_name == None:
            file_name = self.file_name
        if group_name == None:
            if self.group_name == None:
                group_name = self.file_name
            else:
                group_name = self.group_name
        if precision == None:
            precision = self.precision

        # Ensure data dictionary is complete
        try:
            field_order = [key for key in self.data.keys() if key != "ids"]
            if "coordinates" in field_order:
                field_order.insert(0, field_order.pop(field_order.index("coordinates")))
            for field in field_order:
                if len(np.shape(self.data[field]["values"])) == 1:
                    self.data[field]["values"] = np.reshape(self.data[field]["values"], (len(self.data[field]["values"]),1))
            if "ids" not in self.data:
                self.data["ids"] = np.arange(np.shape(self.data[field_order[0]]["values"])[0]) + 1
        except:
            print ("!!! Error exporting ipdata: Data have not been defined !!!")

        # Create file directory
        if not os.path.exists(file_path):
            os.makedirs(file_path)  
        
        # Open ipdata file
        ipdata_path = os.path.join(file_path, file_name + ".ipdata")
        ipdata_file = open(ipdata_path, 'wb');
        ipdata_file.write(" %s\n" % (group_name))
    
        # Write ipdata file
        for i in range(np.shape(self.data["ids"])[0]):
            ipdata_file.write("%10i" % (self.data["ids"][i]))
            for field in field_order:
                for component in range(np.shape(self.data[field]["values"])[1]):
                    ipdata_file.write("\t% *.*f" % (precision+5, precision, self.data[field]["values"][i,component]))
            for field in field_order:
                for component in range(np.shape(self.data[field]["values"])[1]):
                    ipdata_file.write("\t1.00")
            ipdata_file.write("\n")
    
        ipdata_file.close()
                
        return
        
    def export_exdata(self, file_path, file_name = None, group_name = None, precision = None):
        
        if file_name == None:
            file_name = self.file_name
        if group_name == None:
            if self.group_name == None:
                group_name = self.file_name
            else:
                group_name = self.group_name
        if precision == None:
            precision = self.precision

        # Ensure data dictionary is complete
        try:
            field_order = [key for key in self.data.keys() if key != "ids"]
            if "coordinates" in field_order:
                field_order.insert(0, field_order.pop(field_order.index("coordinates")))
            for field in field_order:
                if len(np.shape(self.data[field]["values"])) == 1:
                    self.data[field]["values"] = np.reshape(self.data[field]["values"], (len(self.data[field]["values"]),1))
                if "components" not in self.data[field]:
                    self.data[field]["components"] = []
                    if field == "coordinates":
                        components = ['x','y','z']
                        self.data[field]["components"] = components[:np.shape(self.data["coordinates"]["values"])[1]]
                    else:
                        for c in range(np.shape(self.data[field]["values"])[1]):
                            self.data[field]["components"].append(str(c+1))
            if "ids" not in self.data:
                self.data["ids"] = np.arange(np.shape(self.data[field_order[0]]["values"])[0]) + 1
        except:
            print ("!!! Error exporting exdata: Data have not been defined !!!")
            return
 
        # Create file directory
        if not os.path.exists(file_path):
            os.makedirs(file_path)
            
        # Open exdata file
        exdata_path = os.path.join(file_path, file_name + ".exdata")
        exdata_file = open(exdata_path, 'w')
        exdata_file.write(" Group name: %s\n" % (group_name))
    
        # Write file header
        exdata_file.write(" #Fields=%i\n" % (len(field_order)))
        f, i = 1, 1
        for field in field_order:
            if field == "coordinates":
                field_type = "coordinate"
            else:
                field_type = "field"
            components = len(self.data[field]["components"])
            exdata_file.write(" %i) %s, %s, rectangular cartesian, #Components=%i\n" % (f, field, field_type, components))
            derivatives = 0
            for c in range(components):
                exdata_file.write("   %s.  Value index=%2i, #Derivatives=%2i\n" % (self.data[field]["components"][c], i, derivatives))
                i += 1
            f += 1
        
        # Write exdata file
        for i in range(np.shape(self.data["ids"])[0]):
            exdata_file.write(" Node: %8i\n" % (self.data["ids"][i]))
            for field in field_order:
                components = len(self.data[field]["components"])
                for c in range(components):
                    exdata_file.write("\t% *.*f\n" % (precision+5, precision, self.data[field]["values"][i,c]))

        exdata_file.close()
        
        return

    #%%######################################################################%#
    ## Read node files ########################################################
    # Note: There are assumed to be no versions in the ipnode and exnode files.
    def read_ipnode(self, file_path, file_name = None):
        # Read ipnode file
        if file_name == None:
            file_name = self.file_name
            
        with open(os.path.join(file_path, file_name + ".ipnode"), 'r') as myfile:
            ipnode_file = myfile.read()
        ipnode_file = ipnode_file.split("\n")[1:]
        
        # Get group name
        if self.group_name == None or self.group_name == "":
            self.group_name = ipnode_file[0].split(":")[1].strip()
        ipnode_file = ipnode_file[2:]
        
        # Read file header
        field = {}
        field["nodes"] = int(ipnode_file[0].split(":")[1].strip())
        ipnode_file = ipnode_file[1:]
        
        field["components"] = []
        field["derivivatives"] = {}
        components = int(ipnode_file[0].split(":")[1].strip())
        ipnode_file = ipnode_file[1:]
        for c in range(components):
            if ipnode_file[0].split("?")[-1].strip().lower() == "y":
                print ("!!! Error reading ipnode: Nodes have multiple versions !!!")
                return
            ipnode_file = ipnode_file[1:]
        for c in range(components):
            field["components"].append(str(c+1))
            field["derivivatives"][str(c+1)] = int(ipnode_file[0].split(":")[1].strip())
            ipnode_file = ipnode_file[1:]
        
        # Initiate dictionaries
        self.node["ids"] = np.zeros(field["nodes"], dtype=int)
        self.node["coordinates"] = {}
        self.node["coordinates"]["values"] = np.zeros((field["nodes"], len(field["components"])))
        self.node["coordinates"]["components"] = field["components"][:]
        self.node["coordinates"]["derivatives"] = {}

        # Read first node and define derivatives
        ipnode_file = ipnode_file[1:]
        self.node["ids"][0] = int(ipnode_file[0].split(":")[1].strip())
        
        derivatives = {}
        for c in range(len(field["components"])):
            ipnode_file = ipnode_file[1:]
            component = ipnode_file[0].split("(")[1].split(")")[0]
            c_i = field["components"].index(component)
            self.node["coordinates"]["values"][0,c_i] = float(ipnode_file[0].split(":")[1].strip())
            if field["derivivatives"][component] != 0:
                derivatives[component] = []
                for d in range(field["derivivatives"][component]):
                    ipnode_file = ipnode_file[1:]
                    derivative = "".join([d_i for d_i in ipnode_file[0].split("[")[0] if d_i.isdigit()])
                    derivatives[component].append(derivative)
                    if derivative not in self.node["coordinates"]["derivatives"]:
                        self.node["coordinates"]["derivatives"][derivative] = np.full((field["nodes"], len(field["components"])), np.nan)
                    self.node["coordinates"]["derivatives"][derivative][0,c_i] = float(ipnode_file[0].split(":")[1].strip())
        ipnode_file = ipnode_file[1:]     
        
        # Read remaining nodes
        for i in np.arange(1, field["nodes"]):
            ipnode_file = ipnode_file[1:]
            self.node["ids"][i] = int(ipnode_file[0].split(":")[1].strip())
            
            for c in range(len(field["components"])):
                ipnode_file = ipnode_file[1:]
                component = ipnode_file[0].split("(")[1].split(")")[0]
                c_i = field["components"].index(component)
                self.node["coordinates"]["values"][i,c_i] = float(ipnode_file[0].split(":")[1].strip())
                if field["derivivatives"][component] != 0:                    
                    for derivative in derivatives[component]:
                        ipnode_file = ipnode_file[1:]
                        self.node["coordinates"]["derivatives"][derivative][i,c_i] = float(ipnode_file[0].split(":")[1].strip())
            ipnode_file = ipnode_file[1:]
            
        # Correct coordinate names
        components = ['x','y','z']
        self.node["coordinates"]["components"] = components[:len(field["components"])]
        
        return
    
    def read_exnode(self, file_path, file_name = None):
        # Read exnode file
        if file_name == None:
            file_name = self.file_name
            
        with open(os.path.join(file_path, file_name + ".exnode"), 'r') as myfile:
            exnode_file = myfile.read()
        exnode_file = exnode_file.split("\n")
        
        # Get group name
        if self.group_name == None or self.group_name == "":
            self.group_name = exnode_file[0].split(":")[1].strip()
        exnode_file = exnode_file[1:]

        # Read file header
        fields = {}
        fields["field_order"] = []
        for f in range(int(exnode_file[0].split("=")[1])):
            exnode_file = exnode_file[1:]
            
            # Read field name
            field = exnode_file[0].split(",")[0].split()[1]
            fields["field_order"].append(field)
            fields[field] = {}
            
            # Read field components
            components = int(exnode_file[0].split(",")[-1].split("=")[1].strip())
            fields[field]["components"] = []
            fields[field]["derivatives"] = {}
            for i in range(components):
                exnode_file = exnode_file[1:]
                component = exnode_file[0].split(".")[0].strip()
                fields[field]["components"].append(component)
                fields[field]["derivatives"][component] = []                
                for t in exnode_file[0].split("#"):                    
                    if "Versions" in t.split("="):
                        if int(t.split("=")[1]) != 1:
                            print ("!!! Error reading exnode: Nodes have multiple versions !!!")
                            return
                    if "Derivatives" in t.split("="):
                        derivatives = int(t.split("=")[1].split(",")[0].split()[0])
                        if derivatives != 0:         
                            derivatives = t.split("=")[-1].split()[1].split(",")
                            for d in derivatives:
                                if len(d) > 0:
                                    temp = d.split("/")[1].replace("(","").replace(")","")
                                    fields[field]["derivatives"][component].append("".join([t for t in temp if t.isdigit()]))
        exnode_file = exnode_file[1:]
            
        # Split remaining file into list of nodes
        exnode_file = "\n".join(exnode_file).split("Node:")
        if len(exnode_file[0]) == 1:
            exnode_file = exnode_file[1:]
        
        # Initiate dictionaries
        self.node["ids"] = np.zeros(len(exnode_file), dtype=int)
        for field in fields["field_order"]:
            self.node[field] = {}
            self.node[field]["values"] = np.zeros((len(exnode_file), len(fields[field]["components"])))
            self.node[field]["components"] = []
            self.node[field]["derivatives"] = {}
            for c in fields[field]["components"]:
                if c not in self.node[field]["components"]:
                    self.node[field]["components"].append(c)
                if len(fields[field]["derivatives"][c]) != 0:
                    for d in fields[field]["derivatives"][c]:
                        if d not in self.node[field]["derivatives"]:
                            self.node[field]["derivatives"][d] = np.full((len(exnode_file), len(fields[field]["components"])), np.nan)
        
        # Read nodes
        for i in range(len(exnode_file)):
            temp = exnode_file[0].split()
            self.node["ids"][i] = int(temp[0])
            temp = temp[1:]
            for field in fields["field_order"]:
                for c in range(len(fields[field]["components"])):
                    self.node[field]["values"][i,c] = float(temp[0])
                    temp = temp[1:]
                    component = self.node[field]["components"][c]
                    for d in self.sort_derivatives(self.node[field]["derivatives"].keys()):
                        if d in fields[field]["derivatives"][component]:
                            self.node[field]["derivatives"][d][i,c] = float(temp[0])
                            temp = temp[1:]
            exnode_file = exnode_file[1:]
     
        return

    ## Write node files #######################################################
    def export_ipnode(self, file_path, file_name = None, group_name = None, precision = None):
    # Note: Only the coordinate field is written to the ipnode file. Additional 
    #       fields must be written to an ipfiel file.
    
        if file_name == None:
            file_name = self.file_name
        if group_name == None:
            if self.group_name == None:
                group_name = self.file_name
            else:
                group_name = self.group_name
        if precision == None:
            precision = self.precision
            
        # Ensure node dictionary is complete
        try:
            field_order = [key for key in self.node.keys() if key != "ids"]
            if "coordinates" in field_order:
                field_order.insert(0, field_order.pop(field_order.index("coordinates")))
            for field in field_order:
                if len(np.shape(self.node[field]["values"])) == 1:
                    self.node[field]["values"] = np.reshape(self.node[field]["values"], (len(self.node[field]["values"]),1))
                if "components" not in self.node[field]:
                    self.node[field]["components"] = []
                    if field == "coordinates":
                        components = ['x','y','z']
                        self.node[field]["components"] = components[:np.shape(self.node["coordinates"]["values"])[1]]
                    else:
                        for c in range(np.shape(self.node[field]["values"])[1]):
                            self.node[field]["components"].append(str(c+1))
                if "derivatives" not in self.node[field]:
                    self.node[field]["derivatives"] = {}
            if "ids" not in self.node:
                self.node["ids"] = np.arange(np.shape(self.node[field_order[0]]["values"])[0]) + 1
        except:
            print ("!!! Error exporting ipnode: Nodes have not been defined !!!")
            return

        # Create file directory
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        
        # Open ipnode file
        ipnode_path = os.path.join(file_path, file_name + ".ipnode")
        ipnode_file = open(ipnode_path, 'wb');

        # Write file header
        ipnode_file.write(" CMISS Version 2.1  ipnode File Version 2\n")
        ipnode_file.write(" Heading: %s\n\n" % (group_name))
        
        ipnode_file.write(" The number of nodes is [%6i]: %6i\n" % (len(self.node["ids"]), len(self.node["ids"])))
        ipnode_file.write(" Number of coordinates [3]:%2i\n" % (len(self.node["coordinates"]["components"])))
        for c in range(len(self.node["coordinates"]["components"])):
            ipnode_file.write(" Do you want prompting for different versions of nj=%i [N]? N\n" % (c+1))
        for c in range(len(self.node["coordinates"]["components"])):
            derivatives = 0
            if len(self.node["coordinates"]["derivatives"]) != 0:
                for d in self.sort_derivatives(self.node["coordinates"]["derivatives"].keys()):
                    if not np.isnan(self.node["coordinates"]["derivatives"][d][0,c]):
                        derivatives += 1
            ipnode_file.write(" The number of derivatives for coordinate %i is [0]: %i\n" % (c+1, derivatives))
        
        # Write ipnode file
        for i in range(np.shape(self.node["ids"])[0]):
            ipnode_file.write("\n Node number [%6i]: %6i" % (i + 1, self.node["ids"][i]))
            for c in range(len(self.node["coordinates"]["components"])):
                ipnode_file.write("\n The Xj(%i) coordinate is [ 0.00000E+00]:    % *.*f" % (c+1, precision+5, precision, self.node["coordinates"]["values"][i,c]))
                if len(self.node["coordinates"]["derivatives"]) != 0:
                    for d in self.sort_derivatives(self.node["coordinates"]["derivatives"].keys()):
                        if not np.isnan(self.node["coordinates"]["derivatives"][d][0,c]):
                            if len(d) == 1:
                                ipnode_file.write("\n The derivative wrt direction %s is [ 0.00000E+00]:    % *.*f" % (d, precision+5, precision, self.node["coordinates"]["derivatives"][d][i,c]))
                            if len(d) == 2:
                                ipnode_file.write("\n The derivative wrt directions %s & %s is [ 0.00000E+00]:    % *.*f" % (d[0], d[1], precision+5, precision, self.node["coordinates"]["derivatives"][d][i,c]))
                            if len(d) == 3:
                                ipnode_file.write("\n The derivative wrt directions %s, %s & %s is [ 0.00000E+00]:    % *.*f" % (d[0], d[1], d[2], precision+5, precision, self.node["coordinates"]["derivatives"][d][i,c]))
            ipnode_file.write("\n")
        
        ipnode_file.close()
        
        return

    def export_exnode(self, file_path, file_name = None, group_name = None, precision = None):
        
        if file_name == None:
            file_name = self.file_name
        if group_name == None:
            if self.group_name == None:
                group_name = self.file_name
            else:
                group_name = self.group_name
        if precision == None:
            precision = self.precision
                 
        # Ensure node dictionary is complete
        try:
            field_order = [key for key in self.node.keys() if key != "ids"]
            if "coordinates" in field_order:
                field_order.insert(0, field_order.pop(field_order.index("coordinates")))
            for field in field_order:
                if len(np.shape(self.node[field]["values"])) == 1:
                    self.node[field]["values"] = np.reshape(self.node[field]["values"], (len(self.node[field]["values"]),1))
                if "components" not in self.node[field]:
                    self.node[field]["components"] = []
                    if field == "coordinates":
                        components = ['x','y','z']
                        self.node[field]["components"] = components[:np.shape(self.node["coordinates"]["values"])[1]]
                    else:
                        for c in range(np.shape(self.node[field]["values"])[1]):
                            self.node[field]["components"].append(str(c+1))
                if "derivatives" not in self.node[field]:
                    self.node[field]["derivatives"] = {}
            if "ids" not in self.node:
                self.node["ids"] = np.arange(np.shape(self.node[field_order[0]]["values"])[0]) + 1
        except:
            print ("!!! Error exporting exnode: Nodes have not been defined !!!")
            return
 
        # Create file directory
        if not os.path.exists(file_path):
            os.makedirs(file_path)
            
        # Open exnode file
        exnode_path = os.path.join(file_path, file_name + ".exnode")
        exnode_file = open(exnode_path, 'w')
        exnode_file.write(" Group name: %s\n" % (group_name))
    
        # Write file header
        exnode_file.write(" #Fields=%i\n" % (len(field_order)))
        f, i = 1, 1
        for field in field_order:
            if field == "coordinates":
                field_type = "coordinate"
            else:
                field_type = "field"
            components = len(self.node[field]["components"])
            exnode_file.write(" %i) %s, %s, rectangular cartesian, #Components=%i\n" % (f, field, field_type, components))
            for c in range(components):
                derivatives = 0
                derivative_string = []
                if len(self.node[field]["derivatives"]) != 0:
                    for d in self.sort_derivatives(self.node[field]["derivatives"].keys()):
                        if not np.isnan(self.node[field]["derivatives"][d][0,c]):
                            derivatives += 1
                            if len(d) == 1:
                                derivative_string.append("d/ds" + d)
                            else:
                                derivative_string.append("d%i/%s" % (len(d), "ds" + "ds".join([d_i for d_i in d])))
                if len(derivative_string) == 0:
                    derivative_string = ""
                else:
                    derivative_string = " (" + ",".join(derivative_string) +")"
                exnode_file.write("   %s.  Value index=%2i, #Derivatives=%2i%s\n" % (self.node[field]["components"][c], i, derivatives, derivative_string))
                i += (derivatives + 1)
            f += 1
        
        # Write nodes to file
        for i in range(np.shape(self.node["ids"])[0]):
            exnode_file.write(" Node: %8i\n" % (self.node["ids"][i]))
            for field in field_order:
                components = len(self.node[field]["components"])
                for c in range(components):
                    exnode_file.write("\t% *.*f" % (precision+5, precision, self.node[field]["values"][i,c]))
                    if len(self.node[field]["derivatives"]) != 0:
                        for d in self.sort_derivatives(self.node[field]["derivatives"].keys()):
                            if not np.isnan(self.node[field]["derivatives"][d][0,c]):
                                exnode_file.write("\t% *.*f" % (precision+5, precision, self.node[field]["derivatives"][d][i,c]))
                    exnode_file.write("\n")

        exnode_file.close()
        
        return

    #%%#####################################################################%%#
    ## Read field files #######################################################
    # Note: There are assumed to be no versions in the ipfiel and ipfeld files.
    # Note: These functions have not been tested on field files with derivatives.
    def read_ipfiel(self, file_path, file_name = None, field_name = None):
        # Check that nodes have been defined
        try:
            len(self.node["coordinates"]["values"])
        except:
            print ("!!! Error reading ipfiel: First define nodal values !!!")
            return
        
        # Read ipfiel file
        if file_name == None:
            file_name = self.file_name
        if field_name == None:
            field_name = "general"
        
        with open(os.path.join(file_path, file_name + ".ipfiel"), 'r') as myfile:
            ipfiel_file = myfile.read()
        ipfiel_file = ipfiel_file.split("\n")[1:]
        
        # Get group name
        if self.group_name == None or self.group_name == "":
            self.group_name = ipfiel_file[0].split(":")[1].strip()
        ipfiel_file = ipfiel_file[2:]
        
        # Read file header
        components = int(ipfiel_file[0].split(":")[1].strip())
        ipfiel_file = ipfiel_file[1:]
        
        field = {}
        field["nodes"] = int(ipfiel_file[0].split(":")[1].strip())
        field["components"] = []
        field["derivivatives"] = {}
        
        ipfiel_file = ipfiel_file[components+1:]
        for c in range(components):
            field["components"].append(str(c+1))
            field["derivivatives"][str(c+1)] = int(ipfiel_file[0].split(":")[1].strip())
            ipfiel_file = ipfiel_file[1:]
        
        # Initiate dictionaries
        self.node[field_name] = {}
        self.node[field_name]["values"] = np.zeros((field["nodes"], len(field["components"])))
        self.node[field_name]["components"] = field["components"][:]
        self.node[field_name]["derivatives"] = {}
        
        # Read first node and define derivatives
        derivatives = {}
        ipfiel_file = ipfiel_file[1:]
        for c in range(len(field["components"])):
            ipfiel_file = ipfiel_file[1:]
            component = ipfiel_file[0].split()[3]
            c_i = field["components"].index(component)
            self.node[field_name]["values"][0,c_i] = float(ipfiel_file[0].split(":")[1].strip())
            if field["derivivatives"][component] != 0:
                derivatives[component] = []
                for d in range(field["derivivatives"][component]):
                    ipfiel_file = ipfiel_file[1:]
                    derivative = "".join([d_i for d_i in ipfiel_file[0].split("[")[0] if d_i.isdigit()])
                    derivatives[component].append(derivative)
                    if derivative not in self.node[field_name]["derivatives"]:
                        self.node[field_name]["derivatives"][derivative] = np.full((field["nodes"], len(field["components"])), np.nan)
                    self.node[field_name]["derivatives"][derivative][0,c_i] = float(ipfiel_file[0].split(":")[1].strip())
        ipfiel_file = ipfiel_file[1:]
        
        # Read remaining nodes
        for i in np.arange(1, field["nodes"]):
            ipfiel_file = ipfiel_file[1:]
            for c in range(len(field["components"])):
                ipfiel_file = ipfiel_file[1:]
                component = ipfiel_file[0].split()[3]
                c_i = field["components"].index(component)
                self.node[field_name]["values"][i,c_i] = float(ipfiel_file[0].split(":")[1].strip())
                if field["derivivatives"][component] != 0:
                    for derivative in derivatives[component]:
                        ipfiel_file = ipfiel_file[1:]
                        self.node[field_name]["derivatives"][derivative][i,c_i] = float(ipfiel_file[0].split(":")[1].strip())
            ipfiel_file = ipfiel_file[1:]
        
        return

    def read_ipelfd(self, file_path, file_name = None, field_name = None):
        # Read ipfeld file
        if file_name == None:
            file_name = self.file_name
        if field_name == None:
            field_name = "general"
            
        with open(os.path.join(file_path, file_name + ".ipelfd"), 'r') as myfile:
            ipelfd_file = myfile.read()
        ipelfd_file = ipelfd_file.split("\n")[1:]
        
        # Get group name
        if self.group_name == None or self.group_name == "":
            self.group_name = ipelfd_file[0].split(":")[1].strip()
        ipelfd_file = ipelfd_file[2:]
        
        # Read file header
        elements = int(ipelfd_file[0].split(":")[1].strip())
        ipelfd_file = ipelfd_file[2:]
        
        # Initiate dictionaries
        save_ids = 0
        save_nodes = 0
        if "ids" not in self.elem:
            self.elem["ids"] = np.zeros(elements, dtype=int)
            save_ids = 1
        if field_name not in self.elem:
            self.elem[field_name] = {}
        self.elem[field_name]["bases"] = []

        # Split remaining file into list of elements
        ipelfd_file = "\n".join(ipelfd_file).split("Element")
        if len(ipelfd_file[0]) == 1:
            ipelfd_file = ipelfd_file[1:]

        # Read first elements and define bases
        temp = ipelfd_file[0].split("\n")
        if save_ids == 1:
            self.elem["ids"][0] = int(temp[0].split(":")[1].strip())
        temp = temp[1:]
        
        bases = []
        for t in temp:
            if len(t.strip()) != 0:
                if t.split()[0] == "The":
                    self.elem[field_name]["bases"].append(int(t.split(":")[1].strip()))
                elif t.split()[0] == "Enter":
                        base = [b for b in t.split(":")[0] if b.isdigit()][1]
                        bases.append(int(base))
                        if "nodes" not in self.elem:
                            nodes = map(int, t.split(":")[1].split())
                            self.elem["nodes"] = np.zeros((elements, len(nodes)), dtype=int)
                            self.elem["nodes"][0,:] = nodes
                            save_nodes = 1
                
        self.elem[field_name]["base_order"] = bases[:]
        components = len(self.elem[field_name]["bases"])
        ipelfd_file = ipelfd_file[1:]

        # Read remaining elements
        if save_ids == 1 or save_nodes == 1:            
            for i in range(len(ipelfd_file)):
                temp = ipelfd_file[0].split("\n")
                if save_ids == 1:
                    self.elem["ids"][i+1] = int(temp[0].split(":")[1].strip())
                if save_nodes == 1:
                    self.elem["nodes"][i+1,:] = map(int, temp[components+2].split(":")[1].split())
                ipelfd_file = ipelfd_file[1:]
    
        return
        
    ## Write field files ######################################################
    def export_ipfiel(self, file_path, file_name = None, group_name = None, field_name = None, precision = None):
    # Note: The coordinate field must be written to an ipnode file    
        if file_name == None:
            file_name = self.file_name
        if group_name == None:
            if self.group_name == None:
                group_name = self.file_name
            else:
                group_name = self.group_name
        if field_name == None:
            if len([key for key in self.node.keys() if key not in ["ids", "coordinates"]]) != 0:
                if "general" in self.elem:
                    field_name = "general"
                    print ("!!! Warning: No field name specified. Exporting general field to ipfiel !!!")
                else:
                    print("!!! Error exporting ipfiel: Define field to export !!!")
                    return
            else:
                print ("!!! Error exporting ipfiel: No fields to export !!!")
                return
        if precision == None:
            precision = self.precision
        
        # Ensure field dictionary is complete
        try:
            components = np.shape(self.node[field_name]["values"])[1]
            if "derivatives" not in self.node[field_name]:
                self.node[field_name]["derivatives"] = {}
            if "ids" not in self.node:
                self.node["ids"] = np.arange(np.shape(self.node[field_name]["values"])[0]) + 1
        except:
            print ("!!! Error exporting ipfiel: Nodes have not been defined !!!")
            return

        # Create file directory
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        
        # Open ipfiel file
        ipfiel_path = os.path.join(file_path, file_name + ".ipfiel")
        ipfiel_file = open(ipfiel_path, 'wb');

        # Write file header
        ipfiel_file.write(" CMISS Version 2.1  ipfiel File Version 3\n")
        ipfiel_file.write(" Heading: %s\n\n" % (group_name))

        ipfiel_file.write(" The number of field variables is [3]:%2i\n" % components)
        ipfiel_file.write(" The number of nodes is [%6i]: %6i\n" % (len(self.node["ids"]), len(self.node["ids"])))
        for f in range(components):
            ipfiel_file.write(" Do you want prompting for different versions of field variable %i [N]? N\n" % (f+1))
        for f in range(components):
            derivatives = 0
            if len(self.node[field_name]["derivatives"]) != 0:
                for d in self.sort_derivatives(self.node[field_name]["derivatives"].keys()):
                    derivatives += 1
            ipfiel_file.write(" The number of derivatives for field variable %i is [0]: %i\n" % (f+1, derivatives))
        
        # Write ipfiel file
        for i in range(np.shape(self.node["ids"])[0]):
            ipfiel_file.write("\n Node number [%6i]: %6i" % (i + 1, self.node["ids"][i]))
            for f in range(components):
                ipfiel_file.write("\n The field variable %i value is [ 0.00000E+00]:    % *.*f" % (f+1, precision+5, precision, self.node[field_name]["values"][i,f]))
                if len(self.node[field_name]["derivatives"]) != 0:
                    for d in self.sort_derivatives(self.node[field_name]["derivatives"].keys()):
                        if len(d) == 1:
                            ipfiel_file.write("\n The derivative wrt direction %s is [ 0.00000E+00]:    % *.*f" % (d, precision+5, precision, self.node[field_name]["derivatives"][d][i,f]))
                        if len(d) == 2:
                            ipfiel_file.write("\n The derivative wrt directions %s & %s is [ 0.00000E+00]:    % *.*f" % (d[0], d[1], precision+5, precision, self.node[field_name]["derivatives"][d][i,f]))
                        if len(d) == 3:
                            ipfiel_file.write("\n The derivative wrt directions %s, %s & %s is [ 0.00000E+00]:    % *.*f" % (d[0], d[1], d[2], precision+5, precision, self.node[field_name]["derivatives"][d][i,f]))
            ipfiel_file.write("\n")
        
        ipfiel_file.close()
        
        return

    def export_ipelfd(self, file_path, file_name = None, group_name = None, field_name = None):
    # Note: Elements are assumed to have the same number of coordinates and 
    #       basis functions.

        if file_name == None:
            file_name = self.file_name
        if group_name == None:
            if self.group_name == None:
                group_name = self.file_name
            else:
                group_name = self.group_name
        if field_name == None:
            if "general" in self.elem:
                field_name = "general"
                print ("!!! Warning: No field name specified. Exporting general field to ipelfd !!!")
            else:
                print ("!!! Error exporting ipelfd: Define field to export !!!")
                return
        
        # Check that nodes have been defined
        try:
            len(self.node[field_name]["values"])        
            # Check that elem nodes correspond to node ids
            try:
                len(self.node["ids"])
                elem_nodes = np.unique(np.ndarray.flatten(self.elem["nodes"]))
                elem_nodes.sort()
                check = self.node["ids"] - elem_nodes
                if np.amax(check) != 0:
                    print ("!!! Warning: Node numbers in elements do not match node ids !!!")
            except:
                print ("!!! Warning: Node ids not defined !!!")
        except:
            print ("!!! Warning: Nodal values have not been defined !!!")
        
        # Check that elem dictionary is complete
        try:
            len(self.elem["nodes"])
            if "ids" not in self.elem:
                self.elem["ids"] = np.arange(np.shape(self.elem["nodes"])[0]) + 1
        except:
            print ("!!! Error exporting ipelem: Element nodes not defined !!!")
            return
            
        try:
            len(self.elem[field_name]["bases"])
        except:
            self.elem[field_name]["base_order"] = [1]
            self.elem[field_name]["bases"] = []
            for n in range(np.shape(self.node[field_name]["values"])[1]):
                self.elem[field_name]["bases"].append(1)
            print ("!!! Warning: Basis functions not defined. Using default basis function (1) !!!")

        # Ensure the elements have four nodes (triangle mesh)
        if np.shape(self.elem["nodes"])[1] == 3:
            self.elem["nodes"] = np.append(self.elem["nodes"], np.reshape(self.elem["nodes"][:,2], (len(self.elem["nodes"][:,2]), 1)), axis = 1)
            print ("!!! Warning: Only three nodes defined per element !!!")
        
        # Create file directory
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        
        # Write file header
        ipelfd_path = os.path.join(file_path, file_name + ".ipelfd")
        ipelfd_file = open(ipelfd_path, 'wb')
        ipelfd_file.write(" CMISS Version 2.1  ipelfd File Version 1\n")
        ipelfd_file.write(" Heading: %s\n\n" % group_name)
        
        # Write ipelfd file
        ipelfd_file.write(" The number of elements is [ %5i]: %5i\n\n" % (len(self.elem["ids"]), len(self.elem["ids"])))
        
        for i in range(len(self.elem["ids"])):
            ipelfd_file.write(" Element number [%5i]: %5i\n" % (self.elem["ids"][i], self.elem["ids"][i]))
            for b in range(len(self.elem[field_name]["bases"])):
                ipelfd_file.write(" The basis function type for field variable %i is [1]:  %i\n" % (b+1, self.elem[field_name]["bases"][b]))
            nodes = self.elem["nodes"][i,:]
            for b in range(len(self.elem[field_name]["base_order"])):
                ipelfd_file.write(" Enter the %i numbers for basis %i [prev]:" % (len(nodes), self.elem[field_name]["base_order"][b]))
                for n in range(len(nodes)):
                    ipelfd_file.write(" %5i" % nodes[n])
                ipelfd_file.write("\n")
            ipelfd_file.write("\n")
    
        ipelfd_file.close()
        
        return
        
    #%%#####################################################################%%#
    ## Read elem files ########################################################
    def read_ipelem(self, file_path, file_name = None):
    # Note: There are assumed to be no versions for nodes in the ipelem file.
    # Note: Elements are assumed to have the same number of coordinates and 
    #       basis functions.

        # Open ipelem file
        if file_name == None:
            file_name = self.file_name
            
        with open(os.path.join(file_path, file_name + ".ipelem"), 'r') as myfile:
            ipelem_file = myfile.read()
        ipelem_file = ipelem_file.split("\n")[1:]
        
        # Get group name
        if self.group_name == None or self.group_name == "":
            self.group_name = ipelem_file[0].split(":")[1].strip()
        ipelem_file = ipelem_file[2:]
        
        # Read file header
        elements = int(ipelem_file[0].split(":")[1].strip())
        ipelem_file = ipelem_file[2:]
        
        # Initiate dictionaries
        if "ids" in self.elem and len(self.elem["ids"]) != 0:
            print ("!!! Warning: overwriting elem ids with ipelem file !!!")
        self.elem["ids"] = np.zeros(elements, dtype=int)
        self.elem["coordinates"]["bases"] = []

        # Split remaining file into list of elements
        ipelem_file = "\n".join(ipelem_file).split("Element")
        if len(ipelem_file[0]) == 1:
            ipelem_file = ipelem_file[1:]

        # Read first elements and define bases
        temp = ipelem_file[0].split("\n")
        self.elem["ids"][0] = int(temp[0].split(":")[1].strip())
        temp = temp[1:]
        
        components = int(temp[0].split(":")[1].strip())

        temp = temp[1:]        
        for c in range(components):            
            self.elem["coordinates"]["bases"].append(int(temp[0].split(":")[1].strip()))
            temp = temp[1:]
        
        bases = []
        for i in range(len(temp)):
            if len(temp[i].strip()) != 0:
                base = [b for b in temp[i].split(":")[0] if b.isdigit()][1]
                bases.append(int(base))
                if "nodes" not in self.elem:
                    nodes = map(int, temp[i].split(":")[1].split())
                    self.elem["nodes"] = np.zeros((elements, len(nodes)), dtype=int)
                    self.elem["nodes"][0,:] = nodes
        self.elem["coordinates"]["base_order"] = bases[:]
        ipelem_file = ipelem_file[1:]

        # Read remaining elements
        for i in range(len(ipelem_file)):
            temp = ipelem_file[0].split("\n")
            self.elem["ids"][i+1] = int(temp[0].split(":")[1].strip())
            self.elem["nodes"][i+1,:] = map(int, temp[components+2].split(":")[1].split())
            ipelem_file = ipelem_file[1:]
    
        return    

    def read_exelem(self, file_path, file_name = None):
        # Open ipelem file        
        if file_name == None:
            file_name = self.file_name
        
        with open(os.path.join(file_path, file_name + ".exelem"), 'r') as myfile:
            exelem_file = myfile.read()
        exelem_file = exelem_file.split("\n")
        
        # Get group name
        if self.group_name == None or self.group_name == "":
            self.group_name = exelem_file[0].split(":")[1].strip()
        exelem_file = exelem_file[1:]

        # Split final shape section into header and elements
        exelem_file = "\n".join(exelem_file).split("Shape")[-1].split("Element:")
        header = exelem_file[0].split("\n")[1:]
        elements = exelem_file[1:]
        
        # Read header
        header = header
        interpolation_functions = int(header[0].split("=")[-1])
        header = header[1:]
        interpolations = {}
        for i in range(interpolation_functions):
            function = header[0].split(",")[0].strip()
            interpolations[function] = int(header[0].split(",")[-1].split("=")[-1])
            header = header[1:]
        
        nodes = int(header[0].split("=")[-1])
        header = header[1:]
        
        fields = int(header[0].split("=")[-1])
        header = header[1:]
        
        # Define components and interpolation functions
        for f in range(fields):
            field = header[0].split(")")[1].split(",")[0].strip()
            self.elem[field] = {}
            self.elem[field]["components"] = []
            self.elem[field]["interpolation"] = []
            components = int(header[0].split("=")[-1])
            header = header[1:]
            for c in range(components):
                self.elem[field]["components"].append(header[0].split()[0].split(".")[0])
                if len(self.elem[field]["interpolation"]) == 0:
                    interpolation = header[0].split()[1].split(",")[0]
                    self.elem[field]["interpolation"] = interpolation
                    self.elem[field]["scale_factors"] = interpolations[interpolation]
                header = header[(2 + nodes*3):]
                
        # Read elements
        if "ids" in self.elem and len(self.elem["ids"]) != 0:
            print ("!!! Warning: overwriting elem ids with exelem file !!!")                
        self.elem["ids"] = []
        self.elem["nodes"] = np.zeros([1,nodes], dtype=int)
        for e in range(len(elements)):
            temp = elements[0].split(":")
            self.elem["ids"].append(int(temp[0].split()[0]))
            self.elem["nodes"] = np.vstack((self.elem["nodes"], map(int, temp[-2].split()[:nodes])))
            elements = elements[1:]
        self.elem["nodes"] = self.elem["nodes"][1:,:]                    
        
        return
    
    ## Write elem files #######################################################
    def export_ipelem(self, file_path, file_name = None, group_name = None):
    # Note: There are assumed to be no versions for nodes in the ipelem file.
    # Note: Elements are assumed to have the same number of coordinates and 
    #       basis functions.

        if file_name == None:
            file_name = self.file_name
        if group_name == None:
            if self.group_name == None:
                group_name = self.file_name
            else:
                group_name = self.group_name
        
        # Check that nodes have been defined
        try:
            len(self.node["coordinates"]["values"])        
            # Check that elem nodes correspond to node ids
            try:
                len(self.node["ids"])
                elem_nodes = np.unique(np.ndarray.flatten(self.elem["nodes"]))
                elem_nodes.sort()
                check = self.node["ids"] - elem_nodes
                if np.amax(check) != 0:
                    print ("!!! Warning: Node numbers in elements do not match node ids !!!")
            except:
                print ("!!! Warning: Node ids not defined !!!")
        except:
            print ("!!! Warning: Nodal values have not been defined !!!")
        
        # Check that elem dictionary is complete
        try:
            len(self.elem["nodes"])
            if "ids" not in self.elem:
                self.elem["ids"] = np.arange(np.shape(self.elem["nodes"])[0]) + 1
        except:
            print ("!!! Error exporting ipelem: Element nodes not defined !!!")
            return
            
        try:
            len(self.elem["coordinates"]["bases"])
        except:
            self.elem["coordinates"]["base_order"] = [1]
            self.elem["coordinates"]["bases"] = []
            for n in range(np.shape(self.node["coordinates"]["values"])[1]):
                self.elem["coordinates"]["bases"].append(1)
            print ("!!! Warning: Basis functions not defined. Using default basis function (1) !!!")

        # Ensure the elements have four nodes (triangle mesh)
        if np.shape(self.elem["nodes"])[1] == 3:
            self.elem["nodes"] = np.append(self.elem["nodes"], np.reshape(self.elem["nodes"][:,2], (len(self.elem["nodes"][:,2]), 1)), axis = 1)
            print ("!!! Warning: Only three nodes defined per element !!!")
        
        # Create file directory
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        
        # Write file header
        ipelem_path = os.path.join(file_path, file_name + ".ipelem")
        ipelem_file = open(ipelem_path, 'wb')
        ipelem_file.write(" CMISS Version 2.1  ipelem File Version 2\n")
        ipelem_file.write(" Heading: %s\n\n" % group_name)
        
        # Write ipelem file
        ipelem_file.write(" The number of elements is [1]: %5i\n\n" % len(self.elem["ids"]))
        
        for i in range(len(self.elem["ids"])):
            ipelem_file.write(" Element number [%5i]: %5i\n" % (self.elem["ids"][i], self.elem["ids"][i]))
            ipelem_file.write(" The number of geometric Xj-coordinates is [3]: %i\n" % (len(self.elem["coordinates"]["bases"])))
            for b in range(len(self.elem["coordinates"]["bases"])):
                ipelem_file.write(" The basis function type for geometric variable %i is [1]:  %i\n" % (b+1, self.elem["coordinates"]["bases"][b]))
            nodes = self.elem["nodes"][i,:]
            ipelem_file.write(" Enter the %i global numbers for basis %i:" % (len(nodes), self.elem["coordinates"]["base_order"][0]))
            for n in range(len(nodes)):
                ipelem_file.write(" %5i" % nodes[n])
            ipelem_file.write("\n")
            if len(self.elem["coordinates"]["base_order"]) != 1:
                for b in np.arange(1, len(self.elem["coordinates"]["base_order"])):
                    ipelem_file.write(" Enter the %i numbers for basis %i [prev]:" % (len(nodes), self.elem["coordinates"]["base_order"][b]))
                    for n in range(len(nodes)):
                        ipelem_file.write(" %5i" % nodes[n])
                    ipelem_file.write("\n")
            ipelem_file.write("\n")
    
        ipelem_file.close()
        
        return

    def export_exelem(self, file_path, file_name = None, group_name = None):
    # Note: Elements are assumed to have the same number of coordinates (3) and 
    #       basis functions. Only linear lagrange (3 derivs) and cubic hermite
    #       (7 derivatives) are supported.
    # Note: Scale factors are all set to 1.0.

        if file_name == None:
            file_name = self.file_name
        if group_name == None:
            if self.group_name == None:
                group_name = self.file_name
            else:
                group_name = self.group_name
        
        # Check that nodes have been defined
        try:
            len(self.node["coordinates"]["values"])        
            # Check that elem nodes correspond to node ids
            try:
                len(self.node["ids"])
                elem_nodes = np.unique(np.ndarray.flatten(self.elem["nodes"]))
                elem_nodes.sort()
                check = self.node["ids"] - elem_nodes
                if np.amax(check) != 0:
                    print ("!!! Warning: Node numbers in elements do not match node ids !!!")
            except:
                print ("!!! Warning: Node ids not defined !!!")
        except:
            print ("!!! Error exporting exelem: Nodal values must be defined first !!!")
            return
        
        # Check that elem dictionary is complete
        try:
            len(self.elem["nodes"])
            if "ids" not in self.elem:
                self.elem["ids"] = np.arange(np.shape(self.elem["nodes"])[0]) + 1
        except:
            print ("!!! Error exporting exelem: Element nodes not defined !!!")
            return

        # Ensure the elements have four nodes (triangle mesh)
        if np.shape(self.elem["nodes"])[1] == 3:
            self.elem["nodes"] = np.append(self.elem["nodes"], np.reshape(self.elem["nodes"][:,2], (len(self.elem["nodes"][:,2]), 1)), axis = 1)
            print ("!!! Warning: Only three nodes defined per element !!!")

        # Determine interpolation functions
        nodes = np.shape(self.elem["nodes"])[1]
        if nodes == 2:
            dimensions = 1
        elif nodes == 4:
            dimensions = 2
        elif nodes == 8:
            dimensions = 3
                
        field_order = [key for key in self.node.keys() if key != "ids"]
        if "coordinates" in field_order:
            field_order.insert(0, field_order.pop(field_order.index("coordinates")))
        for field in field_order:
            if field not in self.elem:
                self.elem[field] = {}
            if "components" not in self.elem[field]:
                self.elem[field]["components"] = self.node[field]["components"]
            if dimensions == 1:
                if "1" in self.node[field]["derivatives"]:
                    self.elem[field]["interpolation"] = "c.Hermite"
                    self.elem[field]["scale_factors"] = 4
                else:
                    self.elem[field]["interpolation"] = "l.Lagrange"
                    self.elem[field]["scale_factors"] = 2
            if dimensions == 2:
                if all(d in self.node[field]["derivatives"] for d in ["1", "2"]):
                    self.elem[field]["interpolation"] = "c.Hermite*c.Hermite"
                    self.elem[field]["scale_factors"] = 16
                elif "1" in self.node[field]["derivatives"]:
                    self.elem[field]["interpolation"] = "c.Hermite*l.Lagrange"
                    self.elem[field]["scale_factors"] = 8
                elif "2" in self.node[field]["derivatives"]:
                    self.elem[field]["interpolation"] = "l.Lagrange*c.Hermite"
                    self.elem[field]["scale_factors"] = 8
                else:
                    self.elem[field]["interpolation"] = "l.Lagrange*l.Lagrange"
                    self.elem[field]["scale_factors"] = 4
            if dimensions == 3:
                if all(d in self.node[field]["derivatives"] for d in ["1", "2", "3"]):
                    self.elem[field]["interpolation"] = "c.Hermite*c.Hermite*c.Hermite"
                    self.elem[field]["scale_factors"] = 64
                elif all(d in self.node[field]["derivatives"] for d in ["1", "2"]):
                    self.elem[field]["interpolation"] = "c.Hermite*c.Hermite*l.Lagrange"
                    self.elem[field]["scale_factors"] = 32
                elif all(d in self.node[field]["derivatives"] for d in ["2", "3"]):
                    self.elem[field]["interpolation"] = "l.Lagrange*c.Hermite*c.Hermite"
                    self.elem[field]["scale_factors"] = 32
                elif all(d in self.node[field]["derivatives"] for d in ["1", "3"]):
                    self.elem[field]["interpolation"] = "c.Hermite*l.Lagrange*c.Hermite"
                    self.elem[field]["scale_factors"] = 32
                elif "1" in self.node[field]["derivatives"]:
                    self.elem[field]["interpolation"] = "c.Hermite*l.Lagrange*l.Lagrange"
                    self.elem[field]["scale_factors"] = 16
                elif "2" in self.node[field]["derivatives"]:
                    self.elem[field]["interpolation"] = "l.Lagrange*c.Hermite*l.Lagrange"
                    self.elem[field]["scale_factors"] = 16
                elif "3" in self.node[field]["derivatives"]:
                    self.elem[field]["interpolation"] = "l.Lagrange*l.Lagrange*c.Hermite"
                    self.elem[field]["scale_factors"] = 16
                else:
                    self.elem[field]["interpolation"] = "l.Lagrange*l.Lagrange*l.Lagrange"
                    self.elem[field]["scale_factors"] = 8
            
        interpolations = []
        scale_factors = []
        for field in field_order:
            if self.elem[field]["interpolation"] not in interpolations:
                interpolations.append(self.elem[field]["interpolation"])
            if self.elem[field]["scale_factors"] not in scale_factors:
                scale_factors.append(self.elem[field]["scale_factors"])
                
        # Create lists of lines and faces
        if dimensions > 1:           
            if dimensions == 2:
                lines = np.array([0,0])
                faces = np.array([0,0])
                lines_i = np.array([[1,2],[1,3],[2,4],[3,4]]) - 1
                faces_i = np.array([[1,3],[2,4],[1,2],[3,4]]) - 1
            elif dimensions == 3:
                lines = np.array([0,0])
                faces = np.array([0,0,0,0,0,0,0,0])
                lines_i = np.array([[1,2],[1,3],[1,5],[2,4],[2,6],[3,4],[3,7],[4,8],[5,6],[5,7],[6,8],[7,8]]) - 1
                faces_i = np.array([[1,5,3,7,1,3,5,7],
                                    [2,6,4,8,2,4,6,8],
                                    [1,2,5,6,1,5,2,6],
                                    [3,4,7,8,3,7,4,8],
                                    [1,3,2,4,1,2,3,4],
                                    [5,7,6,8,5,6,7,8]]) - 1

            for e in range(len(self.elem["ids"])):
                for l in range(np.shape(lines_i)[0]):
                    lines = np.vstack((lines, self.elem["nodes"][e,:][lines_i[l,:]]))
                for f in range(np.shape(faces_i)[0]):
                    faces = np.vstack((faces, self.elem["nodes"][e,:][faces_i[f,:]]))
    
            unique_l, unique_l_index = np.unique(lines[1:,:], return_index=True, axis=0)
            lines_u = lines[1:,:][np.sort(unique_l_index)]
            
            unique_f, unique_f_index = np.unique(faces[1:,:], return_index=True, axis=0)
            faces_u = faces[1:,:][np.sort(unique_f_index)]
        
        # Create file directory
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        
        # Write file header
        exelem_path = os.path.join(file_path, file_name + ".exelem")
        exelem_file = open(exelem_path, 'wb')
    
        exelem_file.write(" Group name: %s\n" % group_name)
        
        # Write ipelem file
        if dimensions > 1:
            exelem_file.write(" Shape. Dimension=1\n")
            for l in range(np.shape(lines_u)[0]):
                exelem_file.write(" Element: 0 0 %5i\n" % (l+1))
        if dimensions > 2:
            exelem_file.write(" Shape. Dimension=2\n")
            for f in range(np.shape(faces_u)[0]):
                exelem_file.write(" Element: 0 %5i 0\n" % (f+1))
                exelem_file.write("   Faces:\n")
                l1 = np.where((lines_u == faces_u[f,:2]).all(axis=1))[0][0]
                l2 = np.where((lines_u == faces_u[f,2:4]).all(axis=1))[0][0]
                l3 = np.where((lines_u == faces_u[f,4:6]).all(axis=1))[0][0]
                l4 = np.where((lines_u == faces_u[f,6:]).all(axis=1))[0][0]
                exelem_file.write("   0 0 %5i\n" % (l1+1))
                exelem_file.write("   0 0 %5i\n" % (l2+1))
                exelem_file.write("   0 0 %5i\n" % (l3+1))
                exelem_file.write("   0 0 %5i\n" % (l4+1))

        exelem_file.write(" Shape. Dimension=%i\n" % dimensions)
        exelem_file.write(" #Scale factor sets=%i\n" % len(scale_factors))
        for i in range(len(interpolations)):
            exelem_file.write("   %s, #Scale factors=%i\n" % (interpolations[i], scale_factors[i]))
        exelem_file.write(" #Nodes=%i\n" % nodes)
        exelem_file.write(" #Fields=%i\n" % len(field_order))
        f = 1
        interpolation = None
        for field in field_order:
            if field == "coordinates":
                field_type = "coordinate"
            else:
                field_type = "field"
            values = self.elem[field]["scale_factors"]/nodes
            components = len(self.elem[field]["components"])
            exelem_file.write(" %i) %s, %s, rectangular cartesian, #Components=%i\n" % (f, field, field_type, components))
            for c in range(components):
                if interpolation == None or interpolation == self.elem[field]["interpolation"]:
                    i = 1
                exelem_file.write("   %s.  %s, no modify, standard node based.\n" % (self.elem[field]["components"][c], self.elem[field]["interpolation"]))
                exelem_file.write("      #Nodes=%i\n" % nodes)
                for n in range(nodes):
                    exelem_file.write("      %i.  #Values=%i\n" % (n+1, values))
                    exelem_file.write("       Value indices:  ")
                    for v in range(values):
                        exelem_file.write("%4i" % (v+1))
                    exelem_file.write("\n       Scale factor indices:")
                    for v in range(values):
                        exelem_file.write("%4i" % (i))
                        i += 1
                    exelem_file.write("\n")
            interpolation = self.elem[field]["interpolation"][:]
            f += 1
        
        for e in range(len(self.elem["ids"])):
            nodes = self.elem["nodes"][e,:]
            exelem_file.write(" Element: %12i 0 0\n" % (e+1))
            if dimensions == 2:
                exelem_file.write("   Faces:\n")
                for f in range(np.shape(faces_i)[0]):
                    face = np.where((faces_u == nodes[faces_i[f,:]]).all(axis=1))[0][0]
                    exelem_file.write("     0 0 %5i\n" % (face+1))
            elif dimensions == 3:
                exelem_file.write("   Faces:\n")
                for f in range(np.shape(faces_i)[0]):
                    face = np.where((faces_u == nodes[faces_i[f,:]]).all(axis=1))[0][0]
                    exelem_file.write("     0 %5i 0\n" % (face+1))
            exelem_file.write("   Nodes:\n    ")
            for n in nodes:
                exelem_file.write("     %8i" % (n))
            exelem_file.write("\n   Scale factors:\n  ")
            for s in range(i-1):
                exelem_file.write("  1.0")
                if (s+1) % 5 == 0:
                    exelem_file.write("\n  ")
            exelem_file.write("\n")
        
        exelem_file.close()
            
        return

    #%%#####################################################################%%#
    ## Create default ipfibr file #############################################
    def create_ipfibr(self, file_path, file_name = None, group_name = None, precision = None, nodes = None):
    # Note: This function creates a default ipfibr file, based on nodes already
    #       defined.

        if file_name == None:
            file_name = self.file_name
        if group_name == None:
            if self.group_name == None:
                group_name = self.file_name
            else:
                group_name = self.group_name
        if precision == None:
            precision = self.precision
        
        # Check that nodes have been defined
        try:
            len(self.node["ids"])
        except:
            try:
                self.node["ids"] = np.arange(nodes) + 1
                self.node["coordinates"] = {}
                self.node["coordinates"]["components"] = ["x", "y", "z"]
                print ("!!! Warning: Node ids not defined. Using default node ids and variables !!!")
            except:
                print ("!!! Error exporting ipfibr: Node ids not defined !!!")
                return
        
        # Create file directory
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        
        # Open ipfibr_file file
        ipfibr_path = os.path.join(file_path, file_name + ".ipfibr")
        ipfibr_file = open(ipfibr_path, 'wb');
        
        # Write file header
        ipfibr_file.write(" CMISS Version 2.1  ipfibr File Version 3\n")
        ipfibr_file.write(" Heading: %s\n\n" % (group_name))
        
        ipfibr_file.write(" The number of fibre variables is [1]:%2i\n" % (len(self.node["coordinates"]["components"])))
        ipfibr_file.write(" Specify how fibre angle is defined [1]:\n")
        ipfibr_file.write("   (1) wrt Xi(1) coordinate\n")
        ipfibr_file.write("   (2) wrt Xi(2) coordinate\n    1\n")
        ipfibr_file.write(" Specify whether angles entered in [1]:\n")
        ipfibr_file.write("   (1) degrees\n")
        ipfibr_file.write("   (2) radians\n    1\n")
        
        ipfibr_file.write(" The number of nodes is [%5i]: %5i\n" % (len(self.node["ids"]), len(self.node["ids"])))
        ipfibr_file.write(" Do you want prompting for different versions of the fibre angle [N]? N\n")
        ipfibr_file.write(" Do you want prompting for different versions of the imbrication angle [N]? N\n")
        ipfibr_file.write(" Do you want prompting for different versions of the sheet angle [N]? N\n")
        ipfibr_file.write(" The number of derivatives for the fibre angle is [0]: 0\n")
        ipfibr_file.write(" The number of derivatives for the imbrication angle is [0]: 0\n")
        ipfibr_file.write(" The number of derivatives for the sheet angle is [0]: 0\n")
        
        # Write ipfibr file
        for i in range(np.shape(self.node["ids"])[0]):
            ipfibr_file.write("\n Node number [%5i]: %5i\n" % (self.node["ids"][i], self.node["ids"][i]))
            ipfibr_file.write(" The fibre angle is [ 0.00000E+00]:  0.00000E+00\n")
            ipfibr_file.write(" The imbrication angle is [ 0.00000E+00]:  0.00000E+00\n")
            ipfibr_file.write(" The sheet angle is [ 0.00000E+00]:  0.00000E+00\n")
        
        ipfibr_file.close()        

    ## Create default ipelfb file #############################################
    def create_ipelfb(self, file_path, file_name = None, group_name = None):
    # Note: This function creates a default ipelfb file, based on nodes and 
    #       elems already defined.

        if file_name == None:
            file_name = self.file_name
        if group_name == None:
            if self.group_name == None:
                group_name = self.file_name
            else:
                group_name = self.group_name
        
        # Check that elems have been defined
        try:
            len(self.elem["nodes"])
            try:
                len(self.elem["ids"])
            except:
                self.elem["ids"] = np.arange(np.shape(self.elem["nodes"])[0]) + 1
                print ("!!! Warning: Element ids not defined. Using default element ids !!!")
            # Check that elem nodes correspond to node ids
            try:
                len(self.node["ids"])
                elem_nodes = np.unique(np.ndarray.flatten(self.elem["nodes"]))
                elem_nodes.sort()
                check = self.node["ids"] - elem_nodes
                if np.amax(check) != 0:
                    print ("!!! Warning: Node numbers in elements do not match node ids !!!")
            except:
                print ("!!! Warning: Node ids not defined !!!")
        except:
            print ("!!! Error exporting ipeflb: Elements nodes not defined !!!")
            return
        
        # Check that elem dictionary is complete 
        try:
            len(self.elem["coordinates"]["bases"])
        except:
            self.elem["coordinates"]["base_order"] = [1]
            self.elem["coordinates"]["bases"] = []
            for n in range(np.shape(self.node["coordinates"]["values"])[1]):
                self.elem["coordinates"]["bases"].append(1)
            print ("!!! Warning: Basis functions not defined. Using default basis function (1) !!!")
        
        # Ensure the elements have four nodes (triangle mesh)
        if np.shape(self.elem["nodes"])[1] == 3:
            self.elem["nodes"] = np.append(self.elem["nodes"], np.reshape(self.elem["nodes"][:,2], (len(self.elem["nodes"][:,2]), 1)), axis = 1)
            print ("!!! Warning: Only three nodes defined per element !!!")
        
        # Create file directory
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        
        # Write file header
        ipelfb_path = os.path.join(file_path, file_name + ".ipelfb")
        ipelfb_file = open(ipelfb_path, 'wb')
        ipelfb_file.write(" CMISS Version 2.1  ipelfb File Version 1\n")
        ipelfb_file.write(" Heading: %s\n\n" % group_name)
        
        # Write ipelem file
        ipelfb_file.write(" The number of elements is [%5i]: %5i\n\n" % (len(self.elem["ids"]), len(self.elem["ids"])))
        
        for i in range(len(self.elem["ids"])):
            ipelfb_file.write(" Element number [%5i]: %5i\n" % (self.elem["ids"][i], self.elem["ids"][i]))
            ipelfb_file.write(" The basis function type for the fibre angle is [1]: %i\n" % (self.elem["coordinates"]["bases"][0]))
            ipelfb_file.write(" The basis function type for the imbrication angle is [1]: %i\n" % (self.elem["coordinates"]["bases"][0]))
            ipelfb_file.write(" The basis function type for the sheet angle is [1]: %i\n" % (self.elem["coordinates"]["bases"][0]))
            
            nodes = self.elem["nodes"][i,:]
            ipelfb_file.write(" Enter the %i numbers for basis %i [prev]:" % (len(nodes), self.elem["coordinates"]["base_order"][0]))
            for n in range(len(nodes)):
                ipelfb_file.write(" %5i" % nodes[n])
            ipelfb_file.write("\n")
            if len(self.elem["coordinates"]["base_order"]) != 1:
                for b in np.arange(1, len(self.elem["coordinates"]["base_order"])):
                    ipelfb_file.write(" Enter the %i numbers for basis %i [prev]:" % (len(nodes), self.elem["coordinates"]["base_order"][b]))
                    for n in range(len(nodes)):
                        ipelfb_file.write(" %5i" % nodes[n])
                    ipelfb_file.write("\n")
            ipelfb_file.write("\n")
        
        ipelfb_file.close()

#%%#########################################################################%%#
class cm_lists:    
    ## Create node list #######################################################
    @staticmethod
    def node(file_name, file_path):
        list_path = os.path.join(file_path, file_name + ".txt")
        with open(list_path, 'r') as myfile:
            list_file = myfile.read()

        list_file = list_file.split(",")
   
        node_id = np.array([]).astype(int)
        for i in list_file:
            node = i.split("..")
            if len(node) == 2:
                missing = np.arange(int(node[0]),int(node[1])+1).astype(int)
                node_id = np.append(node_id, missing)
            else:
                node_id = np.append(node_id, int(node[0]))
   
        return node_id

    ## Create elem list #######################################################
    @staticmethod
    def elem(file_name, file_path):
        elem_path = os.path.join(file_path, file_name + ".txt")
        try:
            load_elem = np.loadtxt(elem_path).astype(int)
            if len(np.shape(load_elem)) == 1:
                elem = np.zeros([1,4]).astype(int)
                elem[0:] = load_elem
            else:
                elem = load_elem
            
        except:
            try:
                with open(elem_path, 'r') as myfile:
                    list_file = myfile.readlines()
        
                elem_id = np.zeros([len(list_file),4]).astype(int)
                for i in range(len(list_file)):
                    elem_id[i:] = list_file[i].split(",")
                    
            except:
                print ("!! There is a problem with the element list format !!")
   
        return elem_id

#%%#########################################################################%%#
class cm_mesh(cm):
    '''This class creates a cm object'''
    ## Initiate class #########################################################
    def __init__(self, file_name, group_name = None):
        self.file_name = file_name.split(".")[0]
        self.group_name = group_name
        cm.__init__(self, self.file_name, group_name)
        
    ## Transform mesh #########################################################
    # Note: Derivatives are not modified.
    # Note: Only the coordinate field is transformed.
    def transform(self, tm):
        # Check for derivatives and fields:
        if "derivatives" in self.node["coordinates"]:
            print ("!!! Warning: Derivatives need to be updated !!!")        
        if len([key for key in self.node.keys() if key != "ids" if key != "coordinates"]) != 0:
            print ("!!! Warning: Only the coordinates fields is transformed !!!")
        
        self.node["coordinates"]["values"] = Transform.by_tm(self.node["coordinates"]["values"], tm)
        components = np.shape(self.node["coordinates"]["values"])[1]
        
        fields = [key for key in self.node.keys() if key != "ids" if key != "coordinates"]
        for field in fields:
            if np.shape(self.node[field]["values"])[1] == components:
                self.node[field]["values"] = Transform.by_tm(self.node[field]["values"], tm)
        
        return
    
    ## Remove border elements #############################################
    # Note: This function assumes the mesh is 2D.
    def remove_boundary_elements(self):        
        # Store original node ids
        node_ids = self.node["ids"][:]
        elems = self.elem["nodes"][:]
        
        # Remove nodes that exist in less than four elements
        fields = [key for key in self.node.keys() if key != "ids"]
        total_nodes = np.shape(node_ids)[0]
        for n in range(total_nodes):
            node = node_ids[n]
            e = np.where((elems == node))[0]
            if len(e) < 4:
                e = np.where((self.elem["nodes"] == node))[0]
                self.elem["ids"] = np.delete(self.elem["ids"], e)                
                self.elem["nodes"] = np.delete(self.elem["nodes"], e, axis=0)
                i = np.where((self.node["ids"] == node))[0]
                self.node["ids"] = np.delete(self.node["ids"], i)
                for field in fields:
                    self.node[field]["values"] = np.delete(self.node[field]["values"], i, axis=0)
                    for d in self.node[field]["derivatives"]:
                        self.node[field]["derivatives"][d] = np.delete(self.node[field]["derivatives"][d], i, axis=0)   
        
        # Renumber node and element ids
        node_ids = self.node["ids"]
        total_nodes = len(self.node["ids"])
        self.node["ids"] = np.arange(1, total_nodes+1)
        for n in range(total_nodes):
            indices = [self.elem["nodes"] == node_ids[n]][0]
            self.elem["nodes"][indices] = self.node["ids"][n]
        total_elems = len(self.elem["nodes"])
        self.elem["ids"] = np.arange(1, total_elems+1)        
    
        return
    
    ## Add layers to FE mesh ##################################################
    # Note: This function assumes the input mesh is 2D. 
    # Note: If the input mesh is non-linear, derivatives are removed (to 
    #       maintain interpolation, cross derivatives must be calculated).
    # Note: Only the coordinates field is maintained.
    def add_layers(self, thickness, layers):
        # Expand node arrays        
        final_node = self.node["ids"][-1]
        total_nodes = len(self.node["ids"])
        nodes = {}
        nodes["ids"] = np.pad(self.node["ids"], [(0,total_nodes*layers)], mode='constant')
        nodes["values"] = np.pad(self.node["coordinates"]["values"], [(0,total_nodes*layers),(0,0)], mode='constant')
        
        # Expand elem arrays
        elems = {}
        elems[1] = np.pad(self.elem["nodes"], [(0,0),(0,4)], mode='constant')
        for i in np.arange(2, layers+1):
            elems[i] = np.zeros(np.shape(elems[1]), dtype = int)
        
        # Extrude through nodes
        for n in range(total_nodes):
            n0 = self.node["coordinates"]["values"][n,:]
            d1 = self.node["coordinates"]["derivatives"]["1"][n,:]
            d2 = self.node["coordinates"]["derivatives"]["2"][n,:]
            d12 = np.cross(d1,d2)
            d12 = d12/norm(d12)
            
            i0 = n + 1
            for l in np.arange(1, layers+1):
                i1 = final_node*l + n + 1
                n1 = n0 + d12 * thickness * l
                nodes["ids"][i1-1] = i1
                nodes["values"][i1-1,:] = n1
        
                indices_1 = [elems[l] == nodes["ids"][i0-1]][0]
                indices_0 = np.zeros(np.shape(indices_1), dtype = bool)
                indices_0[:,:4] = indices_1[:,4:]
                indices_0[:,4:] = indices_1[:,:4]
                elems[l][indices_0] = i1
                if l < layers:
                    elems[l+1][indices_1] = i1
                i0 = i1

        # Collate elem arrays
        elements = np.array([[0,0,0,0,0,0,0,0]], dtype = int)
        for e in np.arange(1, len(elems)+1):
            elements = np.vstack((elements, elems[e]))
        elements = elements[1:,:]
        
        # Update cm object
        self.node["ids"] = nodes["ids"][:]
        self.node["coordinates"]["values"] = nodes["values"][:]
        self.elem["ids"] = np.arange(0, np.shape(elements)[0]) + 1
        self.elem["nodes"] = elements[:]
        
        # Remove derivatives and additional fields
        del self.node["coordinates"]["derivatives"]
        fields = [key for key in self.node.keys() if key != "ids" if key != "coordinates"]
        for field in fields:
            del self.node[field]
        
        return
        
#%%#########################################################################%%#
        
class cm_fit:
    ## Write fit files ########################################################
    @staticmethod   
    def create_files(file_name, group_name, template_name, bone, file_path):
        if not os.path.exists(file_path):
            os.makedirs(file_path)
        regions = 1
        data = 1000
        nodes = 2500
        dofs = 12
        elems = 2500
        faces = elems
        lines = elems*2
        
        ## Write fit com file #################################################
        fitcom_path = os.path.join(file_path, "1_fit.com")
        fitcom_file = open(fitcom_path, 'wb')
        fitcom_file.write("#!/usr/bin/perl\n")
        fitcom_file.write("############################################\n")
        fitcom_file.write("# Perform initial course fit\n")
        fitcom_file.write("fem def para;r;fit;\n")
        fitcom_file.write("fem def coord 3,1;\n")
        fitcom_file.write("fem def base;r;bicubic_hermite;\n\n")
        fitcom_file.write("fem def node;r;%s;\n" % (template_name))
        fitcom_file.write("fem def elem;r;%s;\n\n" % (template_name))
        fitcom_file.write("fem update node derivative 1 linear;\n")
        fitcom_file.write("fem update node derivative 2 linear;\n\n")      
        fitcom_file.write("#######################\n")
        fitcom_file.write("fem def data;r;%s;\n\n" % (file_name))
        fitcom_file.write("fem def xi;c close;\n")
        fitcom_file.write("fem list data error;\n")
        fitcom_file.write("fem def field;d;\n")
        fitcom_file.write("fem def elem;d field;\n\n")
        fitcom_file.write("#######################\n")
        fitcom_file.write("fem def fit;r;fit_1 geometry;\n\n")
        fitcom_file.write("fem fit;\n")
        fitcom_file.write("fem update node fit;\n")
        fitcom_file.write("fem def xi;c close;\n")
        fitcom_file.write("fem li data error;\n\n")
        fitcom_file.write("#######################\n")
        fitcom_file.write("fem exp data;%s as %s offset 0;\n\n" % (file_name, group_name))
        fitcom_file.write("fem def node;w;%s;\n" % (file_name))
        fitcom_file.write("fem def elem;w;%s;\n\n" % (file_name))	
        fitcom_file.write("############################################\n")
        fitcom_file.write("# Perform second close fit (256 elements)\n")
        fitcom_file.write("fem reallocate;\n")
        fitcom_file.write("fem def para;r;fit;\n")
        fitcom_file.write("fem def coord 3,1;\n")
        fitcom_file.write("fem def base;r;bicubic_hermite;\n\n")
        fitcom_file.write("fem def node;r;%s;\n" % (file_name))
        fitcom_file.write("fem def elem;r;%s;\n\n" % (file_name))
        fitcom_file.write("######################\n")
        fitcom_file.write("fem ref xi 1 at 0.5;\n")
        fitcom_file.write("fem ref xi 2 at 0.5;\n\n")
        fitcom_file.write("fem ref xi 1 at 0.5;\n")
        fitcom_file.write("fem ref xi 2 at 0.5;\n\n")
        if bone == "Rad": 
            fitcom_file.write("fem ref xi 1 at 0.5;\n")
            fitcom_file.write("fem ref xi 2 at 0.5;\n\n")       			
        fitcom_file.write("######################\n")
        fitcom_file.write("fem def data;r;%s;\n\n" % (file_name))
        fitcom_file.write("fem def xi;c close;\n")
        fitcom_file.write("fem list data error;\n")
        fitcom_file.write("fem def field;d;\n")
        fitcom_file.write("fem def elem;d field;\n\n")
        fitcom_file.write("######################\n")
        fitcom_file.write("fem def fit;r;fit_2 geometry;\n\n")
        fitcom_file.write("fem fit;\n")
        fitcom_file.write("fem update node fit;\n")
        fitcom_file.write("fem def xi;c close;\n")
        fitcom_file.write("fem li data error;\n\n")
        fitcom_file.write("######################\n")
        fitcom_file.write("fem def node;w;%s;\n" % (file_name))
        fitcom_file.write("fem def elem;w;%s;\n\n" % (file_name))		
        fitcom_file.write("###########################################\n")
        if bone == "Uln":
            fitcom_file.write("# Refine mesh\n")
        if bone == "Rad": 
            fitcom_file.write("# Refine and perform third close fit (1024 elements)\n")
        fitcom_file.write("fem reallocate;\n")
        fitcom_file.write("fem def para;r;fit;\n")
        fitcom_file.write("fem def coord 3,1;\n")
        fitcom_file.write("fem def base;r;bicubic_hermite;\n\n")
        fitcom_file.write("fem def node;r;%s;\n" % (file_name))
        fitcom_file.write("fem def elem;r;%s;\n\n" % (file_name))
        if bone == "Uln":
            fitcom_file.write("fem ref xi 1 at 0.5;\n")
            fitcom_file.write("fem ref xi 2 at 0.5;\n")
        fitcom_file.write("fem ref xi 1 at 0.5;\n")
        fitcom_file.write("fem ref xi 2 at 0.5;\n\n")            
        if bone == "Rad":
            fitcom_file.write("######################\n")
            fitcom_file.write("fem def data;r;%s;\n\n" % (file_name))
            fitcom_file.write("fem def xi;c close;\n")
            fitcom_file.write("fem list data error;\n")
            fitcom_file.write("fem def field;d;\n")
            fitcom_file.write("fem def elem;d field;\n\n")
            fitcom_file.write("######################\n")
            fitcom_file.write("fem def fit;r;fit_3 geometry;\n\n")
            fitcom_file.write("fem fit;\n")
            fitcom_file.write("fem update node fit;\n")
            fitcom_file.write("fem def xi;c close;\n")
            fitcom_file.write("fem li data error;\n\n")
            fitcom_file.write("######################\n")
        fitcom_file.write("fem def node;w;%s;\n" % (file_name))
        fitcom_file.write("fem def elem;w;%s;\n\n" % (file_name))       
        fitcom_file.write("fem exp node;%s as %s;\n" % (file_name, group_name))
        fitcom_file.write("fem exp elem;%s as %s;\n\n" % (file_name, group_name))	
        fitcom_file.write("###########################################\n")
        fitcom_file.write("fem quit;\n")
        fitcom_file.close()
        
        ## Write first ipfit file #############################################
        ipfit_path = os.path.join(file_path, "fit_1.ipfit")
        ipfit_file = open(ipfit_path, 'wb')
        ipfit_file.write("CMISS Version 1.21 ipfit File Version 2\n")
        ipfit_file.write(" Heading: %s\n\n" % (group_name)) 
        ipfit_file.write(" Specify whether problem is (1) linear or (2) nonlinear [1]: 1\n")
        ipfit_file.write(" Specify the number of fitting problems [1]: 1\n")
        ipfit_file.write(" Specify #geometric vars to be fitted for problem 1 [3]: 3\n")
        ipfit_file.write(" Specify field# to store var 1 of problem 1 [1]: 1\n")
        ipfit_file.write(" Specify field# to store var 2 of problem 1 [2]: 2\n")
        ipfit_file.write(" Specify field# to store var 3 of problem 1 [3]: 3\n")
        ipfit_file.write(" Specify geometric# to be fitted for var 1 of problem 1 [1]: 1\n")
        ipfit_file.write(" Specify geometric# to be fitted for var 2 of problem 1 [2]: 2\n")
        ipfit_file.write(" Specify geometric# to be fitted for var 3 of problem 1 [3]: 3\n")
        ipfit_file.write(" Enter smoothing type [0]:\n")
        ipfit_file.write("   (0) None\n")
        ipfit_file.write("   (1) Sobolev on field\n")
        ipfit_file.write("   (2) Sobolev on deviation from initial field\n")
        ipfit_file.write("   (3) Strain energy\n")
        ipfit_file.write("    1\n")
        if bone == "Rad":   
            ipfit_file.write(" Enter element #s/name [EXIT]: 1..4\n")
        elif bone == "Uln":
            ipfit_file.write(" Enter element #s/name [EXIT]: 1..8\n")
        else:
            print ("!! There is a problem with the bone name !!")
            return
        ipfit_file.write(" The 5 weights on derivs wrt Xi_1/_11/_2/_22/_12 are [prev]:\n")
        ipfit_file.write("	0.1 0.05 0.1 0.05 0.01\n")
        ipfit_file.write(" Enter element #s/name [EXIT]: 0\n\n")
        ipfit_file.write(" Do you want the global matrices stored as sparse matrices [Y]? Y\n")
        ipfit_file.write(" Do you want to calculate the sparsity pattern for the global matrices [Y]? Y\n\n")
        ipfit_file.write(" Do you want to enter the coupling coefficients [N]? N\n\n")
        ipfit_file.write(" Enter node #s/name to fix [EXIT]: 0\n")
        ipfit_file.write(" Specify type of linear solution procedure [1]:\n")
        ipfit_file.write("   (1)  LU Decomposition\n")
        ipfit_file.write("   (2)  Single Value Decomposition\n")
        ipfit_file.write("   (3)  Least Squares\n")
        ipfit_file.write("   (4)  Cholesky Decomposition\n")
        ipfit_file.write("   (5)  Jacobi Iteration\n")
        ipfit_file.write("   (6)  Succesive Over Relaxation\n")
        ipfit_file.write("   (7)  Incomplete LU Decomposition(0)\n")
        ipfit_file.write("   (8)  Incomplete LU Decomposition(1)\n")
        ipfit_file.write("   (9)  Conjugate Gradient\n")
        ipfit_file.write("   (10) Biconjugate Gradient Stabilised\n")
        ipfit_file.write("   (11) Generalised Minimum Residual\n")
        ipfit_file.write("   (12) Black Box Multigrid (BOXMG)\n")
        ipfit_file.write("   (13) Preconditioned CG from BOXMG\n")
        ipfit_file.write("   (14) Algebraic Multigrid (AMG1R6)\n")
        ipfit_file.write("    1\n")
        ipfit_file.write(" Do you want the solution matrix stored as a sparse matrix [Y]? Y\n")
        ipfit_file.write(" Specify the LU solver [2]:\n")
        ipfit_file.write("   (1) SuperLU\n")
        ipfit_file.write("   (2) Umfpack\n")
        ipfit_file.write("    2\n")
        ipfit_file.write(" Specify the pivot threshold [0.1D0]: 0.10000D+00\n")
        ipfit_file.write(" Specify option for linear solution [0]:\n")
        ipfit_file.write("   (0) No output\n")
        ipfit_file.write("   (1) Timing and error check\n")
        ipfit_file.write("   (2) & Solver output\n")
        ipfit_file.write("   (3) & Global solution matrices\n")
        ipfit_file.write("   (4) & Global stiffness matrices\n")
        ipfit_file.write("   (5) & Element matrices\n")
        ipfit_file.write("    0\n")
        ipfit_file.close()

        ## Write second ipfit file ############################################
        ipfit_path = os.path.join(file_path, "fit_2.ipfit")
        ipfit_file = open(ipfit_path, 'wb')
        ipfit_file.write("CMISS Version 1.21 ipfit File Version 2\n")
        ipfit_file.write(" Heading: %s\n\n" % (group_name)) 
        ipfit_file.write(" Specify whether problem is (1) linear or (2) nonlinear [1]: 1\n")
        ipfit_file.write(" Specify the number of fitting problems [1]: 1\n")
        ipfit_file.write(" Specify #geometric vars to be fitted for problem 1 [3]: 3\n")
        ipfit_file.write(" Specify field# to store var 1 of problem 1 [1]: 1\n")
        ipfit_file.write(" Specify field# to store var 2 of problem 1 [2]: 2\n")
        ipfit_file.write(" Specify field# to store var 3 of problem 1 [3]: 3\n")
        ipfit_file.write(" Specify geometric# to be fitted for var 1 of problem 1 [1]: 1\n")
        ipfit_file.write(" Specify geometric# to be fitted for var 2 of problem 1 [2]: 2\n")
        ipfit_file.write(" Specify geometric# to be fitted for var 3 of problem 1 [3]: 3\n")
        ipfit_file.write(" Enter smoothing type [0]:\n")
        ipfit_file.write("   (0) None\n")
        ipfit_file.write("   (1) Sobolev on field\n")
        ipfit_file.write("   (2) Sobolev on deviation from initial field\n")
        ipfit_file.write("   (3) Strain energy\n")
        ipfit_file.write("    1\n")
        if bone == "Rad":   
            ipfit_file.write(" Enter element #s/name [EXIT]: 1..256\n")
            ipfit_file.write(" The 5 weights on derivs wrt Xi_1/_11/_2/_22/_12 are [prev]:\n")
            ipfit_file.write("	0.1 0.05 0.1 0.05 0.01\n")
            ipfit_file.write(" Enter element #s/name [EXIT]: 3,5,7,17,19,21,23,65,67,69,71,81,83,85,172,174,176,186,188,190,192,234,236,238,240,250,252,254\n")
            ipfit_file.write(" The 5 weights on derivs wrt Xi_1/_11/_2/_22/_12 are [prev]:\n")
            ipfit_file.write("	1.0 1.00 0.1 0.05 0.01\n")
            ipfit_file.write(" Enter element #s/name [EXIT]: 2,9..10,33..34,41..42,88,95..96,119..120,127..130,137..138,161..162,169,215..216,223..224,247..248,255\n")
            ipfit_file.write(" The 5 weights on derivs wrt Xi_1/_11/_2/_22/_12 are [prev]:\n")
            ipfit_file.write("	0.1 0.05 1.0 1.00 0.01\n")
            ipfit_file.write(" Enter element #s/name [EXIT]: 1,87,170,256\n")
            ipfit_file.write(" The 5 weights on derivs wrt Xi_1/_11/_2/_22/_12 are [prev]:\n")
            ipfit_file.write("	1.0 1.00 1.0 1.00 1.00\n")
        elif bone == "Uln":
            ipfit_file.write(" Enter element #s/name [EXIT]: 1..128\n")
            ipfit_file.write(" The 5 weights on derivs wrt Xi_1/_11/_2/_22/_12 are [prev]:\n")
            ipfit_file.write("	0.1 0.05 0.1 0.05 0.01\n")
        else:
            print ("!! There is a problem with the bone name !!")
            return
        ipfit_file.write(" Enter element #s/name [EXIT]: 0\n\n")
        ipfit_file.write(" Do you want the global matrices stored as sparse matrices [Y]? Y\n")
        ipfit_file.write(" Do you want to calculate the sparsity pattern for the global matrices [Y]? Y\n\n")
        ipfit_file.write(" Do you want to enter the coupling coefficients [N]? N\n\n")
        ipfit_file.write(" Enter node #s/name to fix [EXIT]: 0\n")
        ipfit_file.write(" Specify type of linear solution procedure [1]:\n")
        ipfit_file.write("   (1)  LU Decomposition\n")
        ipfit_file.write("   (2)  Single Value Decomposition\n")
        ipfit_file.write("   (3)  Least Squares\n")
        ipfit_file.write("   (4)  Cholesky Decomposition\n")
        ipfit_file.write("   (5)  Jacobi Iteration\n")
        ipfit_file.write("   (6)  Succesive Over Relaxation\n")
        ipfit_file.write("   (7)  Incomplete LU Decomposition(0)\n")
        ipfit_file.write("   (8)  Incomplete LU Decomposition(1)\n")
        ipfit_file.write("   (9)  Conjugate Gradient\n")
        ipfit_file.write("   (10) Biconjugate Gradient Stabilised\n")
        ipfit_file.write("   (11) Generalised Minimum Residual\n")
        ipfit_file.write("   (12) Black Box Multigrid (BOXMG)\n")
        ipfit_file.write("   (13) Preconditioned CG from BOXMG\n")
        ipfit_file.write("   (14) Algebraic Multigrid (AMG1R6)\n")
        ipfit_file.write("    1\n")
        ipfit_file.write(" Do you want the solution matrix stored as a sparse matrix [Y]? Y\n")
        ipfit_file.write(" Specify the LU solver [2]:\n")
        ipfit_file.write("   (1) SuperLU\n")
        ipfit_file.write("   (2) Umfpack\n")
        ipfit_file.write("    2\n")
        ipfit_file.write(" Specify the pivot threshold [0.1D0]: 0.10000D+00\n")
        ipfit_file.write(" Specify option for linear solution [0]:\n")
        ipfit_file.write("   (0) No output\n")
        ipfit_file.write("   (1) Timing and error check\n")
        ipfit_file.write("   (2) & Solver output\n")
        ipfit_file.write("   (3) & Global solution matrices\n")
        ipfit_file.write("   (4) & Global stiffness matrices\n")
        ipfit_file.write("   (5) & Element matrices\n")
        ipfit_file.write("    0\n")
        ipfit_file.close()

        ## Write third ipfit file ############################################
        if bone == "Rad":
            ipfit_path = os.path.join(file_path, "fit_3.ipfit")
            ipfit_file = open(ipfit_path, 'wb')
            ipfit_file.write("CMISS Version 1.21 ipfit File Version 2\n")
            ipfit_file.write(" Heading: %s\n\n" % (group_name)) 
            ipfit_file.write(" Specify whether problem is (1) linear or (2) nonlinear [1]: 1\n")
            ipfit_file.write(" Specify the number of fitting problems [1]: 1\n")
            ipfit_file.write(" Specify #geometric vars to be fitted for problem 1 [3]: 3\n")
            ipfit_file.write(" Specify field# to store var 1 of problem 1 [1]: 1\n")
            ipfit_file.write(" Specify field# to store var 2 of problem 1 [2]: 2\n")
            ipfit_file.write(" Specify field# to store var 3 of problem 1 [3]: 3\n")
            ipfit_file.write(" Specify geometric# to be fitted for var 1 of problem 1 [1]: 1\n")
            ipfit_file.write(" Specify geometric# to be fitted for var 2 of problem 1 [2]: 2\n")
            ipfit_file.write(" Specify geometric# to be fitted for var 3 of problem 1 [3]: 3\n")
            ipfit_file.write(" Enter smoothing type [0]:\n")
            ipfit_file.write("   (0) None\n")
            ipfit_file.write("   (1) Sobolev on field\n")
            ipfit_file.write("   (2) Sobolev on deviation from initial field\n")
            ipfit_file.write("   (3) Strain energy\n")
            ipfit_file.write("    1\n")
            ipfit_file.write(" Enter element #s/name [EXIT]: 1..1024\n")
            ipfit_file.write(" The 5 weights on derivs wrt Xi_1/_11/_2/_22/_12 are [prev]:\n")
            ipfit_file.write("	0.1 0.05 0.1 0.05 0.01\n")
            ipfit_file.write(" Enter element #s/name [EXIT]: 0\n\n")
            ipfit_file.write(" Do you want the global matrices stored as sparse matrices [Y]? Y\n")
            ipfit_file.write(" Do you want to calculate the sparsity pattern for the global matrices [Y]? Y\n\n")
            ipfit_file.write(" Do you want to enter the coupling coefficients [N]? N\n\n")            
            ipfit_file.write(" Enter node #s/name to fix [EXIT]: 7..9,24..25,78..81,282..289,1074..1089\n")
            ipfit_file.write(" Are any variables for variable 1 of fit variable 1 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 1 of fit variable 1 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 1 of fit variable 1 derivative 2 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 1 of fit variable 1 derivative 3 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 1 of fit variable 1 derivative 4 fixed [N]?: y\n")
            ipfit_file.write(" Are any variables for variable 2 of fit variable 1 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 2 of fit variable 1 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 2 of fit variable 1 derivative 2 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 2 of fit variable 1 derivative 3 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 2 of fit variable 1 derivative 4 fixed [N]?: y\n")
            ipfit_file.write(" Are any variables for variable 3 of fit variable 1 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 3 of fit variable 1 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 3 of fit variable 1 derivative 2 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 3 of fit variable 1 derivative 3 fixed [N]?: y\n")
            ipfit_file.write(" Is variable 3 of fit variable 1 derivative 4 fixed [N]?: y\n\n")            
            ipfit_file.write(" Enter node #s/name to fix [EXIT]: 0\n")
            ipfit_file.write(" Specify type of linear solution procedure [1]:\n")
            ipfit_file.write("   (1)  LU Decomposition\n")
            ipfit_file.write("   (2)  Single Value Decomposition\n")
            ipfit_file.write("   (3)  Least Squares\n")
            ipfit_file.write("   (4)  Cholesky Decomposition\n")
            ipfit_file.write("   (5)  Jacobi Iteration\n")
            ipfit_file.write("   (6)  Succesive Over Relaxation\n")
            ipfit_file.write("   (7)  Incomplete LU Decomposition(0)\n")
            ipfit_file.write("   (8)  Incomplete LU Decomposition(1)\n")
            ipfit_file.write("   (9)  Conjugate Gradient\n")
            ipfit_file.write("   (10) Biconjugate Gradient Stabilised\n")
            ipfit_file.write("   (11) Generalised Minimum Residual\n")
            ipfit_file.write("   (12) Black Box Multigrid (BOXMG)\n")
            ipfit_file.write("   (13) Preconditioned CG from BOXMG\n")
            ipfit_file.write("   (14) Algebraic Multigrid (AMG1R6)\n")
            ipfit_file.write("    1\n")
            ipfit_file.write(" Do you want the solution matrix stored as a sparse matrix [Y]? Y\n")
            ipfit_file.write(" Specify the LU solver [2]:\n")
            ipfit_file.write("   (1) SuperLU\n")
            ipfit_file.write("   (2) Umfpack\n")
            ipfit_file.write("    2\n")
            ipfit_file.write(" Specify the pivot threshold [0.1D0]: 0.10000D+00\n")
            ipfit_file.write(" Specify option for linear solution [0]:\n")
            ipfit_file.write("   (0) No output\n")
            ipfit_file.write("   (1) Timing and error check\n")
            ipfit_file.write("   (2) & Solver output\n")
            ipfit_file.write("   (3) & Global solution matrices\n")
            ipfit_file.write("   (4) & Global stiffness matrices\n")
            ipfit_file.write("   (5) & Element matrices\n")
            ipfit_file.write("    0\n")
            ipfit_file.close()
        
        ## Write ippara file ##################################################
        ippara_path = os.path.join(file_path, "fit.ippara")
        ippara_file = open(ippara_path, 'wb')
        ippara_file.write("CMISS Version 2.1  ippara File Version 1\n")
        ippara_file.write(" Heading:\n\n")
        ippara_file.write(" Max# auxiliary parameters          (NAM)[1]:         5\n")
        ippara_file.write(" Max# basis functions               (NBM)[1]:        42\n")
        ippara_file.write(" Max# var. types for a  dep. var.   (NCM)[1]:         2\n")
        ippara_file.write(" Max# data points                   (NDM)[1]:{:>10}\n".format(data))
        ippara_file.write(" Max# elements                      (NEM)[1]:{:>10}\n".format(elems*regions))
        ippara_file.write(" Max# elements in a region       (NE_R_M)[1]:{:>10}\n".format(elems))
        ippara_file.write(" Max# adjacent elements in Xi      (NEIM)[1]:       120\n")
        ippara_file.write(" Max# global face segments          (NFM)[1]:{:>10}\n".format(faces*regions))
        ippara_file.write(" Max# faces in a region          (NF_R_M)[1]:{:>10}\n".format(faces))
        ippara_file.write(" Max# local Voronoi faces         (NFVCM)[1]:         6\n")
        ippara_file.write(" Max# Gauss points per element      (NGM)[1]:        81\n")
        ippara_file.write(" Max# dependent variables           (NHM)[1]:         6\n")
        ippara_file.write(" Max# local Xi coordinates          (NIM)[1]:         3\n")
        ippara_file.write(" Max# global reference coordinates  (NJM)[1]:        15\n")
        ippara_file.write(" Max# derivatives per variable      (NKM)[1]:         8\n")
        ippara_file.write(" Max# global line segments          (NLM)[1]:{:>10}\n".format(lines*regions))
        ippara_file.write(" Max# lines in a region          (NL_R_M)[1]:{:>10}\n".format(lines))
        ippara_file.write(" Max# material parameters           (NMM)[1]:        35\n")
        ippara_file.write(" Max# element nodes                 (NNM)[1]:        64\n")
        ippara_file.write(" Max# degrees of freedom            (NOM)[1]:{:>10}\n".format(nodes*dofs*regions))
        ippara_file.write(" Max# global nodes                  (NPM)[1]:{:>10}\n".format(nodes*regions))
        ippara_file.write(" Max# global nodes in a region   (NP_R_M)[1]:{:>10}\n".format(nodes))
        ippara_file.write(" Max# global grid points            (NQM)[1]:      4225\n")
        ippara_file.write(" Max# signal sets                  (NSSM)[1]:         1\n")
        ippara_file.write(" Max# grid degrees of freedom      (NYQM)[1]:         1\n")
        ippara_file.write(" Max# regions                       (NRM)[1]:{:>10}\n".format(regions))
        ippara_file.write(" Max# element dofs per variable     (NSM)[1]:        64\n")
        ippara_file.write(" Max# face dofs per variable       (NSFM)[1]:        16\n")
        ippara_file.write(" Max# eigenvalues                   (NTM)[1]:        20\n")
        ippara_file.write(" Max# time samples                 (NTSM)[1]:       200\n")
        ippara_file.write(" Max# derivatives up to 2nd order   (NUM)[1]:        11\n")
        ippara_file.write(" Max# Voronoi boundary nodes      (NVCBM)[1]:       100\n")
        ippara_file.write(" Max# Voronoi cells                (NVCM)[1]:       200\n")
        ippara_file.write(" Max# versions of a variable        (NVM)[1]:        16\n")
        ippara_file.write(" Max# workstations                  (NWM)[1]:         3\n")
        ippara_file.write(" Max# problem types                 (NXM)[1]:         3\n")
        ippara_file.write(" Max# mesh dofs                     (NYM)[1]:{:>10}\n".format(nodes*dofs*regions))
        ippara_file.write(" Max# mesh dofs in a region      (NY_R_M)[1]:{:>10}\n".format(nodes*dofs))
        ippara_file.write(" Max# dimension of GD           (NZ_GD_M)[1]:    500000\n")
        ippara_file.write(" Max# dimension of GK           (NZ_GK_M)[1]:  20000000\n")
        ippara_file.write(" Max# dimension of GKK         (NZ_GKK_M)[1]:  10969344\n")
        ippara_file.write(" Max# dimension of GM           (NZ_GM_M)[1]:   1600000\n")
        ippara_file.write(" Max# dimension of GMM         (NZ_GMM_M)[1]:   1600000\n")
        ippara_file.write(" Max# dimension of GQ           (NZ_GQ_M)[1]:    500000\n")
        ippara_file.write(" Max# dimension of ISC_GD      (NISC_GDM)[1]:    500000\n")
        ippara_file.write(" Max# dimension of ISR_GD      (NISR_GDM)[1]:      4000\n")
        ippara_file.write(" Max# dimension of ISC_GK      (NISC_GKM)[1]:   5000000\n")
        ippara_file.write(" Max# dimension of ISR_GK      (NISR_GKM)[1]:   1000000\n")
        ippara_file.write(" Max# dimension of ISC_GKK    (NISC_GKKM)[1]:    500000\n")
        ippara_file.write(" Max# dimension of ISR_GKK    (NISR_GKKM)[1]:   1000000\n")
        ippara_file.write(" Max# dimension of ISC_GM      (NISC_GMM)[1]:   5000000\n")
        ippara_file.write(" Max# dimension of ISR_GM      (NISR_GMM)[1]:      4000\n")
        ippara_file.write(" Max# dimension of ISC_GMM    (NISC_GMMM)[1]:    500000\n")
        ippara_file.write(" Max# dimension of ISR_GMM    (NISR_GMMM)[1]:      4000\n")
        ippara_file.write(" Max# dimension of ISC_GQ      (NISC_GQM)[1]:    500000\n")
        ippara_file.write(" Max# dimension of ISR_GQ      (NISR_GQM)[1]:      4000\n")
        ippara_file.write(" Max# size of Minos arrays    (NZ_MINOSM)[1]:       125\n")
        ippara_file.write(" Max# basis function families      (NBFM)[1]:        13\n")
        ippara_file.write(" Max# nonlin. optim.n constraints  (NCOM)[1]:         0\n")
        ippara_file.write(" Max# data points in one element   (NDEM)[1]:      4000\n")
        ippara_file.write(" Max# dipoles in a region      (NDIPOLEM)[1]:        20\n")
        ippara_file.write(" Max# time points for a dipole (NDIPTIMM)[1]:       500\n")
        ippara_file.write(" Max# elements along a line        (NELM)[1]:        16\n")
        ippara_file.write(" Max# elements a node can be in    (NEPM)[1]:        20\n")
        ippara_file.write(" Max# segments                  (NGRSEGM)[1]:         6\n")
        ippara_file.write(" Max# variables per grid point     (NIQM)[1]:         6\n")
        ippara_file.write(" Max# cell state variables        (NIQSM)[1]:         0\n")
        ippara_file.write(" Max# variables for fibre extens(NIFEXTM)[1]:         8\n")
        ippara_file.write(" Max# variables per mesh dof       (NIYM)[1]:        16\n")
        ippara_file.write(" Max# variables / mesh dof(fix) (NIYFIXM)[1]:         5\n")
        ippara_file.write(" Max# vars. at each gauss point   (NIYGM)[1]:         6\n")
        ippara_file.write(" Max# vars. at face gauss points (NIYGFM)[1]:         0\n")
        ippara_file.write(" Max# linear optimis.n constraints (NLCM)[1]:         1\n")
        ippara_file.write(" Max# auxiliary grid parameters   (NMAQM)[1]:         6\n")
        ippara_file.write(" Max# cell material parameters     (NMQM)[1]:         0\n")
        ippara_file.write(" Max# optimisation variables       (NOPM)[1]:       125\n")
        ippara_file.write(" Max size fractal tree order array (NORM)[1]:        20\n")
        ippara_file.write(" Max# soln dofs for mesh dof       (NOYM)[1]:         1\n")
        ippara_file.write(" Max# domain nodes for BE problems (NPDM)[1]:         1\n")
        ippara_file.write(" Max# grid points per element      (NQEM)[1]:        81\n")
        ippara_file.write(" Max# non-zeros in grid matrix row (NQGM)[1]:        22\n")
        ippara_file.write(" Max# cell integer variables       (NQIM)[1]:         1\n")
        ippara_file.write(" Max# cell real variables          (NQRM)[1]:         1\n")
        ippara_file.write(" Max# spatial var cell int vars  (NQISVM)[1]:         1\n")
        ippara_file.write(" Max# spatial var cell real vars (NQRSVM)[1]:         1\n")
        ippara_file.write(" Max# number of grid schemes      (NQSCM)[1]:         9\n")
        ippara_file.write(" Max# cell variants                (NQVM)[1]:         0\n")
        ippara_file.write(" Max# rows and columns (sb 2)      (NRCM)[1]:         2\n")
        ippara_file.write(" Max# optimisation residuals       (NREM)[1]:       400\n")
        ippara_file.write(" Max# time points          (NTIMEPOINTSM)[1]:         1\n")
        ippara_file.write(" Max# time variables         (NTIMEVARSM)[1]:         1\n")
        ippara_file.write(" Max# mesh dofs for soln dof       (NYOM)[1]:         4\n")
        ippara_file.write(" Max# rows in a problem          (NYROWM)[1]:{:>10}\n".format(nodes*dofs*regions))
        ippara_file.write(" Max image cell array dimension (NIMAGEM)[1]:         1\n")
        ippara_file.write(" Size of transfer matrix  (NY_TRANSFER_M)[1]:       200\n")
        ippara_file.write(" Max# mesh dofs map to 1 mesh dof  (NYYM)[1]:        10\n")
        ippara_file.write(" USE_BEM       (0 or 1)[1]: 1\n")
        ippara_file.write(" USE_CELL      (0 or 1)[1]: 0\n")
        ippara_file.write(" USE_DATA      (0 or 1)[1]: 1\n")
        ippara_file.write(" USE_DIPOLE    (0 or 1)[1]: 1\n")
        ippara_file.write(" USE_GAUSS_PT_MATERIALS  (0 or 1)[0]: 0\n")
        ippara_file.write(" USE_GRAPHICS  (0 or 1)[1]: 1\n")
        ippara_file.write(" USE_GRID      (0 or 1)[1]: 1\n")
        ippara_file.write(" USE_LUNG      (0 or 1)[1]: 1\n")
        ippara_file.write(" USE_MAGNETIC  (0 or 1)[0]: 0\n")
        ippara_file.write(" USE_MAPS      (0 or 1)[0]: 1\n")
        ippara_file.write(" USE_MINOS     (0 or 1)[1]: 0\n")
        ippara_file.write(" USE_NLSTRIPE  (0 or 1)[1]: 1\n")
        ippara_file.write(" USE_NONLIN    (0 or 1)[1]: 1\n")
        ippara_file.write(" USE_NPSOL     (0 or 1)[1]: 1\n")
        ippara_file.write(" USE_SPARSE    (0 or 1)[1]: 1\n")
        ippara_file.write(" USE_TRANSFER  (0 or 1)[1]: 1\n")
        ippara_file.write(" USE_TRIANGLE  (0 or 1)[1]: 0\n")
        ippara_file.write(" USE_VORONOI   (0 or 1)[1]: 0\n")
        ippara_file.write(" USE_TIME      (0 or 1)[1]: 1\n")             
        ippara_file.close()

        ## Write ipbase files #################################################
        ipbase_path = os.path.join(file_path, "bicubic_hermite.ipbase")
        ipbase_file = open(ipbase_path, 'wb')
        ipbase_file.write("CMISS Version 1.21 ipbase File Version 2\n")
        ipbase_file.write(" Heading:\n\n")
        ipbase_file.write(" Enter the number of types of basis function [1]:  1\n\n")
        ipbase_file.write(" For basis function type 1 the type of nodal interpolation is [1]:\n")
        ipbase_file.write("   (0) Auxiliary basis only\n")
        ipbase_file.write("   (1) Lagrange/Hermite tensor prod\n")
        ipbase_file.write("   (2) Simplex/Serendipity/Sector\n")
        ipbase_file.write("   (3) B-spline tensor product\n")
        ipbase_file.write("   (4) Fourier Series/Lagrange/Hermite tensor prod\n")
        ipbase_file.write("   (5) Boundary Element Lagrange/Hermite tensor pr.\n")
        ipbase_file.write("   (6) Boundary Element Simplex/Serendipity/Sector\n")
        ipbase_file.write("   (7) Extended Lagrange (multigrid collocation)\n")
        ipbase_file.write("    1\n")
        ipbase_file.write(" Enter the number of Xi-coordinates [1]: 2\n\n")
        ipbase_file.write(" The interpolant in the Xi(1) direction is [1]:\n")
        ipbase_file.write("   (1) Linear Lagrange\n")
        ipbase_file.write("   (2) Quadratic Lagrange\n")
        ipbase_file.write("   (3) Cubic Lagrange\n")
        ipbase_file.write("   (4) Quadratic Hermite\n")
        ipbase_file.write("   (5) Cubic Hermite\n")
        ipbase_file.write("    5\n")
        ipbase_file.write(" Enter the number of Gauss points in the Xi(1) direction [3]: 3\n\n")
        ipbase_file.write(" The interpolant in the Xi(2) direction is [1]:\n")
        ipbase_file.write("   (1) Linear Lagrange\n")
        ipbase_file.write("   (2) Quadratic Lagrange\n")
        ipbase_file.write("   (3) Cubic Lagrange\n")
        ipbase_file.write("   (4) Quadratic Hermite\n")
        ipbase_file.write("   (5) Cubic Hermite\n")
        ipbase_file.write("    5\n")
        ipbase_file.write(" Enter the number of Gauss points in the Xi(2) direction [3]: 3\n")
        ipbase_file.write(" Do you want to set cross derivatives to zero [N]? N\n")
        ipbase_file.write(" Enter the node position indices [11211222]:  1 1 2 1 1 2 2 2\n")
        ipbase_file.write(" Enter the derivative order indices [11211222]:  1 1 2 1 1 2 2 2\n")
        ipbase_file.write(" Enter the number of auxiliary element parameters [0]:  0\n\n")
        ipbase_file.write(" For basis function type 1 scale factors are [6]:\n")
        ipbase_file.write("   (1) Unit\n")
        ipbase_file.write("   (2) Read in - Element based\n")
        ipbase_file.write("   (3) Read in - Node based\n")
        ipbase_file.write("   (4) Calculated from angle change\n")
        ipbase_file.write("   (5) Calculated from arc length\n")
        ipbase_file.write("   (6) Calculated from arithmetic mean arc length\n")
        ipbase_file.write("   (7) Calculated from harmonic mean arc length\n")
        ipbase_file.write("    6\n")
        ipbase_file.close()
        
        ## Write view com files ###############################################
        viewcom_path = os.path.join(file_path, "0_view.com")
        viewcom_file = open(viewcom_path, 'wb')
        viewcom_file.write("#!/usr/bin/perl\n")
        viewcom_file.write("$i = 10000\n")
        viewcom_file.write("$int = 10000\n\n")
        viewcom_file.write("###########################\n")
        viewcom_file.write("sub load_mesh\n")
        viewcom_file.write("{\n${file} = \"$_[0]\";\n")
        viewcom_file.write("${group} = \"$_[1]\";\n")
        viewcom_file.write("${mat} = \"$_[2]\";\n\n")
        viewcom_file.write("gfx read node ${file}.exnode;\n")
        viewcom_file.write("gfx read elem ${file}.exelem;\n\n")
        viewcom_file.write("$offset = $i;\n")
        viewcom_file.write("gfx change_id group ${group} node_offset $offset line_offset $offset face_offset $offset element_offset $offset;\n\n")
        viewcom_file.write("gfx mod g_elem ${group} node_points glyph sphere general size \"0.5*0.5*0.5\" centre 0,0,0 font default select_on invisible material ${mat} selected_material default_selected;\n")
        viewcom_file.write("gfx mod g_elem ${group} cylinders constant_radius 0.05 select_on material ${mat} selected_material default_selected render_shaded;\n")
        viewcom_file.write("gfx mod g_elem ${group} surfaces select_on material ${mat} data xi spectrum default selected_material default_selected render_shaded;\n")
        viewcom_file.write("# gfx mod g_elem ${group} element_points glyph diamond general size \"1*1*1\" centre 0,0,0 font default use_elements cell_centres discretization \"1*1*1\" native_discretization NONE select_on material default selected_material default_selected;\n\n")
        viewcom_file.write("$i = $i + $int;\n}\n\n")
        viewcom_file.write("sub load_data\n")
        viewcom_file.write("{\n${file} = \"$_[0]\";\n")
        viewcom_file.write("${group} = \"$_[1]\";\n")
        viewcom_file.write("${mat} = \"$_[2]\";\n\n")
        viewcom_file.write("gfx read data ${file}.exdata;\n\n")
        viewcom_file.write("$offset = $i;\n")
        viewcom_file.write("gfx change_id group ${group} data_offset $offset;\n\n")
        viewcom_file.write("gfx mod g_elem ${group} data_points glyph sphere general size \"0.4*0.4*0.4\" centre 0,0,0 font default select_on material ${mat} selected_material default_selected;\n\n")
        viewcom_file.write("$i = $i + $int;\n}\n\n")
        viewcom_file.write("###########################\n\n")
        viewcom_file.write("load_mesh(\"%s\", \"%s\", \"red\")\n" % (template_name, template_name))
        viewcom_file.write("load_mesh(\"%s\", \"%s\", \"gold\")\n" % (file_name, group_name))
        viewcom_file.write("load_data(\"%s\", \"%s\", \"silver\")\n" % (file_name, group_name))
        viewcom_file.write("load_data(\"%s\", \"%s\", \"muscle\")\n\n" % (file_name + "_all", group_name))
        viewcom_file.write("###########################\n")
        viewcom_file.write("gfx cre win;\n\n")
        viewcom_file.write("gfx cre axes len 200;\n")
        viewcom_file.write("gfx dra axes;\n")
        viewcom_file.close()
        
        return
    
#%%#########################################################################%%#
