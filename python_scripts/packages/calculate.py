# -*- coding: utf-8 -*-
"""
Created by Desney Greybe, 2018.

"""
import os
import numpy as np
import csv

#%%#########################################################################%%#
class Calculate:
    ## Calculate rms distance #################################################
    @staticmethod
    def rms_distance(data_1, data_2):
        # Calculate the rms distance between two data sets
        diff = data_1 - data_2
        rms = np.linalg.norm(diff, axis=1)

        return rms

    ## Calculate common statistics ############################################
    @staticmethod
    def statistics_2(data, percentiles):
        mean = np.mean(data)
        stdev = np.std(data)
        minimum = np.min(data)
        maximum = np.max(data)
        percentile = np.percentile(data, percentiles)

        return np.hstack(([mean, stdev, minimum, maximum], percentile))
    
    @staticmethod
    def statistics(data, percentiles = None):
        stats = {}
        
        # Averages
        stats["RMSE"] = np.sqrt(np.mean(data**2))
        stats["MAE"] = np.mean(abs(data))
        stats["Mean"] = np.mean(data)
        
        # Spread
        stats["SD"] = np.std(data)

        # Min/Max
        stats["Min"] = np.min(data)
        stats["Max"] = np.max(data)

        # Percentiles
        if percentiles != None:
            percentile = np.percentile(data, percentiles)
            for p in range(len(percentiles)):
                stats[str(percentiles[p]) + "th"] = percentile[p]

        return stats
    
#%%#########################################################################%%#
class Write:
    ## Write rms distance files ###############################################
    @staticmethod
    def rms_distance(file_name, bone_name, data, percentiles, data_path, i):
        # Open text file
        rms_path = os.path.join(data_path, file_name + ".txt")
        rms_file = open(rms_path, 'a')
        
        # Write headers
        if i==0:
            rms_file.write("Case,Mean,SD,Min,Max")
            for p in percentiles:
                rms_file.write(",Perc " + str(p))
            rms_file.write("\n")
        
        # Write data
        rms_file.write(bone_name)
        for d in range(len(data)):
            rms_file.write(",%6.6f" % (data[d]))
        rms_file.write("\n")
            
        # Close text file
        rms_file.close()
            
        return
        
#%%#########################################################################%%#
        