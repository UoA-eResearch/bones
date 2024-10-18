"""
Code to find the landmarks of the pelvis, femur, and tibfib and align them in the same global CS
1.  Rough alignment of meshes from manually identified bone landmarks to the coordinate system defined by ISB (Wu et. al,
    2005) as follows:
    For the pelvis: origin = mid-point of the ASIS'
        x = line parallel to the line connecting the mid-point of the PSIS' and
            the ASIS line which is perpendicular to the z line (posterior -> anterior)
        y = line perpendicular to both x and z axis and pointing upwards (distal -> proximal)
        z = line parallel to a line connecting the left and right ASIS going
            from the origin and pointing to the right (of the bone) (lateral -> medial)
    For the femur:  origin = the centre of the femur head
        x = line perpendicular to the y and z axis (posterior -> anterior)
        y = from the mid-point between the two epicondyles to the femoral head centre (distal -> proximal)
        z = line parallel to the line connecting the lateral and medial epicondyle (lateral -> medial)
    For the tibia-fibula:   origin = the inter-malleolar point
        x = line perpendicular to the torsional plane (plane containing the malleoli and the intercondyle point) and
            pointing anteriorly (posterior -> anterior)
        y = line perpendicular to both the x and z axis and pointing superiorly (distal -> proximal)
        z = line parallel to the line connecting the medial and lateral malleolus (lateral -> medial)
2.  Determination of bony landmarks
3.  Re-alignment of bones where:
        - origin at mid-ASIS
        - coincidence of acetabular centre and femoral head centre
        - offset of 5mm between the most distal point of the femur and the most proximal point of the tibia
author: Laura Carman 22/06/22
"""
import os
from find_landmarks_and_align_functions import PelvisPLY, FemurPLY, TibfibPLY


############## TO CHANGE ###############################
# Directory parameters
#os.chdir('F:\ASM')
#os.chdir("/Volumes/Laura2TB/ASM")
os.chdir("/Users/lauracarman/Desktop/single_bone_ssm/scaled_meshes_paed")
#os.chdir('/Users/lauracarman/Desktop/Lower_limb_models/PLB-03') #set working directory
root_dir = os.getcwd()
fit_fold = os.path.join(root_dir, 'cc_aligned') #where are the meshes, should be split into pelvis, femur and tibfib meshes
align_fold = os.path.join(root_dir, 'aligned_meshes_CS_3') #where will the aligned meshes go
back_transform_landmarks = False
if not os.path.exists(align_fold):
    os.makedirs(align_fold)
data_path = os.path.join(root_dir, 'data_cc_aligned_meshes') #where will the data go (landmarks, measurements, and transformation matrix files)
if not os.path.exists(data_path):
    os.makedirs(data_path)
name = 'Left_side_pred' #name to append to the measurements file
########################################################

landmarks_path = os.path.join(data_path, "landmarks")
if not os.path.exists(landmarks_path):
    os.makedirs(landmarks_path)
measurements_path = os.path.join(data_path, "measurements")
if not os.path.exists(measurements_path):
    os.makedirs(measurements_path)
#list_fold = os.path.join(root_dir, "landmark_lists")
list_fold = '/Users/lauracarman/Desktop/single_bone_ssm/landmark_lists_2'
tm_fold = os.path.join(data_path, "tms")
if not os.path.exists(tm_fold):
    os.makedirs(tm_fold)

# manually identified landmarks for ASM_Sacrum
#pel_landmarks = [55159, 17482, 46474, 27577] # LASIS, RASIS, LPSIS, RPSIS
#fem_landmarks = [20902, 3694] # lat_epi, med_epi
#tib_landmarks = [12304, 17820, 39741, 23671] # lat_mall, med_mall, lat_con, med_con

#pel_landmarks = [55165, 9610, 46474, 27577] # LASIS, RASIS, LPSIS, RPSIS
#fem_landmarks = [22012, 4057] # lat_epi, med_epi
#tib_landmarks = [1605, 17820, 39744, 23748] # lat_mall, med_mall, lat_con, med_con

pel_landmarks = [5159, 1723, 7861, 532] # LASIS, RASIS, LPSIS, RPSIS
fem_landmarks = [3606, 126] # lat_epi, med_epi
tib_landmarks = [421, 2829, 4179, 4201] # lat_mall, med_mall, lat_con, med_con

#pel_landmarks = [54970, 17479, 46513, 2434] # LASIS, RASIS, LPSIS, RPSIS
#fem_landmarks = [21454, 3667] # lat_epi, med_epi
#tib_landmarks = [11296, 16181, 18908, 26659] # lat_mall, med_mall, lat_con, med_con

file_path_lengths = os.path.join(measurements_path, "measurements.txt")
if os.path.exists(file_path_lengths):
    os.remove(file_path_lengths)
file_fit_lengths = open(file_path_lengths, "w+")
file_fit_lengths.write("Case_name,ASIS_width,LHJ_diameter,PSIS_width,RHJ_diameter,Pelvis_depth,AA,BA,NSA,mLDFA,fem_head_diameter,epicon_width,fem_length,TT,mMPTA,condylar_width,malleolar_width,tibial_length\n")
cases = [f for f in os.listdir(fit_fold) if "Pelvis" in f and '.ply' in f and '._' not in f]
all_cases = [f for f in os.listdir(fit_fold) if f.endswith('.ply') and '._' not in f]
for case in sorted(cases[:]):
    case_num = case[:-34] #for mean mesh
    #case_num = case[:-20]
    #case_num = case[:-26] #for aligned meshes
    #case_num = case[:-17] # for fitted meshes
    #suff = '_rbfreg_rigidreg.ply' # for aligned meshes
    #suff = '_rbfreg.ply' # for fitted meshes
    if case_num + 'Left_femur_generically_scaled_mean.ply' in all_cases: #"Left_femur" + suff in all_cases: #'_Left_femur.ply' in all_cases:
        x=1
    else:
        continue
    if case_num + 'Left_tibfib_generically_scaled_mean.ply' in all_cases: # "Left_tibfib" + suff in all_cases: #'_Left_tibfib.ply' in all_cases:
        x=2
    else:
        continue
    print("Starting: " + case)
    # set up text files
    file_path_angles = os.path.join(measurements_path, "angles.txt")
    if os.path.exists(file_path_angles):
        os.remove(file_path_angles)

    # align and take measurements of the pelvis
    Pelvis = PelvisPLY(case)
    Pelvis.load_vertices(os.path.join(fit_fold,case))
    Pelvis.load_landmark_lists(list_fold)
    #tm_1 = Pelvis.init_align(pel_landmarks)
    #tm_1 = Pelvis.init_align_ISB(pel_landmarks)
    Pelvis.find_landmarks()
    #tm_2 = Pelvis.align()
    #tm_pelvis = tm_2.dot(tm_1)
    #tm_3 = Pelvis.align_ISB()
    #tm_pelvis = tm_3.dot(tm_1)
    Pelvis.find_landmarks()
    Pelvis.length_measurements()
    Pelvis.update_acs()
    #Pelvis.save_ply(align_fold)

    # align and take measurements of the femur
    #case_femur = case[:-10]+"Femur.ply" #for mean mesh
    #case_femur = case_num + "Left_femur"+suff
    #case_femur = case_num + "_Left_femur.ply"
    case_femur = case_num + "Left_femur_generically_scaled_mean.ply"
    Femur = FemurPLY(case_femur)
    Femur.load_vertices(os.path.join(fit_fold, case_femur))
    Femur.load_landmark_lists(list_fold)
    #tm_1 = Femur.init_align(fem_landmarks)
    #tm_1 = Femur.init_align_ISB(fem_landmarks)
    Femur.find_landmarks()
    #tm_2 = Femur.align(Pelvis)
    #tm_femur = tm_2.dot(tm_1)
    #tm_2 = Femur.align_ISB(Pelvis)
    #tm_femur = tm_2.dot(tm_1)
    Femur.find_landmarks()
    Femur.measurements()
    Femur.update_acs()
    #Femur.save_ply(align_fold)

    # align and take measurements of the tibfib
    #case_tibfib = case[:-10] + "Tibfib.ply" #for mean mesh
    #case_tibfib = case_num + "Left_tibfib"+suff
    #case_tibfib = case_num + "_Left_tibfib.ply"
    case_tibfib = case_num + "Left_tibfib_generically_scaled_mean.ply"
    Tibfib = TibfibPLY(case_tibfib)
    Tibfib.load_vertices(os.path.join(fit_fold, case_tibfib))
    Tibfib.load_landmark_lists(list_fold)

    #tm_1 = Tibfib.init_align(tib_landmarks)
    #tm_1 = Tibfib.init_align_ISB(tib_landmarks)
    Tibfib.find_landmarks()
    #tm_2 = Tibfib.align(Femur.landmarks)
    #tm_tibfib = tm_2.dot(tm_1)
    #tm_2 = Tibfib.align_ISB(Femur)
    #tm_tibfib = tm_2.dot(tm_1)
    Tibfib.find_landmarks()
    Tibfib.update_acs()
    #tm_3 = Tibfib.reset_knee_gap(Femur)
    #tm_tifib = tm_3.dot(tm_tibfib)
    Tibfib.find_landmarks()
    Tibfib.update_acs()
    #tm_4 = Tibfib.reset_knee_gap_2(Femur)
    #tm_tibfib = tm_4.dot(tm_3).dot(tm_2).dot(tm_1)
    Tibfib.find_landmarks()
    Tibfib.measurements()
    #Tibfib.save_ply(align_fold)


    # save landmarks and measurements
    file_path_landmarks = os.path.join(landmarks_path, case_num + "landmarks.txt")
    file_path_mocap_landmarks = os.path.join(landmarks_path, case_num + "landmarks_mocap.txt")
    if os.path.exists(file_path_landmarks):
        os.remove(file_path_landmarks)
    if os.path.exists(file_path_mocap_landmarks):
        os.remove(file_path_mocap_landmarks)
    file_fit_landmarks = open(file_path_landmarks, "w+")
    file_fit_landmarks.write("Landmark,x,y,z,ID\n")
    file_fit_lengths.write(str(case_num))
    file_fit_mocap_landmarks = open(file_path_mocap_landmarks, "w+")
    ASIS = Pelvis.landmarks["LASIS"]["coords"]
    RASIS = Pelvis.landmarks['RASIS']['coords']
    PSIS = Pelvis.landmarks['LPSIS']['coords']
    RPSIS = Pelvis.landmarks['RPSIS']['coords']
    file_fit_mocap_landmarks.write(
        str('left-ASIS') + ' ' + str(ASIS[0]) + ' ' + str(ASIS[1]) + ' ' + str(ASIS[2]) + '\n')
    file_fit_mocap_landmarks.write(
        str('right-ASIS') + ' ' + str(RASIS[0]) + ' ' + str(RASIS[1]) + ' ' + str(RASIS[2]) + '\n')
    file_fit_mocap_landmarks.write(
        str('left-PSIS') + ' ' + str(PSIS[0]) + ' ' + str(PSIS[1]) + ' ' + str(PSIS[2]) + '\n')
    file_fit_mocap_landmarks.write(
        str('right-PSIS') + ' ' + str(RPSIS[0]) + ' ' + str(RPSIS[1]) + ' ' + str(RPSIS[2]) + '\n')
    lat = Femur.landmarks["lat_epicon"]["coords"]
    med = Femur.landmarks["med_epicon"]["coords"]
    file_fit_mocap_landmarks.write(
        str('left-LEC') + ' ' + str(lat[0]) + ' ' + str(lat[1]) + ' ' + str(lat[2]) + '\n')
    file_fit_mocap_landmarks.write(
        str('left-MEC') + ' ' + str(med[0]) + ' ' + str(med[1]) + ' ' + str(med[2]) + '\n')
    lattf = Tibfib.landmarks['lateral_malleolus']['coords']
    medtf = Tibfib.landmarks['medial_malleolus']['coords']
    file_fit_mocap_landmarks.write(str('left-malleolus_lat') + ' ' + str(lattf[0]) + ' ' + str(lattf[1]) + ' ' + str(lattf[2]) + '\n')
    file_fit_mocap_landmarks.write(str('left-malleolus_med') + ' ' + str(medtf[0]) + ' ' + str(medtf[1]) + ' ' + str(medtf[2]) + '\n')
    #for putting the node ID in the local coordinate system add 55758 to the femur node number and 86442 to the tibfib node number

    for ld in sorted(Pelvis.landmarks):
        if ld == "lengths":
            for ll in sorted(Pelvis.landmarks[ld]):
                file_fit_lengths.write(',' + str(Pelvis.landmarks[ld][ll]))
        elif ld == "angles":
            for ll in sorted(Pelvis.landmarks[ld]):
                file_fit_lengths.write(',' + str(Pelvis.landmarks[ld][ll]))
        else:
            coords = Pelvis.landmarks[ld]['coords']
            try:
                ID = Pelvis.landmarks[ld]['ID']
                file_fit_landmarks.write(
                    str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_landmarks.write(str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')

    for ld in sorted(Femur.landmarks):

        if ld == "lengths":
            for ll in sorted(Femur.landmarks[ld]):
                file_fit_lengths.write(',' + str(Femur.landmarks[ld][ll]))
        elif ld == "angles":
            for ll in sorted(Femur.landmarks[ld]):
                file_fit_lengths.write(',' + str(Femur.landmarks[ld][ll]))
        else:
            coords = Femur.landmarks[ld]['coords']
            try:
                ID = Femur.landmarks[ld]['ID']
                file_fit_landmarks.write(str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_landmarks.write(
                    str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')
    for ld in sorted(Tibfib.landmarks):
        if ld == "lengths":
            for ll in sorted(Tibfib.landmarks[ld]):
                file_fit_lengths.write(',' + str(Tibfib.landmarks[ld][ll]))
            file_fit_lengths.write('\n')
        elif ld == "angles":
            for ll in sorted(Tibfib.landmarks[ld]):
                file_fit_lengths.write(',' + str(Tibfib.landmarks[ld][ll]))
        else:
            coords = Tibfib.landmarks[ld]['coords']
            try:
                ID = Tibfib.landmarks[ld]['ID']
                file_fit_landmarks.write(
                    str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_landmarks.write(str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')
