"""
Code to find the landmarks of the pelvis, femur, and tibfib and align them in the same global CS
1.  Rough alignment of meshes from manually identified bone landmarks to the coordinate system as follows:
    For the pelvis:
        o: placed at the midpoint of the left ASIS and right ASIS of the pelvis
        z: a line connecting the left ASIS and right ASIS of the pelvis
        n1: normal to the plane connecting the midpoint of the left and right PSIS, the left ASIS, and the right ASIS. This is an intermittent axis to ensure orthogonality.
        x: a line perpendicular to z and n1
        y: a line perpendicular to x and z

    For the femur:
        o: placed at the midpoint of the femoral epicondyles
        y: a line connecting the midpoint of the femoral epicondyles and the femoral head centre, pointing proximally.
        n1: normal to the plane connecting the femoral head centre, the medial epicondyle, and the lateral epicondyle. This is an intermittent axis to ensure orthogonality.
        z: a line perpendicular to y and n1
        x: a line perpendicular to y and z

    For the tibia-fibula:
        o: placed at the midpoint of the lateral and medial malleoli
        z: a line connecting the lateral and medial tibial condyles
        n1: normal to the plane connecting the midpoint of the malleoli, the medial condyle, and the lateral condyle. This is an intermittent axis to ensure orthogonality.
        y: a line perpendicular to z and n1
        x: a line perpendicular to y and z

2.  Determination of bony landmarks
3.  Re-alignment of all bones where:
        - origin at mid-ASIS
        - coincidence of acetabular centre and femoral head centre
        - offset of 5mm between the most distal point of the femur and the most proximal point of the tibia
author: Laura Carman
"""
import os
from find_landmarks_and_align_functions_hemi import PelvisPLY, FemurPLY, TibfibPLY
import write_results_functions as wrf

############## TO CHANGE ###############################
# Directory parameters
os.chdir("C:/Users/lcar475/Documents/Combined_shape_model/landmarks_all_cases")
root_dir = os.getcwd()
fit_fold = os.path.join(root_dir, 'fitted_meshes_CT_scan_CS_hemi') #where are the meshes, should be split into pelvis, femur and tibfib meshes
align_fold = os.path.join(root_dir, 'aligned_meshes_CS_3_hemi') #where will the aligned meshes go
back_transform_landmarks = True
if not os.path.exists(align_fold):
    os.makedirs(align_fold)
data_path = os.path.join(root_dir, 'data_hemi') #where will the data go (landmarks, measurements, and transformation matrix files)
if not os.path.exists(data_path):
    os.makedirs(data_path)
name = 'hemi' #name to append to the measurements file

#manually selected landmarks for initial alignment
pel_landmarks = [548, 3250] # ASIS, PSIS
fem_landmarks = [3626, 214] # lat_epi, med_epi (left side)
tib_landmarks = [421, 2829, 3403, 4201] # lat_mall, med_mall, lat_con, med_con (left side)

list_fold = 'C:/Users/lcar475/Documents/Combined_shape_model/landmark_lists_2' #where the selected node lists are

suff_num = 16 # 18 for fitted meshes, 11 for mean mesh # number to remove from pelvis name to get case name
suffix = '_rbfreg.ply' # '_rbfreg.ply' for fitted meshes, '.ply' for mean mesh # what to add on after the case name
########################################################

landmarks_path = os.path.join(data_path, "landmarks")
if not os.path.exists(landmarks_path):
    os.makedirs(landmarks_path)
measurements_path = os.path.join(data_path, "measurements")
if not os.path.exists(measurements_path):
    os.makedirs(measurements_path)

tm_fold = os.path.join(data_path, "tms")
if not os.path.exists(tm_fold):
    os.makedirs(tm_fold)

file_path_lengths = os.path.join(measurements_path, "measurements.txt")
if os.path.exists(file_path_lengths):
    os.remove(file_path_lengths)
file_fit_lengths = open(file_path_lengths, "w+")

pel_case = [f for f in os.listdir(fit_fold) if "Left_pelvis" in f and '.ply' in f and '._' not in f]
all_cases = [f for f in os.listdir(fit_fold) if f.endswith('.ply') and '._' not in f]
turn = 1
for case in sorted(pel_case[:]):
    case_num = case[:-suff_num]
    if case_num + '_Right_pelvis' + suffix in all_cases:
        x=0
    else:
        continue
    if case_num + '_Left_femur' + suffix in all_cases:
        x=1
    else:
        continue
    if case_num + '_Left_tibfib' + suffix in all_cases:
        x=2
    else:
        continue
    if case_num + '_Right_femur' + suffix in all_cases:
        x = 3
    else:
        continue
    if case_num + '_Right_tibfib' + suffix in all_cases:
        x = 4
    else:
        continue
    print("Starting: " + case)
    # set up text files
    file_path_angles = os.path.join(measurements_path, "angles.txt")
    if os.path.exists(file_path_angles):
        os.remove(file_path_angles)

    # align and take measurements of the left and right hemi of the pelvis
    Pelvis = PelvisPLY(case)
    Pelvis.load_vertices(os.path.join(fit_fold, case))
    Pelvis.load_landmark_lists(list_fold)
    #load right hemi
    case_pelvis_r = case_num + '_Right_pelvis' + suffix
    Right_pelvis = PelvisPLY(case_pelvis_r)
    Right_pelvis.load_vertices(os.path.join(fit_fold, case_pelvis_r))
    Right_pelvis.load_landmark_lists(list_fold)
    tm_pelvis = Pelvis.find_lm_and_align(Right_pelvis,pel_landmarks,align_fold)

    # align and take measurements femurs
    case_femur_l = case_num + "_Left_femur" + suffix
    Left_femur = FemurPLY(case_femur_l, side = 'left')
    Left_femur.load_vertices(os.path.join(fit_fold, case_femur_l))
    Left_femur.load_landmark_lists(list_fold)
    tm_femur_l = Left_femur.find_lm_and_align(fem_landmarks,Pelvis, align_fold)

    case_femur_r = case_num + "_Right_femur" + suffix
    Right_femur = FemurPLY(case_femur_r, side='right')
    Right_femur.load_vertices(os.path.join(fit_fold, case_femur_r))
    Right_femur.load_landmark_lists(list_fold)
    tm_femur_r = Right_femur.find_lm_and_align(fem_landmarks,Right_pelvis, align_fold)

    # align and take measurements of the tibfib
    case_tibfib = case_num + "_Left_tibfib" + suffix
    Left_tibfib = TibfibPLY(case_tibfib, 'left')
    Left_tibfib.load_vertices(os.path.join(fit_fold, case_tibfib))
    Left_tibfib.load_landmark_lists(list_fold)
    tm_tibfib_l = Left_tibfib.find_lm_and_align(tib_landmarks,Left_femur,align_fold)

    case_tibfib = case_num + "_Right_tibfib" + suffix
    Right_tibfib = TibfibPLY(case_tibfib, 'right')
    Right_tibfib.load_vertices(os.path.join(fit_fold, case_tibfib))
    Right_tibfib.load_landmark_lists(list_fold)
    tm_tibfib_r = Right_tibfib.find_lm_and_align(tib_landmarks,Right_femur,align_fold)

    if turn == 1:
        wrf.setup_measurements_file_header(file_fit_lengths,Pelvis,Left_femur,Left_tibfib)

    if back_transform_landmarks:
        Pelvis.back_transform_landmarks(tm_pelvis)
        Right_pelvis.back_transform_landmarks(tm_pelvis)
        Left_femur.back_transform_landmarks(tm_femur_l)
        Left_tibfib.back_transform_landmarks(tm_tibfib_l)
        Right_femur.back_transform_landmarks(tm_femur_r)
        Right_tibfib.back_transform_landmarks(tm_tibfib_r)
        original_landmarks_path = os.path.join(data_path, 'original_landmarks')
        if not os.path.exists(original_landmarks_path):
            os.makedirs(original_landmarks_path)
        file_path_original_landmarks = os.path.join(original_landmarks_path, case_num + "_original_landmarks.txt")
        file_path_original_mocap_landmarks = os.path.join(original_landmarks_path,
                                                          case_num + "_original_landmarks_mocap.txt")
        if os.path.exists(file_path_original_landmarks):
            os.remove(file_path_original_landmarks)
        if os.path.exists(file_path_original_mocap_landmarks):
            os.remove(file_path_original_mocap_landmarks)

        wrf.write_original_landmarks_file(file_path_original_landmarks,file_path_original_mocap_landmarks,Pelvis,Right_pelvis,Left_femur,Right_femur,Left_tibfib,Right_tibfib)

    # save landmarks and measurements
    file_path_landmarks = os.path.join(landmarks_path, case_num + "_landmarks.txt")
    file_path_mocap_landmarks = os.path.join(landmarks_path, case_num + "_landmarks_mocap.txt")
    if os.path.exists(file_path_landmarks):
        os.remove(file_path_landmarks)
    if os.path.exists(file_path_mocap_landmarks):
        os.remove(file_path_mocap_landmarks)
    wrf.write_landmarks_file(case_num,file_path_landmarks,file_path_mocap_landmarks,file_fit_lengths,Pelvis,Right_pelvis,Left_femur,Right_femur,Left_tibfib,Right_tibfib)
    turn = 2