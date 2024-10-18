'''
functions to write the results of bone measurement calculations
'''

def setup_measurements_file_header(file_fit_lengths,pel,fem,tf):
    file_fit_lengths.write('Case_name')
    for length in sorted(pel.landmarks["lengths"]):
        file_fit_lengths.write(',' + str(length))
    for angle in sorted(fem.landmarks['angles']):
        file_fit_lengths.write(',' + 'left_' + str(angle))
    for length in sorted(fem.landmarks['lengths']):
        file_fit_lengths.write(',' + 'left_' + str(length))
    for angle in sorted(tf.landmarks['angles']):
        file_fit_lengths.write(',' + 'left_' + str(angle))
    for length in sorted(tf.landmarks['lengths']):
        file_fit_lengths.write(',' + 'left_' + str(length))
    for angle in sorted(fem.landmarks['angles']):
        file_fit_lengths.write(',' + 'left_' + str(angle))
    for length in sorted(fem.landmarks['lengths']):
        file_fit_lengths.write(',' + 'left_' + str(length))
    for angle in sorted(tf.landmarks['angles']):
        file_fit_lengths.write(',' + 'left_' + str(angle))
    for length in sorted(tf.landmarks['lengths']):
        file_fit_lengths.write(',' + 'left_' + str(length))

    file_fit_lengths.write('\n')
    return

def write_landmarks_file(case_num,file_path_landmarks,file_path_mocap_landmarks,file_fit_lengths,Pelvis,Right_pelvis,Left_femur,Right_femur,Left_tibfib,Right_tibfib):
    file_fit_landmarks = open(file_path_landmarks, "w+")
    file_fit_landmarks.write("Landmark,x,y,z,ID\n")
    file_fit_lengths.write(str(case_num))
    file_fit_mocap_landmarks = open(file_path_mocap_landmarks, "w+")
    ASIS = Pelvis.landmarks["LASIS"]["coords"]
    RASIS = Right_pelvis.landmarks['RASIS']['coords']
    PSIS = Pelvis.landmarks['LPSIS']['coords']
    RPSIS = Right_pelvis.landmarks['RPSIS']['coords']
    file_fit_mocap_landmarks.write(
        str('left-ASIS') + ' ' + str(ASIS[0]) + ' ' + str(ASIS[1]) + ' ' + str(ASIS[2]) + '\n')
    file_fit_mocap_landmarks.write(
        str('right-ASIS') + ' ' + str(RASIS[0]) + ' ' + str(RASIS[1]) + ' ' + str(RASIS[2]) + '\n')
    file_fit_mocap_landmarks.write(
        str('left-PSIS') + ' ' + str(PSIS[0]) + ' ' + str(PSIS[1]) + ' ' + str(PSIS[2]) + '\n')
    file_fit_mocap_landmarks.write(
        str('right-PSIS') + ' ' + str(RPSIS[0]) + ' ' + str(RPSIS[1]) + ' ' + str(RPSIS[2]) + '\n')
    lat = Left_femur.landmarks["lat_epicon"]["coords"]
    med = Left_femur.landmarks["med_epicon"]["coords"]
    file_fit_mocap_landmarks.write(
        str('left-LEC') + ' ' + str(lat[0]) + ' ' + str(lat[1]) + ' ' + str(lat[2]) + '\n')
    file_fit_mocap_landmarks.write(
        str('left-MEC') + ' ' + str(med[0]) + ' ' + str(med[1]) + ' ' + str(med[2]) + '\n')
    lattf = Left_tibfib.landmarks['lateral_malleolus']['coords']
    medtf = Left_tibfib.landmarks['medial_malleolus']['coords']
    file_fit_mocap_landmarks.write(
        str('left-malleolus_lat') + ' ' + str(lattf[0]) + ' ' + str(lattf[1]) + ' ' + str(lattf[2]) + '\n')
    file_fit_mocap_landmarks.write(
        str('left-malleolus_med') + ' ' + str(medtf[0]) + ' ' + str(medtf[1]) + ' ' + str(medtf[2]) + '\n')
    lat = Right_femur.landmarks["lat_epicon"]["coords"]
    med = Right_femur.landmarks["med_epicon"]["coords"]
    file_fit_mocap_landmarks.write(
        str('right-LEC') + ' ' + str(lat[0]) + ' ' + str(lat[1]) + ' ' + str(lat[2]) + '\n')
    file_fit_mocap_landmarks.write(
        str('right-MEC') + ' ' + str(med[0]) + ' ' + str(med[1]) + ' ' + str(med[2]) + '\n')
    lattf = Right_tibfib.landmarks['lateral_malleolus']['coords']
    medtf = Right_tibfib.landmarks['medial_malleolus']['coords']
    file_fit_mocap_landmarks.write(
        str('right-malleolus_lat') + ' ' + str(lattf[0]) + ' ' + str(lattf[1]) + ' ' + str(lattf[2]) + '\n')
    file_fit_mocap_landmarks.write(
        str('right-malleolus_med') + ' ' + str(medtf[0]) + ' ' + str(medtf[1]) + ' ' + str(medtf[2]) + '\n')

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
                file_fit_landmarks.write(
                    str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')

    for ld in sorted(Right_pelvis.landmarks):
        if ld == "lengths":
            for ll in sorted(Pelvis.landmarks[ld]):
                file_fit_lengths.write(',' + str(Pelvis.landmarks[ld][ll]))
        elif ld == "angles":
            for ll in sorted(Pelvis.landmarks[ld]):
                file_fit_lengths.write(',' + str(Pelvis.landmarks[ld][ll]))
        else:
            coords = Right_pelvis.landmarks[ld]['coords']
            try:
                ID = Right_pelvis.landmarks[ld]['ID']
                file_fit_landmarks.write(
                    str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_landmarks.write(
                    str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')

    for ld in sorted(Left_femur.landmarks):

        if ld == "lengths":
            for ll in sorted(Left_femur.landmarks[ld]):
                file_fit_lengths.write(',' + str(Left_femur.landmarks[ld][ll]))
        elif ld == "angles":
            for ll in sorted(Left_femur.landmarks[ld]):
                file_fit_lengths.write(',' + str(Left_femur.landmarks[ld][ll]))
        else:
            coords = Left_femur.landmarks[ld]['coords']
            try:
                ID = Left_femur.landmarks[ld]['ID']
                file_fit_landmarks.write(
                    str(ld) + '_left' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_landmarks.write(
                    str(ld) + '_left' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')
    for ld in sorted(Left_tibfib.landmarks):
        if ld == "lengths":
            for ll in sorted(Left_tibfib.landmarks[ld]):
                file_fit_lengths.write(',' + str(Left_tibfib.landmarks[ld][ll]))
            file_fit_lengths.write('\n')
        elif ld == "angles":
            for ll in sorted(Left_tibfib.landmarks[ld]):
                file_fit_lengths.write(',' + str(Left_tibfib.landmarks[ld][ll]))
        else:
            coords = Left_tibfib.landmarks[ld]['coords']
            try:
                ID = Left_tibfib.landmarks[ld]['ID']
                file_fit_landmarks.write(
                    str(ld) + '_left' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_landmarks.write(
                    str(ld) + '_left' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')

    for ld in sorted(Right_femur.landmarks):

        if ld == "lengths":
            for ll in sorted(Right_femur.landmarks[ld]):
                file_fit_lengths.write(',' + str(Right_femur.landmarks[ld][ll]))
        elif ld == "angles":
            for ll in sorted(Right_femur.landmarks[ld]):
                file_fit_lengths.write(',' + str(Right_femur.landmarks[ld][ll]))
        else:
            coords = Right_femur.landmarks[ld]['coords']
            try:
                ID = Right_femur.landmarks[ld]['ID']
                file_fit_landmarks.write(
                    str(ld) + '_right' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_landmarks.write(
                    str(ld) + '_right' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + '\n')
    for ld in sorted(Right_tibfib.landmarks):
        if ld == "lengths":
            for ll in sorted(Right_tibfib.landmarks[ld]):
                file_fit_lengths.write(',' + str(Right_tibfib.landmarks[ld][ll]))
            file_fit_lengths.write('\n')
        elif ld == "angles":
            for ll in sorted(Right_tibfib.landmarks[ld]):
                file_fit_lengths.write(',' + str(Right_tibfib.landmarks[ld][ll]))
        else:
            coords = Right_tibfib.landmarks[ld]['coords']
            try:
                ID = Right_tibfib.landmarks[ld]['ID']
                file_fit_landmarks.write(
                    str(ld) + '_right' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_landmarks.write(
                    str(ld) + '_right' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + '\n')

def write_landmarks_file_2(case_num,file_path_landmarks,file_path_mocap_landmarks,file_fit_lengths,Pelvis,femur,tibfib,side):
    file_fit_landmarks = open(file_path_landmarks, "w+")
    file_fit_landmarks.write("Landmark,x,y,z,ID\n")
    file_fit_lengths.write(str(case_num))
    file_fit_mocap_landmarks = open(file_path_mocap_landmarks, "w+")
    if side == 'Left':
        ASIS = Pelvis.landmarks["LASIS"]["coords"]
        PSIS = Pelvis.landmarks['LPSIS']['coords']
        file_fit_mocap_landmarks.write(
            str('left-ASIS') + ' ' + str(ASIS[0]) + ' ' + str(ASIS[1]) + ' ' + str(ASIS[2]) + '\n')
        file_fit_mocap_landmarks.write(
            str('left-PSIS') + ' ' + str(PSIS[0]) + ' ' + str(PSIS[1]) + ' ' + str(PSIS[2]) + '\n')
    else:
        ASIS = Pelvis.landmarks["RASIS"]["coords"]
        PSIS = Pelvis.landmarks['RPSIS']['coords']
        file_fit_mocap_landmarks.write(
            str('right-ASIS') + ' ' + str(ASIS[0]) + ' ' + str(ASIS[1]) + ' ' + str(ASIS[2]) + '\n')
        file_fit_mocap_landmarks.write(
            str('right-PSIS') + ' ' + str(PSIS[0]) + ' ' + str(PSIS[1]) + ' ' + str(PSIS[2]) + '\n')


    lat = femur.landmarks["lat_epicon"]["coords"]
    med = femur.landmarks["med_epicon"]["coords"]
    lattf = tibfib.landmarks['lateral_malleolus']['coords']
    medtf = tibfib.landmarks['medial_malleolus']['coords']

    if side == "Left":
        file_fit_mocap_landmarks.write(
            str('left-LEC') + ' ' + str(lat[0]) + ' ' + str(lat[1]) + ' ' + str(lat[2]) + '\n')
        file_fit_mocap_landmarks.write(
            str('left-MEC') + ' ' + str(med[0]) + ' ' + str(med[1]) + ' ' + str(med[2]) + '\n')
        file_fit_mocap_landmarks.write(
            str('left-malleolus_lat') + ' ' + str(lattf[0]) + ' ' + str(lattf[1]) + ' ' + str(lattf[2]) + '\n')
        file_fit_mocap_landmarks.write(
            str('left-malleolus_med') + ' ' + str(medtf[0]) + ' ' + str(medtf[1]) + ' ' + str(medtf[2]) + '\n')
    else:
        file_fit_mocap_landmarks.write(
            str('right-LEC') + ' ' + str(lat[0]) + ' ' + str(lat[1]) + ' ' + str(lat[2]) + '\n')
        file_fit_mocap_landmarks.write(
            str('right-MEC') + ' ' + str(med[0]) + ' ' + str(med[1]) + ' ' + str(med[2]) + '\n')
        file_fit_mocap_landmarks.write(
            str('right-malleolus_lat') + ' ' + str(lattf[0]) + ' ' + str(lattf[1]) + ' ' + str(lattf[2]) + '\n')
        file_fit_mocap_landmarks.write(
            str('right-malleolus_med') + ' ' + str(medtf[0]) + ' ' + str(medtf[1]) + ' ' + str(medtf[2]) + '\n')

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
                file_fit_landmarks.write(
                    str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')

    for ld in sorted(femur.landmarks):

        if ld == "lengths":
            for ll in sorted(femur.landmarks[ld]):
                file_fit_lengths.write(',' + str(femur.landmarks[ld][ll]))
        elif ld == "angles":
            for ll in sorted(femur.landmarks[ld]):
                file_fit_lengths.write(',' + str(femur.landmarks[ld][ll]))
        else:
            coords = femur.landmarks[ld]['coords']
            try:
                ID = femur.landmarks[ld]['ID']
                file_fit_landmarks.write(
                    str(ld) + '_' + side + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_landmarks.write(
                    str(ld) + '_' + side + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')
    for ld in sorted(tibfib.landmarks):
        if ld == "lengths":
            for ll in sorted(tibfib.landmarks[ld]):
                file_fit_lengths.write(',' + str(tibfib.landmarks[ld][ll]))
            file_fit_lengths.write('\n')
        elif ld == "angles":
            for ll in sorted(tibfib.landmarks[ld]):
                file_fit_lengths.write(',' + str(tibfib.landmarks[ld][ll]))
        else:
            coords = tibfib.landmarks[ld]['coords']
            try:
                ID = tibfib.landmarks[ld]['ID']
                file_fit_landmarks.write(
                    str(ld) + '_' + side + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_landmarks.write(
                    str(ld) + '_' + side + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')

def write_original_landmarks_file(file_path_original_landmarks,file_path_original_mocap_landmarks,Pelvis,Right_pelvis,Left_femur,Right_femur,Left_tibfib,Right_tibfib):
    file_fit_original_landmarks = open(file_path_original_landmarks, "w+")
    file_fit_original_landmarks.write("Landmark,x,y,z,ID\n")
    file_fit_original_mocap_landmarks = open(file_path_original_mocap_landmarks, "w+")
    ASIS = Pelvis.original_landmarks["LASIS"]["coords"]
    RASIS = Right_pelvis.original_landmarks['RASIS']['coords']
    PSIS = Pelvis.original_landmarks['LPSIS']['coords']
    RPSIS = Right_pelvis.original_landmarks['RPSIS']['coords']
    file_fit_original_mocap_landmarks.write(
        str('left-ASIS') + ' ' + str(ASIS[0]) + ' ' + str(ASIS[1]) + ' ' + str(ASIS[2]) + '\n')
    file_fit_original_mocap_landmarks.write(
        str('right-ASIS') + ' ' + str(RASIS[0]) + ' ' + str(RASIS[1]) + ' ' + str(RASIS[2]) + '\n')
    file_fit_original_mocap_landmarks.write(
        str('left-PSIS') + ' ' + str(PSIS[0]) + ' ' + str(PSIS[1]) + ' ' + str(PSIS[2]) + '\n')
    file_fit_original_mocap_landmarks.write(
        str('right-PSIS') + ' ' + str(RPSIS[0]) + ' ' + str(RPSIS[1]) + ' ' + str(RPSIS[2]) + '\n')
    lat = Left_femur.original_landmarks["lat_epicon"]["coords"]
    med = Left_femur.original_landmarks["med_epicon"]["coords"]
    file_fit_original_mocap_landmarks.write(
        str('left-LEC') + ' ' + str(lat[0]) + ' ' + str(lat[1]) + ' ' + str(lat[2]) + '\n')
    file_fit_original_mocap_landmarks.write(
        str('left-MEC') + ' ' + str(med[0]) + ' ' + str(med[1]) + ' ' + str(med[2]) + '\n')
    lattf = Left_tibfib.original_landmarks['lateral_malleolus']['coords']
    medtf = Left_tibfib.original_landmarks['medial_malleolus']['coords']
    file_fit_original_mocap_landmarks.write(
        str('left-malleolus_lat') + ' ' + str(lattf[0]) + ' ' + str(lattf[1]) + ' ' + str(lattf[2]) + '\n')
    file_fit_original_mocap_landmarks.write(
        str('left-malleolus_med') + ' ' + str(medtf[0]) + ' ' + str(medtf[1]) + ' ' + str(medtf[2]) + '\n')
    lat = Right_femur.original_landmarks["lat_epicon"]["coords"]
    med = Right_femur.original_landmarks["med_epicon"]["coords"]
    file_fit_original_mocap_landmarks.write(
        str('right-LEC') + ' ' + str(lat[0]) + ' ' + str(lat[1]) + ' ' + str(lat[2]) + '\n')
    file_fit_original_mocap_landmarks.write(
        str('right-MEC') + ' ' + str(med[0]) + ' ' + str(med[1]) + ' ' + str(med[2]) + '\n')
    lattf = Right_tibfib.original_landmarks['lateral_malleolus']['coords']
    medtf = Right_tibfib.original_landmarks['medial_malleolus']['coords']
    file_fit_original_mocap_landmarks.write(
        str('right-malleolus_lat') + ' ' + str(lattf[0]) + ' ' + str(lattf[1]) + ' ' + str(lattf[2]) + '\n')
    file_fit_original_mocap_landmarks.write(
        str('right-malleolus_med') + ' ' + str(medtf[0]) + ' ' + str(medtf[1]) + ' ' + str(medtf[2]) + '\n')

    # for putting the node ID in the local coordinate system add 55758 to the femur node number and 86442 to the tibfib node number

    for ld in sorted(Pelvis.original_landmarks):
        if ld == "lengths":
            continue
        elif ld == "angles":
            continue
        else:
            coords = Pelvis.original_landmarks[ld]['coords']
            try:
                ID = Pelvis.original_landmarks[ld]['ID']
                file_fit_original_landmarks.write(
                    str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + ',' + str(
                        ID) + '\n')
            except KeyError:
                file_fit_original_landmarks.write(
                    str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')

    for ld in sorted(Right_pelvis.landmarks):
        if ld == "lengths":
            continue
        elif ld == "angles":
            continue
        else:
            coords = Right_pelvis.landmarks[ld]['coords']
            try:
                ID = Right_pelvis.landmarks[ld]['ID']
                file_fit_original_landmarks.write(
                    str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_original_landmarks.write(
                    str(ld) + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(coords[2]) + '\n')

    for ld in sorted(Left_femur.original_landmarks):
        if ld == "lengths":
            continue
        elif ld == "angles":
            continue
        else:
            coords = Left_femur.original_landmarks[ld]['coords']
            try:
                ID = Left_femur.original_landmarks[ld]['ID']
                file_fit_original_landmarks.write(
                    str(ld) + '_left' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_original_landmarks.write(
                    str(ld) + '_left' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + '\n')
    for ld in sorted(Left_tibfib.original_landmarks):
        if ld == "lengths":
            continue
        elif ld == "angles":
            continue
        else:
            coords = Left_tibfib.original_landmarks[ld]['coords']
            try:
                ID = Left_tibfib.original_landmarks[ld]['ID']
                file_fit_original_landmarks.write(
                    str(ld) + '_left' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_original_landmarks.write(
                    str(ld) + '_left' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + '\n')

    for ld in sorted(Right_femur.original_landmarks):
        if ld == "lengths":
            continue
        elif ld == "angles":
            continue
        else:
            coords = Right_femur.original_landmarks[ld]['coords']
            try:
                ID = Right_femur.original_landmarks[ld]['ID']
                file_fit_original_landmarks.write(
                    str(ld) + '_right' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_original_landmarks.write(
                    str(ld) + '_right' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + '\n')
    for ld in sorted(Right_tibfib.original_landmarks):
        if ld == "lengths":
            continue
        elif ld == "angles":
            continue
        else:
            coords = Right_tibfib.original_landmarks[ld]['coords']
            try:
                ID = Right_tibfib.original_landmarks[ld]['ID']
                file_fit_original_landmarks.write(
                    str(ld) + '_right' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + ',' + str(ID) + '\n')
            except KeyError:
                file_fit_original_landmarks.write(
                    str(ld) + '_right' + ',' + str(coords[0]) + ',' + str(coords[1]) + ',' + str(
                        coords[2]) + '\n')