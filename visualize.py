import os
from pathlib import Path
import numpy as np
import pyvista as pv

class SetVisibilityCallback:
    """Helper callback to keep a reference to the actor being modified."""

    def __init__(self, actor) -> None:
        self.actor = actor

    def __call__(self, state):
        self.actor.SetVisibility(state)

def create_checkboxes():
    checkbox_widget_size = 20
    # set some buffer between the widgets
    buffer_size = 5
    X_pos_checkbox = 5
    Y_pos_checkbox = 5
    X_pos_text = X_pos_checkbox + checkbox_widget_size + buffer_size
    Y_pos_text = 5
    # Make a callback for landmarks checkbox
    p.add_text('Landmarks on/off',position=(X_pos_text,Y_pos_text),color=other_text_color,font_size=10)
    callback = SetVisibilityCallback(landmark_actor)
    p.add_checkbox_button_widget(
            callback,
            value=True,
            position=(X_pos_checkbox,Y_pos_checkbox),
            size=20,
            border_size=1,
            color_on='green',
            color_off="red",
            background_color="grey",
        )
    
    # Make a callback for Measurements checkbox
    p.add_text('Measurements on/off',position=(30,30),color=other_text_color,font_size=10)
    callback = SetVisibilityCallback(lines_actor)
    p.add_checkbox_button_widget(
            callback,
            value=True,
            position=(5.0,  30.0),
            size=20,
            border_size=1,
            color_on='green',
            color_off="red",
            background_color="grey",
        )

    # Make a callback for Angle lines checkbox
    p.add_text('Angle Lines on/off',position=(30,55),color=other_text_color,font_size=10)
    angles_callback = SetVisibilityCallback(angles_data_actor)
    p.add_checkbox_button_widget(
            angles_callback,
            value=True,
            position=(5,55),
            size=20,
            border_size=1,
            color_on='green',
            color_off="red",
            background_color="grey",
        )
    
    def angles_lines_toggle_vis(flag):
        for actor in angles_lines_actor_arr:
            actor.SetVisibility(flag)
    
    # Make a callback for Angle data checkbox
    p.add_text('Angle Data on/off',position=(30,80),color=other_text_color,font_size=10) #change x value
    p.add_checkbox_button_widget(
            angles_lines_toggle_vis,
            value=True,
            position=(5,80),
            size=20,
            border_size=1,
            color_on='green',
            color_off="red",
            background_color="grey",
        )
    
    # Make a callback for lines checkbox
    p.add_text('Measurement Lines on/off',position=(30,105),color=other_text_color,font_size=10)
    # Since there are multiple lines, we need to toggle the visibility of all of them
    def lines_toggle_vis(flag):
        for actor in lines_actor_arr:
            actor.SetVisibility(flag)

    _ = p.add_checkbox_button_widget(lines_toggle_vis, 
                                    value=True,
                                    position=(5,105),
                                    size=20,
                                    border_size=1,
                                    color_on='green',
                                    color_off="red",
                                    background_color="grey",
                                    )

    # Make a callback for bones checkbox
    p.add_text('Bones on/off',position=(30,130),color=other_text_color,font_size=10)
    def bones_toggle_vis(flag):
        for actor in bones_mesh_actor_arr:
            actor.SetVisibility(flag)

    _ = p.add_checkbox_button_widget(bones_toggle_vis, 
                                    value=True,
                                    position=(5,130),
                                    size=20,
                                    border_size=1,
                                    color_on='green',
                                    color_off="red",
                                    background_color="grey",
                                    )

    # Define callbacks for camera position toggles
    def toggle_camera_pos_frontal_saggital(flag):
        if flag:
            p.view_zy(negative=True)
        else:
            p.view_xy()

    def toggle_camera_pos_frontal_axial(flag):
        if flag:
            p.view_zy(negative=True)
        else:
            p.view_zx()

    p.add_text('Frontal/Axial view',position=(30,155),color=other_text_color,font_size=10)
    _ = p.add_checkbox_button_widget(toggle_camera_pos_frontal_axial,
                                    value=True,
                                    position=(5,155),
                                    size=20,
                                    border_size=1,
                                    color_on='green',
                                    color_off="red",
                                    background_color="grey")
    
    p.add_text('Frontal/Saggital view',position=(30,180),color=other_text_color,font_size=10)
    _ = p.add_checkbox_button_widget(toggle_camera_pos_frontal_saggital,
                                    value=True,
                                    position=(5,180),
                                    size=20,
                                    border_size=1,
                                    color_on='green',
                                    color_off="red",
                                    background_color="grey")

def process_landmarks(landmarks_labels, landmarks_points, landmark_filtered_labels):
    plot_landmarks_lbls = []
    plot_landmarks_points = []
    for lbl in landmark_filtered_labels:
        idx = landmarks_labels.index(lbl)
        pnt = landmarks_points[idx]
        plot_landmarks_lbls.append(lbl)
        plot_landmarks_points.append(pnt)

    return plot_landmarks_lbls, plot_landmarks_points
    
def process_measurements_lines(plot_landmarks_lbls, plot_landmarks_points, measurements_labels, measurements_data, measurement_filtered_labels):   
    # Construct the line co-ordinates to draw on the mesh
    lines = []
    lines_labels = []
    lines_data = []
    lines_points = []

    for lbl in measurement_filtered_labels:
        idx = measurements_labels.index(lbl)
        if lbl == 'ASIS_width':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('RASIS')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('LASIS')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[2] += -80
            lines_points.append(temp)

        if lbl == 'PSIS_width':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('RPSIS')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('LPSIS')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[2] += -45
            lines_points.append(temp)
        
        if lbl == 'pelvis_depth':
            # get corresponding points from landmarks
            #Caclulate mid point of ASIS and PSIS since they are not being displayed
            #p1 = plot_landmarks_points[plot_landmarks_lbls.index('ASIS_mid')]
            #p2 = plot_landmarks_points[plot_landmarks_lbls.index('PSIS_mid')]

            p1 = plot_landmarks_points[plot_landmarks_lbls.index('RASIS')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('LASIS')]
            # get mid point of line
            ASIS_mid = [(p1[0] + p2[0])/2, (p1[1] + p2[1])/2, (p1[2] + p2[2])/2]
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('RPSIS')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('LPSIS')]
            # get mid point of line
            PSIS_mid = [(p1[0] + p2[0])/2, (p1[1] + p2[1])/2, (p1[2] + p2[2])/2]
            
            lines.append(pv.Line(ASIS_mid,PSIS_mid))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = ASIS_mid.copy()
            temp[0] += -80
            lines_points.append(temp)

        if lbl == 'left_epicon_width':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('med_epicon_left')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('lat_epicon_left')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[2] += -40
            lines_points.append(temp)

        if lbl == 'right_epicon_width':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('med_epicon_right')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('lat_epicon_right')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[2] += 40
            lines_points.append(temp)
        
        if lbl == 'left_femoral_length':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('lat_epicon_left')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('lat_great_trochant_left')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[1] += 80
            lines_points.append(temp) 

        if lbl == 'right_femoral_length':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('lat_epicon_right')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('lat_great_trochant_right')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[1] += 80
            lines_points.append(temp) 

        if lbl == 'left_condylar_width':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('medial_condyle_left')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('lateral_condyle_left')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[2] += -40
            lines_points.append(temp)

        if lbl == 'right_condylar_width':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('medial_condyle_right')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('lateral_condyle_right')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[2] += 40
            lines_points.append(temp)

        if lbl == 'left_malleolar_width':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('medial_malleolus_left')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('lateral_malleolus_left')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[2] += -15
            temp[0] += -15
            lines_points.append(temp)

        if lbl == 'right_malleolar_width':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('medial_malleolus_right')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('lateral_malleolus_right')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[2] += 15
            temp[0] += -18
            lines_points.append(temp)

        if lbl == 'left_tibial_length':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('lateral_malleolus_left')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('lateral_condyle_left')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[1] += 80
            lines_points.append(temp)

        if lbl == 'right_tibial_length':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('lateral_malleolus_right')]
            p2 = plot_landmarks_points[plot_landmarks_lbls.index('lateral_condyle_right')]
            lines.append(pv.Line(p1,p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')

            #add offset to make display label more visible
            temp = p1.copy()
            temp[1] += 80
            lines_points.append(temp)

        # Femoral head diameter
        if lbl == 'left_FHC_diameter':
            # get corresponding points from landmarks
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('FHC_left')]
            dia = float(measurements_data[idx][:6])
            line_p1 = p1.copy()
            line_p2 = p1.copy()
            line_p1[0] -= dia/2
            line_p2[0] += dia/2
            lines.append(pv.Line(line_p1, line_p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')
            temp = p1.copy()
            temp[0] += 10
            lines_points.append(temp)
            #p.add_mesh(pv.Sphere(dia/2, p1), color='orange', show_edges=False, opacity=0.9)

        if lbl == 'left_FHC_diameter':
            p1 = plot_landmarks_points[plot_landmarks_lbls.index('FHC_right')]
            dia = float(measurements_data[idx][:6])
            line_p1 = p1.copy()
            line_p2 = p1.copy()
            line_p1[0] -= dia/2
            line_p2[0] += dia/2
            lines.append(pv.Line(line_p1, line_p2))
            lines_labels.append(measurements_labels[idx])
            lines_data.append(lbl + ':' + measurements_data[idx][:6] + 'mm')
            temp = p1.copy()
            temp[0] += 10
            lines_points.append(temp)
            #p.add_mesh(pv.Sphere(dia/2, p1), color='orange', show_edges=False, opacity=0.9)

    return lines, lines_labels, lines_data, lines_points

def process_angles(landmarks_labels, landmarks_points, measurements_labels, measurements_data, measurement_filtered_labels):
    def find_equidistant_points_diverging(center, p1, p2, p3, p4, distance):
        def calculate_direction(start, end):
            # Calculate the direction vector from start to end
            direction = end - start
            # Normalize the direction vector
            direction = direction / np.linalg.norm(direction)
            return direction
        p1 = np.array(p1)
        p2 = np.array(p2)
        p3 = np.array(p3)
        p4 = np.array(p4)
        direction1 = calculate_direction(p1, p2)
        direction2 = calculate_direction(p3, p4)

        # Normalize the direction vectors
        direction1 = direction1 / np.linalg.norm(direction1)
        direction2 = direction2 / np.linalg.norm(direction2)
        
        # Calculate the points along each line at the specified distance
        point1 = center + distance * direction1
        point2 = center + distance * direction2
        
        return point1, point2

    def extend_line(p1, p2, distance):
        # Convert points to numpy arrays
        p1 = np.array(p1)
        p2 = np.array(p2)
        
        # Calculate the direction vector from p1 to p2
        direction = p2 - p1
        
        # Normalize the direction vector to unit length
        unit_direction = direction / np.linalg.norm(direction)
        
        # Calculate the new external point
        new_point = p2 + unit_direction * distance
        
        return tuple(new_point)

    angles_lines = []
    angles_labels = []
    angles_data = []
    angles_points = []
    arcs = []
    for lbl in measurement_filtered_labels:
        idx = measurements_labels.index(lbl)
        if lbl == 'left_AA':
            # get corresponding points from landmarks
            p1 = landmarks_points[landmarks_labels.index('FHC_left')]
            p2 = landmarks_points[landmarks_labels.index('neck_centre_left')]
            p3 = landmarks_points[landmarks_labels.index('med_bot_con_left')]
            p4 = landmarks_points[landmarks_labels.index('lat_bot_con_left')]
             
            pe = extend_line(p1, p2, 550)
            pe2 = extend_line(p3, p4, 200)

            angles_lines.append(pv.Line(p2,pe))
            angles_lines.append(pv.Line(p4,pe2))
            
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_data.append(lbl + ':' + measurements_data[idx][:6] + '°')

            # Create arc
            arc_a, arc_b = find_equidistant_points_diverging(pe2,pe2,p4,pe,p2,50)
            arcs.append(pv.CircularArc(arc_a,arc_b,pe2))

            # Add offset for better visualization
            arc_off = (arc_a[0], arc_a[1] + 10, arc_a[2] + 50)
            angles_points.append(arc_off)

        if lbl == 'right_AA':
            # get corresponding points from landmarks
            p1 = landmarks_points[landmarks_labels.index('FHC_right')]
            p2 = landmarks_points[landmarks_labels.index('neck_centre_right')]
            p3 = landmarks_points[landmarks_labels.index('med_bot_con_right')]
            p4 = landmarks_points[landmarks_labels.index('lat_bot_con_right')]
             
            pe = extend_line(p1, p2, 550)
            pe2 = extend_line(p3, p4, 250)

            angles_lines.append(pv.Line(p2,pe))
            angles_lines.append(pv.Line(p4,pe2))
            
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_data.append(lbl + ':' + measurements_data[idx][:6] + '°')

            # Create arc
            arc_a, arc_b = find_equidistant_points_diverging(pe2,pe2,p4,pe,p2,50)
            arcs.append(pv.CircularArc(arc_a,arc_b,pe2))

            # Add offset for better visualization
            arc_off = (arc_a[0], arc_a[1] + 10, arc_a[2] - 20)
            angles_points.append(arc_off)

        if lbl == 'left_NSA':
            # get corresponding points from landmarks
            p1 = landmarks_points[landmarks_labels.index('FHC_left')]
            p2 = landmarks_points[landmarks_labels.index('neck_centre_left')]
            p3 = landmarks_points[landmarks_labels.index('shaft_centre_left')]
            p4 = landmarks_points[landmarks_labels.index('epicon_mid_left')]

            angles_lines.append(pv.Line(p1,p2))
            angles_lines.append(pv.Line(p2,p4))
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_data.append(lbl + ':' + measurements_data[idx][:6] + '°')

            # Create arc
            arc_a, arc_b = find_equidistant_points_diverging(p2,p2,p1,p2,p3,3)
            arcs.append(pv.CircularArc(arc_a,arc_b,p2, negative=False))
            # Add offset for better visualization
            arc_off = (arc_a[0]-2, arc_a[1] - 5, arc_a[2] - 1)
            angles_points.append(arc_off)

        if lbl == 'right_NSA':
            # get corresponding points from landmarks
            p1 = landmarks_points[landmarks_labels.index('FHC_right')]
            p2 = landmarks_points[landmarks_labels.index('neck_centre_right')]
            p3 = landmarks_points[landmarks_labels.index('shaft_centre_right')]
            p4 = landmarks_points[landmarks_labels.index('epicon_mid_right')]

            angles_lines.append(pv.Line(p1,p2))
            angles_lines.append(pv.Line(p2,p4))
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_data.append(lbl + ':' + measurements_data[idx][:6] + '°')

            # Create arc
            arc_a, arc_b = find_equidistant_points_diverging(p2,p2,p1,p2,p3,3)
            arcs.append(pv.CircularArc(arc_a,arc_b,p2, negative=False))
            # Add offset for better visualization
            arc_off = (arc_a[0]-2, arc_a[1] - 5, arc_a[2] - 1)
            angles_points.append(arc_off)

        if lbl == 'left_mLDFA':
            # get corresponding points from landmarks
            p1 = landmarks_points[landmarks_labels.index('med_bot_con_left')]
            p2 = landmarks_points[landmarks_labels.index('lat_bot_con_left')]
            # mid point of p1 and p2
            p3 = [(p1[0] + p2[0])/2, (p1[1] + p2[1])/2, (p1[2] + p2[2])/2]
            p4 = landmarks_points[landmarks_labels.index('FHC_left')]
            #p4 = landmarks_points[landmarks_labels.index('epicon_mid_left')]

            angles_lines.append(pv.Line(p1,p2))
            angles_lines.append(pv.Line(p3,p4))
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_data.append(lbl + ':' + measurements_data[idx][:6] + '°')

            # Create arc
            arc_a, arc_b = find_equidistant_points_diverging(p3,p1,p2,p3,p4,3)
            arcs.append(pv.CircularArc(arc_b,arc_a,p3))
            # Add offset for better visualization
            arc_off = (arc_a[0], arc_a[1] + 2, arc_a[2] - 1)
            angles_points.append(arc_off)

        if lbl == 'right_mLDFA':
            # get corresponding points from landmarks
            p1 = landmarks_points[landmarks_labels.index('med_bot_con_right')]
            p2 = landmarks_points[landmarks_labels.index('lat_bot_con_right')]
            # mid point of p1 and p2
            p3 = [(p1[0] + p2[0])/2, (p1[1] + p2[1])/2, (p1[2] + p2[2])/2]
            p4 = landmarks_points[landmarks_labels.index('FHC_right')]
            #p4 = landmarks_points[landmarks_labels.index('epicon_mid_right')]

            angles_lines.append(pv.Line(p1,p2))
            angles_lines.append(pv.Line(p3,p4))
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_data.append(lbl + ':' + measurements_data[idx][:6] + '°')

            # Create arc
            arc_a, arc_b = find_equidistant_points_diverging(p3,p2,p1,p3,p4,3)
            arcs.append(pv.CircularArc(arc_b,arc_a,p3))
            # Add offset for better visualization
            arc_off = (arc_a[0], arc_a[1] + 2, arc_a[2] - 3)
            angles_points.append(arc_off)

        if lbl == 'left_TT':
            # get corresponding points from landmarks
            p1 = landmarks_points[landmarks_labels.index('medial_condyle_left')]
            p2 = landmarks_points[landmarks_labels.index('lateral_condyle_left')]
            #p3 = landmarks_points[landmarks_labels.index('medial_malleolus_left')]
            #p4 = landmarks_points[landmarks_labels.index('lateral_malleolus_left')]
            p4 = landmarks_points[landmarks_labels.index('medial_malleolus_left')]

            angles_lines.append(pv.Line(p1,p2))
            angles_lines.append(pv.Line(p2,p4))
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_data.append(lbl + ':' + measurements_data[idx][:6] + '°')

            # Create arc
            arc_a, arc_b = find_equidistant_points_diverging(p2,p2,p1,p2,p4,3)
            arcs.append(pv.CircularArc(arc_a,arc_b,p2))
            # Add offset for better visualization
            arc_off = (arc_a[0], arc_a[1] - 2, arc_a[2] + 5)
            angles_points.append(arc_off)
        
        if lbl == 'right_TT':
            # get corresponding points from landmarks
            p1 = landmarks_points[landmarks_labels.index('medial_condyle_right')]
            p2 = landmarks_points[landmarks_labels.index('lateral_condyle_right')]
            #p3 = landmarks_points[landmarks_labels.index('medial_malleolus_left')]
            #p4 = landmarks_points[landmarks_labels.index('lateral_malleolus_left')]
            p4 = landmarks_points[landmarks_labels.index('medial_malleolus_right')]

            angles_lines.append(pv.Line(p1,p2))
            angles_lines.append(pv.Line(p2,p4))
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_data.append(lbl + ':' + measurements_data[idx][:6] + '°')

            # Create arc
            arc_a, arc_b = find_equidistant_points_diverging(p2,p2,p1,p2,p4,3)
            arcs.append(pv.CircularArc(arc_a,arc_b,p2))
            # Add offset for better visualization
            arc_off = (arc_a[0], arc_a[1] - 2, arc_a[2] + 2)
            angles_points.append(arc_off)

        if lbl == 'left_mMPTA':
            # get corresponding points from landmarks
            p1 = landmarks_points[landmarks_labels.index('tib_plateau_med_left')]
            p2 = landmarks_points[landmarks_labels.index('tib_plateau_lat_left')]
            # mid point of p1 and p2
            p3 = [(p1[0] + p2[0])/2, (p1[1] + p2[1])/2, (p1[2] + p2[2])/2]
            p4 = landmarks_points[landmarks_labels.index('int_mal_left')]
            #p4 = landmarks_points[landmarks_labels.index('con_mid_left')]

            angles_lines.append(pv.Line(p1,p2))
            angles_lines.append(pv.Line(p3,p4))
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_data.append(lbl + ':' + measurements_data[idx][:6] + '°')

            # Create arc
            arc_a, arc_b = find_equidistant_points_diverging(p3,p2,p1,p3,p4,3)
            arcs.append(pv.CircularArc(arc_a,arc_b,p3))
            # Add offset for better visualization
            arc_off = (arc_a[0], arc_a[1] - 2, arc_a[2])
            angles_points.append(arc_off)

        if lbl == 'right_mMPTA':
            # get corresponding points from landmarks
            p1 = landmarks_points[landmarks_labels.index('tib_plateau_med_right')]
            p2 = landmarks_points[landmarks_labels.index('tib_plateau_lat_right')]
            # mid point of p1 and p2
            p3 = [(p1[0] + p2[0])/2, (p1[1] + p2[1])/2, (p1[2] + p2[2])/2]
            p4 = landmarks_points[landmarks_labels.index('int_mal_right')]
            #p4 = landmarks_points[landmarks_labels.index('con_mid_left')]

            angles_lines.append(pv.Line(p1,p2))
            angles_lines.append(pv.Line(p3,p4))
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_labels.append(measurements_data[idx][:6] + '°')
            angles_data.append(lbl + ':' + measurements_data[idx][:6] + '°')

            # Create arc
            arc_a, arc_b = find_equidistant_points_diverging(p3,p2,p1,p3,p4,3)
            arcs.append(pv.CircularArc(arc_a,arc_b,p3))
            # Add offset for better visualization
            arc_off = (arc_a[0], arc_a[1] - 2, arc_a[2])
            angles_points.append(arc_off)

    return angles_lines, angles_labels, angles_data, angles_points, arcs

patient_number = 'TD5'

root_dir = Path('.\\data_mean\\')

# get landmarks and measurement data file paths
l_path = Path('data_mean/' + 'landmarks/' + patient_number + '_landmarks.txt')

#Patient number is missing from the measurements file name hence excluded
m_path = Path('data_mean/' + 'measurements/' + 'measurements.txt')

# load landmarks and measurements
landmarks = np.loadtxt(l_path, dtype=object)
measurements = np.loadtxt(m_path, dtype=object)

measurements_labels = measurements[0].split(',')
measurements_data = measurements[1] + measurements[2]
measurements_data = measurements_data.split(',')

# using the processed and aligned meshes after running the find_landmarks_and_align.py script
mesh_dir = Path('.\\aligned_meshes_CS_3\\')

meshes = []
for file in os.listdir(mesh_dir):
    if file.endswith('.ply') and patient_number in file:
        meshes.append(pv.read(mesh_dir / file))

landmarks_points = []
landmarks_labels = []

#row 0 is the header
for row in landmarks[1:]:    
    splits = row.split(',')
    landmarks_points.append([float(x) for x in splits[1:4]])
    landmarks_labels.append(splits[0])

angles = []
angles_labels = []
angles_data = []
angles_points = []

landmark_filtered_labels = np.loadtxt('landmarks_labels_filtered.txt', dtype=str)
measurement_filtered_labels = np.loadtxt('measurements_labels_filtered.txt', dtype=str)

plot_landmarks_lbls, plot_landmarks_points = process_landmarks(landmarks_labels, 
                                                               landmarks_points, 
                                                               landmark_filtered_labels)

lines, lines_labels, lines_data, lines_points = process_measurements_lines(plot_landmarks_lbls, 
                                                                           plot_landmarks_points, 
                                                                           measurements_labels, 
                                                                           measurements_data, 
                                                                           measurement_filtered_labels)

angles_lines, angles_labels, angles_data, angles_points, arcs = process_angles(landmarks_labels, 
                                                                                landmarks_points, 
                                                                                measurements_labels, 
                                                                                measurements_data, 
                                                                                measurement_filtered_labels)

p = pv.Plotter(lighting = 'light kit', theme=pv.set_plot_theme('default'))
label_text_color = 'white'
other_text_color = 'black'

# Bones meshes
bones_mesh_actor_arr = [] 
for mesh in meshes:
    bones_mesh_actor_arr.append(p.add_mesh(mesh, color='white', show_edges=False, opacity=0.99))

# Plots landmarks
landmark_actor = p.add_point_labels(plot_landmarks_points, plot_landmarks_lbls, text_color=label_text_color, shape_color = 'green', background_color = 'green', font_size = 12, always_visible=True, render_points_as_spheres=True, point_size=8, pickable=True)
lines_actor_arr = []

# Plots lines
for idx, line in enumerate(lines):
    lines_actor_arr.append(p.add_mesh(line, color='red', label=lines_labels[idx], line_width=3))

# Plots line measurements and angle measurements
lines_actor = p.add_point_labels(lines_points, lines_data, font_size = 12, always_visible=False, text_color = label_text_color, background_color = 'blue', shape_color = 'blue', render_points_as_spheres=False)

# Plots angles lines
angles_lines_actor_arr = []
angles_data_actor_arr = []
for idx, angles_line in enumerate(angles_lines):
    angles_lines_actor_arr.append(p.add_mesh(angles_line, color='orange', label=angles_labels[idx], line_width=3))

angles_data_actor =  p.add_point_labels(angles_points, angles_data, font_size = 12, always_visible=True, text_color = label_text_color, background_color = 'blue', shape_color = 'blue', render_points_as_spheres=False)

for arc in arcs:
    p.add_mesh(arc, color='black', line_width=3)

# Data table
data_table_text = ''

data_table_text += '---- Measurements ----\n\n'
for item in lines_data:
    splits = item.split(':')
    data_table_text += splits[0] + ' = ' + splits[1] + '\n'

data_table_text += '\n ---- Angles ----\n\n'
for item in angles_data:
    splits = item.split(':')
    data_table_text += splits[0] + ' = ' + splits[1] + '\n'

actor = p.add_text(data_table_text,position='upper_right',color=other_text_color,font_size=10)

#Set intial view to frontal view
p.view_zy(negative=True)
p.add_axes(labels_off=False)

create_checkboxes()

p.show()
