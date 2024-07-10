# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 12:25:05 2023

@author: JTurecek
"""

import cv2
import os
import numpy as np
import matplotlib.pyplot as plt

from pywinauto.findwindows import find_window
from win32gui import SetForegroundWindow

#### WARNING --- MUST TYPE THIS INTO COMMAND LINE BEFORE BEGINNING:
####                     %matplotlib qt
###################################

label =[]

##################################################################################
#inputs:
label.append('swing_acrylic')  #behavior label
#file_no = np.arange(147,164)
file_no = np.arange(37,42)
# file_no = ['0034','0035','0036','0037','0038']
folder_path = 'Z:/GintyLab/Turecek/data/2023_03_23_AB_proprioceptor'
    
##################################################################################


file_noSr = ['0' + str(num) for num in file_no]

avi_files = [filename for filename in os.listdir(folder_path)
             if filename.endswith('.avi') and any(number in filename for number in file_noSr)
             and 'COMBINED' in filename]

avi_files = sorted(avi_files)

os.chdir(folder_path)


for filename in avi_files:


    # Open the video file
    cap = cv2.VideoCapture(filename)

    # os.chdir('C:/temp/2023_03_06_AB_HRA')

    # # Open the video file
    # cap = cv2.VideoCapture('2023_03_06_0061_C2-03062023153411-0000.avi')

    # Create filename for output list
    
    fileAdd = []
    fileAdd.append(filename)
    keep_frames = [];
    labelAdd = []
    labelAdd.append(label)

    length = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))

    # create video window with trackbar
    def onChange(trackbarValue):
        cap.set(cv2.CAP_PROP_POS_FRAMES,trackbarValue)
        err,img = cap.read()
        cv2.imshow(filename, img)
        pass

    cv2.namedWindow(filename)
    cv2.createTrackbar( 'frame', filename, 0, length, onChange )

    onChange(1)
    cv2.waitKey()
    frame = cv2.getTrackbarPos('frame',filename)

    # plt.figure(figsize=(12, 1))
    # 

    # create capture frame plot and set location
    mngr = plt.get_current_fig_manager()
    mngr.window.setGeometry(50,50,800, 100)

    # set limits, hide X and Y axes label marks
    plt.xlim(0, 3150)
    ax = plt.gca()
    ax.yaxis.set_tick_params(labelleft=False)
    ax.set_yticks([])
    plt.show()
    plt.ion()

    # fig, ax = plt.subplots()
    # plt.imshow(cv2.cvtColor(image, cv2.COLOR_BGR2RGB))
    
    p = 1
    
    while p == 1:
        cap.set(cv2.CAP_PROP_POS_FRAMES,frame)
        frame = cv2.getTrackbarPos('frame',filename)
        err,img = cap.read()
        # ax1.imshow(cv2.cvtColor(img, cv2.COLOR_BGR2RGB))
        cv2.imshow(filename, img)
        
        SetForegroundWindow(find_window(title=filename))
        
        k = cv2.waitKey(0) 
    
        if k == ord('g'):
            F = frame
            print(F)
            plt.plot(F,1, color='red',marker=11)
            plt.text(F,1,F)
            plt.draw()
            keep_frames.append(F)
            if len(keep_frames)>1:
                fileAdd.append(fileAdd[0])
                labelAdd.append(label)
            
            
        if k == ord('s'):
            frame += 1
            cv2.setTrackbarPos('frame',filename, frame)
            
        if k == ord('S'):
            frame += 2
            cv2.setTrackbarPos('frame',filename, frame)
        
        if k == ord('a'):
            frame -= 1
            cv2.setTrackbarPos('frame',filename, frame)
            
        if k == ord('A'):
            frame -= 2
            cv2.setTrackbarPos('frame',filename, frame)
            
        if k == ord('w'):
            frame += 10
            cv2.setTrackbarPos('frame',filename, frame)
            
        if k == ord('W'):
            frame += 100
            cv2.setTrackbarPos('frame',filename, frame)
            
        if k == ord('q'):
            frame -= 10
            cv2.setTrackbarPos('frame',filename, frame)
            
        if k == ord('Q'):
            frame -= 100
            cv2.setTrackbarPos('frame',filename, frame)
    
        if k == ord('x'):
            cap.release()
            cv2.destroyAllWindows()
            plt.close()
            if len(keep_frames)>0:
                data = np.column_stack([fileAdd, labelAdd, keep_frames])
                datafile_path = filename[0:16] + label[0] + '.txt'
                np.savetxt(datafile_path , data, fmt='"%s"')
            p = 2



home = os.path.expanduser('~')
os.chdir(home)