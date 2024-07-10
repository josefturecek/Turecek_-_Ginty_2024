# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 07:14:04 2023

@author: Josef
"""

import cv2
import os

folder_path = "Z:/GintyLab/Turecek/data/2024_04_23_kx132" 

os.chdir(folder_path)

def __draw_label(img, text, pos, bg_color):
   font_face = cv2.FONT_HERSHEY_SIMPLEX
   scale = 1
   color = (0, 0, 0)
   thickness = cv2.FILLED
   margin = 5
   txt_size = cv2.getTextSize(text, font_face, scale, thickness)
   end_x = pos[0] + txt_size[0][0] + margin
   end_y = pos[1] - txt_size[0][1] - margin
   cv2.rectangle(img, pos, (end_x, end_y), bg_color, thickness)
   cv2.putText(img, text, pos, font_face, scale, color, 1, cv2.LINE_AA)


# get all filenames in folder
filenames = os.listdir(folder_path)

filenames = sorted(filenames)

# filter filenames by extension and string
avi_filenames = [filename for filename in filenames if filename.endswith(".avi") and "C2" in filename and not "SYNCD" in filename]

img_array = []

for filename in avi_filenames:
    video = cv2.VideoCapture(filename)
    
    video.set(cv2.CAP_PROP_POS_FRAMES, 5)
    ret, frame = video.read()
    height,width,layers=frame.shape
    
    __draw_label(frame, filename, (20,20), (255,255,255))
    
    # img = cv2.imread(frame)
        
    # cv2.imshow('frame',frame)
    # cv2.waitKey(0)
    img_array.append(frame)
    
fourcc = cv2.VideoWriter_fourcc('M','J','P','G')
out = cv2.VideoWriter('2023_04_23_compilation.avi', fourcc, 3, (width, height))
 
for i in range(len(img_array)):
    out.write(img_array[i])
out.release()

home = os.path.expanduser('~')
os.chdir(home)

# img_array = []
# for filename in glob.glob('C:/New folder/Images/*.jpg'):
#     img = cv2.imread(filename)
#     height, width, layers = img.shape
#     size = (width,height)
#     img_array.append(img)
 
 
# create jpgs 

# for file in os.listdir("C:/temp/2023_03_06_AB_HRA"):
#     if file.endswith(".avi"):
#         if f.read(".avi"):
#         print(file)
        
#         import os
# for root, dirs, files in os.walk("/mydir"):
#     for file in files:
#         if file.endswith(".txt"):
#              print(os.path.join(root, file))