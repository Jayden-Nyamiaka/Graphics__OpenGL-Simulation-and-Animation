from difflib import Differ
import sys


def find_dif(fpath1, fpath2):
    diff = False
    with open(fpath1) as file_1, open(fpath2) as file_2:
        differ = Differ()
    
        for line in differ.compare(file_1.readlines(), file_2.readlines()):
            if (line[0] != " "):
                diff = True
                print(line)
    return diff

       
actual_dir = "interpolated_frames/"
out_dir = "output/"
name = "bunny"
ext = ".obj"
 
for i in range(0, 20):
    if (i % 5 == 0):
        i += 1
        continue
    
    frame_name = name + str(format(i, "02")) + ext
    actual_frame = actual_dir + frame_name
    out_frame = out_dir + frame_name
    
    print(i, ": ", actual_frame, " ", out_frame)
    same = find_dif(actual_frame, out_frame)
        
