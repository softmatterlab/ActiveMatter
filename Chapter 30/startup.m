close all; clear

global usecamera_g video_g filename_g

addpath('subs/');

usecamera_g = 0;  % if 1: get image from camera
                  % if 0: get image from file

if usecamera_g
                 
    instrreset % disconnect all (serial) devices

    addpath('subs/')

    % initialize camera
    video_g = [];
    confCamera;
    preview(video_g);

else
    
     filename_g = 'Cameron_pendantdrop_SDM - 0_Date2021-02-18_Time16-06-12.png';
    
end