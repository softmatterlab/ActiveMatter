close all; clear

global usecamera_g video_g folder_g

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
    
    % initialize Arduino
    confArduino;

else
    
    folder_g = 'test_figures/water in semi-purified HD/';
    
end