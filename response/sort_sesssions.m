

% This script is used to sort out the stimulation sessions that contained
% an artefact. For usage:

% Run the syn_analysis script with the config.debug_images option set to 1
% this will create a average window plot for every stimulation sesssion.
% place this script into the folder and run it. Then you can rate every
% response as you like. In the synanalysis you can afterwards sort out
% sessions wich you rated with a certain number.. 

% I usually saved my flags in a handsorted.mat in the results/response
% folder


% sort out the bad results (noisy channels, artefact and co)
% 1 is good
% 2 is with an artefact
% 3 is no response

pics = dir('*.tiff');

for i = 1:length(pics)
    im = imread(pics(i).name);
    figure(1);
    image(im);
    flags(i) = input('wenn ja dann 1: ');
end