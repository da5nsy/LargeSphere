% This script calls BaselineDataGenerator with a range of different values
% for slider values and the 'current centres' (these are normally randomly
% distributed between 0.25 and 0.75).

clear, clc, close all

slider_vals = ...
       [500,    500;
        0,      0;
        1000,   0;
        0,      1000;
        1000,   1000];
    
current_centres = ...
       [0.5,    0.5;
        0.25,   0.25;
        0.75,   0.25;
        0.25,   0.75;
        0.75,   0.75];
    
%%

for sv = 1:5
    for cc = 1:5
        BaselineDataGenerator(slider_vals(sv,1),slider_vals(sv,2),current_centres(cc,1),current_centres(cc,2));
    end
end

    
    