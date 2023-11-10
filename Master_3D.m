%% modelling cross sections
addpath(genpath('Y:\DASH_v1.0\'));
casefolder = 'Y:\DASH_v1.0\hsf_01\'; 
casefoldername = 'hsf_01';
outputFolder_simulation = [casefolder,'Simulation_1\'];
parameterfilename = 'INPUTS_parameters.m';


run('VerticalSectionTopoLine.m');
run('HorizontalCountourLine.m');
run('CalTrueDip4CrossSection_xz.m');

%% Generate 3D images
H = 218; % the height of image in the number of pixels
ImageDis = 995; % the width of image in the number of pixels
INum = 200; % the number of cross sections
trimEdge = 10; % the number of pixels that have been excluded for generating 3D images

run('image2data_3D_CSleft_new.m');
run('image2data_3D_CSright_new.m');
run('image2data_3D_horizontal_new.m');

%% Output modelling results to .txt file (ASCII or GSLib )
run('Data2TextFile.m');

%% create fence diagram
run('INPUTS_parameters.m');
crosssectionBetweenstep = 20; % specify the distance between cross sections
az = -38; % azimuth of graph
el = 50; % elevation of graph 
z_axis_exaggerated = 1.5; % vertical exaggeration, times of the x or y axis
run('fenceDiagram.m');