clear;
clc;
close all;

%% create model
model = createpde('structural','modal-solid');

%% set dimension
XWidth = 10; %width of the plate
YLength = XWidth; %length of the plate
ZHeight = 1; %height of the plate

%% set number of elements
numXelements = 32; %number of elements along the plate width
numYelements = 32; %number of elements along the plate length
numZelements = 20; %number of elements along the plate height

%% set mesh grid
[xg, yg, zg] = meshgrid(0:XWidth/numXelements:XWidth, ...
                        0:YLength/numYelements:YLength, ...
                        0:ZHeight/numZelements:ZHeight);
                    
%% definitions for our mesh nodes and our mesh elements
x = xg(:); %x nodes
y = yg(:); %y nodes
z = zg(:); %z nodes
K = convhull(x, y, z); %mesh elements
nodes = [x'; y'; z'];
elements = K';

%% generate the geometry of the plate from the mesh that we have created

geometryFromMesh(model, nodes, elements); %create our geometry
pdegplot(model, 'FaceLabels', 'on', 'FaceAlpha', 0.5); %plot the geometry with 
%labels on the faces for setting boundary conditions

%% material properties
structuralProperties(model, 'YoungsModulus', 200e9, ...
                             'PoissonsRatio', 0.3, ...
                             'MassDensity', 8000);

%% Boundary conditions
structuralBC(model, 'Face', [2,3,4,5],'Constraint','free');

%% add the final mesh to the model.
generateMesh(model, 'Hmin', ZHeight*2);
figure;
pdeplot3D(model); %plot the mesh
title('Mesh with Quadratic Tetrahedral Elements');

%% specify the range of frequencies
maxFreq = 1000; %maximum frequency in Hz
result = solve(model, 'FrequencyRange', [-0.1 maxFreq*2*pi]);
freqHz = result.NaturalFrequencies/(2*pi); %convert to Hz

%% display our results in a tabular format.
tfreqHz = table(freqHz);
tfreqHz.Properties.VariableNames = {'Computed'};
disp(tfreqHz);

%% Plot
h = figure;
h.Position = [100, 100, 900, 600]; %set the location and size of the figure
numToPrint = 16; %number of modes to plot
for ii = 7:numToPrint
 subplot(5, 2, ii-6);
 pdeplot3D(model, 'ColorMapData', result.ModeShapes.uz(:, ii));
 axis equal;
 title(sprintf(['Mode=%d, z-displacement\n',...
 'Frequency(Hz): FEM=%g'],...
 ii, freqHz(ii)));

end


