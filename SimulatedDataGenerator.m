%% Ground-truth data simulator for MAS
%% Copyright 2022 John S H Danial
%% Department of Chemistry, Univerity of Cambridge

%% clears all variables and sets path - DO NOT EDIT
clear all;
path = uigetdir;

%% camera parameters
cameraFrameSize = 1000;

%% sample parameters - EDITABLE
stoichiometry = [500 500 500 500 500];
numStructures = sum(stoichiometry);
sampleType = 'Unknown';
lateralPrecision = 0.7;

%% calling FrameGenerator
[frameCoorX,frameCoorY,intensity] = FrameGenerator(numStructures,stoichiometry,cameraFrameSize,sampleType);

%% calling StructureGenerator
rawImage = ImageGenerator(frameCoorX,frameCoorY,intensity,lateralPrecision,cameraFrameSize);

%% exporting movies
imwrite(uint16(rawImage),fullfile(path,'image.tif'),'WriteMode','overwrite');

%%==========FrameGenerator===========%%
function [frameCoorX,frameCoorY,intensity] = FrameGenerator(numStructures,stoichiometry,cameraFrameSize,sampleType)

% initializing
photonCount = 2;
varPhotonCount = 0;
marginSize = 150;

% calculating x and y coordinates of each structure
frameCoorX = marginSize + ((cameraFrameSize - 2 * marginSize) * rand(1,numStructures));
frameCoorY = marginSize + ((cameraFrameSize - 2 * marginSize) * rand(1,numStructures));

% calculating number of units in each structure
numUnits = [];
multFactor = [];
for unitId = 1 : length(stoichiometry)
    numUnits = [numUnits unitId .* ones(1,stoichiometry(unitId))];
    if strcmp(sampleType,'Unknown')
        multFactor = [multFactor ((200 * (unitId - 1) + 100) / 30) * photonCount .* ones(1,stoichiometry(unitId))];
    else
        multFactor = [multFactor photonCount .* ones(1,stoichiometry(unitId))];
    end
end

% creating random intensities
intensity = multFactor .* normrnd(double(photonCount),varPhotonCount * double(photonCount),1,length(frameCoorX));
end

%%==========ImageGenerator===========%%
function rawImage = ImageGenerator(frameCoorX,frameCoorY,intensity,lateralPrecision,cameraFrameSize)

% initializing
roiRadius = 5;

% creating kernel
for xVal = 1 : (roiRadius * 2) + 1
    for yVal = 1 : (roiRadius * 2) + 1
        kernel(yVal,xVal) = double(exp(-((xVal - (roiRadius + 1)) ^ 2 + (yVal - (roiRadius + 1)) ^ 2) /...
            (2 * ((lateralPrecision) ^ 2))));
    end
end

% creating an empty image array
rawImage = zeros(cameraFrameSize,cameraFrameSize);

% copying frameCoor to currentFrameCoor for future manipulations
currentFrameCoorX = frameCoorX;
currentFrameCoorY = frameCoorY;

% generating gaussian profiles
for unitId = 1 : size(frameCoorX,2)

    % adding kernel to movie
    try
        xVal = round(currentFrameCoorX(unitId)) - roiRadius : round(currentFrameCoorX(unitId)) + roiRadius;
        yVal = round(currentFrameCoorY(unitId)) - roiRadius : round(currentFrameCoorY(unitId)) + roiRadius;
        rawImage(yVal,xVal) = rawImage(yVal,xVal) + intensity(unitId) .* kernel;
    catch
    end
end

% creating noise
noise = gamrnd(0.15,0.9,cameraFrameSize,cameraFrameSize);
noise = noise .* (noise <= 1);
rawImage = rawImage + noise;
end