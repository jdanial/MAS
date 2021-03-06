function returnFlag = detectFn(app)
% detectFn() -
% detect diffraction-limited blobs in Field of View in a 4D stack of LSM
% confocal images.
%
% Syntax -
% detectFn(app).
%
% Parameters -
% - app: MAS UI class

%% initializing returnFlag
returnFlag = false;

%% issuing initial error statements
if isempty(app.param.paths.calibrationAndUnknownData)
    returnFlag = true;
    app.msgBox.Value = sprintf('%s','Error: no data path selected.');
    return;
end

%% checking availability of calibration and unknown data
cd(fullfile(app.param.paths.calibrationAndUnknownData,'Calibration'));
dataFileList = [dir('*.lsm'); dir('*.tif')];
if isempty(dataFileList)
    returnFlag = true;
    app.msgBox.Value = sprintf('%s','Error: no calibration data in selected patb.');
    return;
end
cd(fullfile(app.param.paths.calibrationAndUnknownData,'Unknown'));
dataFileList = [dir('*.lsm'); dir('*.tif')];
if isempty(dataFileList)
    returnFlag = true;
    app.msgBox.Value = sprintf('%s','Error: no unknown data in selected patb.');
    return;
end

%% displaying progress
app.msgBox.Value = sprintf('%s','Detection started.');

%% looping though both folders, calibrtion (1) and unknown data (2)
for roundId = 1 : 2

    %% initializing struct
    app.data = struct();

    %% reading input files
    if roundId == 1
        cd(fullfile(app.param.paths.calibrationAndUnknownData,'Calibration'));
    else
        cd(fullfile(app.param.paths.calibrationAndUnknownData,'Unknown'));
    end
    inputFiles = [dir('*.lsm'); dir('*.tif')];

    %% reading number of files
    numFiles = numel(inputFiles);

    %% looping through files
    for fileId = 1 : numFiles

        % registering state
        app.data.file(fileId).state = 'detected';

        % registering type
        if roundId == 1
            app.data.file(fileId).type = 'Calibration';
        else
            app.data.file(fileId).type = 'Unknown';
        end

        %% extracting file path
        if numFiles == 1
            fileName = inputFiles.name;
            fileFolder = inputFiles.folder;
        else
            fileName = inputFiles(fileId).name;
            fileFolder = inputFiles(fileId).folder;
        end
        filePath = fullfile(fileFolder,fileName);
        app.data.file(fileId).name = fileName;

        %% reading first frame
        try
            [~,~,ext] = fileparts(app.data.file(fileId).name);
            if ext == ".lsm"
                image = lsmread(filePath);
            else
                frameId = 1;
                for tId = 1 : app.param.detection.timeLength
                    for zId = 1 : app.param.detection.depths
                        for cId = 1 : app.param.detection.colors
                            image(tId,cId,zId,:,:) = imread(filePath,frameId);
                            frameId = frameId + 1;
                        end
                    end
                end
            end
        catch
            returnFlag = true;
            app.msgBox.Value = sprintf('%s',['Error: cannot read image file (' fileName ').']);
            return;
        end

        %% adding image to particleData struct
        app.data.file(fileId).image = double(image);

        %% setting up progress
        app.msgBox.Value = sprintf('%s',['Processing stack in ' app.data.file(fileId).type ' file ' num2str(fileId) ' out of ' num2str(numFiles)]);
        drawnow;

        %% calling StackProcessor
        StackProcessor(app,fileId);

        %% setting up progress
        app.msgBox.Value = sprintf('%s',['Detecting particles in ' app.data.file(fileId).type ' file ' num2str(fileId)]);
        drawnow;

        %% calling ParticleDetector
        ParticleDetector(app,fileId);
    end

    %% exporting generic .m files
    exportFn(app);

    %% exporting files relevant to this function
    SpecificExport(app);
end

%% displaying progress
app.msgBox.Value = sprintf('%s','Detection complete.');
end

%%====================StackProcessor=====================%%
function StackProcessor(app,fileId)

% looping through time
for tId = 1 : size(app.data.file(fileId).image,1)

    if app.param.detection.depths > 1
        
        % calculating max projection
        app.data.file(fileId).time(tId).maxImage = squeeze(max(app.data.file(fileId).image(tId,1,:,:,:)));
    else
        app.data.file(fileId).time(tId).maxImage = squeeze(app.data.file(fileId).image(tId,1,:,:,:));
    end
end
end

%%====================ParticleDetector=====================%%
function ParticleDetector(app,fileId)

% extracting parameters
roiRadius = app.param.detection.roiRadius;

% looping through time
for tId = 1 : size(app.data.file(fileId).image,1)

    % extracting image
    imageRaw = int16(app.data.file(fileId).time(tId).maxImage);

    % checking image type
    if strcmp(app.data.file(fileId).type,'Calibration')
        imageDOG = imageRaw;
        threshold = 2;
    else
        sigmaLow = 1;
        sigmaHigh = 2;
        threshold = 1;

        % applying DOG filter to image
        filterSize = [roiRadius * 2 + 1,roiRadius * 2 + 1];
        kernelLow = fspecial('gaussian',filterSize,sigmaLow);
        kernelHigh = fspecial('gaussian',filterSize,sigmaHigh);
        gaussLow = imfilter(imageRaw,kernelLow,'replicate');
        gaussHigh = imfilter(imageRaw,kernelHigh,'replicate');
        imageDOG = int16(gaussLow - gaussHigh);
    end

    % segmenting cell image
    imageRawTemp = imageRaw .* int16(imageDOG > threshold);

    % changing to imageRaw to double
    imageRaw = double(imageRaw);

    % checking if image does not have any particles
    if ~sum(imageRawTemp,'all') == 0

        % thresholding image
        thresholdImage = imregionalmax(imageRawTemp,8);

        % extracting connected objects from threshold image
        [centroidTemp_y,centroidTemp_x] = ind2sub(size(thresholdImage),find(thresholdImage == true));
        centroids = [centroidTemp_x centroidTemp_y];

        % accepting or rejecting particles based on classification
        trueParticlePosCount = 0;
        for particleId = 1 : size(centroids,1)
            if centroids(particleId,1) > 2 + roiRadius && ...
                    centroids(particleId,1) < size(imageRaw,1) - roiRadius - 1 && ...
                    centroids(particleId,2) > 2 + roiRadius && ...
                    centroids(particleId,2) < size(imageRaw,2) - roiRadius - 1
                trueParticlePosCount = trueParticlePosCount + 1;
                app.data.file(fileId).time(tId).particle(trueParticlePosCount).state = 'accepted';
                app.data.file(fileId).time(tId).particle(trueParticlePosCount).centroid.x = centroids(particleId,1);
                app.data.file(fileId).time(tId).particle(trueParticlePosCount).centroid.y = centroids(particleId,2);
                background = [];
                for rowId = round(centroids(particleId,2)) - roiRadius - 1 : round(centroids(particleId,2)) + roiRadius + 1
                    for colId = round(centroids(particleId,1)) - roiRadius - 1 : round(centroids(particleId,1)) + roiRadius + 1
                        if (rowId < round(centroids(particleId,2)) - roiRadius || rowId > round(centroids(particleId,2)) + roiRadius) && ...
                                (colId < round(centroids(particleId,1)) - roiRadius || colId > round(centroids(particleId,1)) + roiRadius)
                            background = [background imageRaw(rowId,colId)];
                        end
                    end
                end
                app.data.file(fileId).time(tId).particle(trueParticlePosCount).background = ...
                    max(background);
                app.data.file(fileId).time(tId).particle(trueParticlePosCount).intensity = max(...
                    imageRaw(round(centroids(particleId,2)) - roiRadius : round(centroids(particleId,2)) + roiRadius,...
                    round(centroids(particleId,1)) - roiRadius : round(centroids(particleId,1)) + roiRadius),[],'all');
                app.data.file(fileId).time(tId).particle(trueParticlePosCount).intensity = ...
                    app.data.file(fileId).time(tId).particle(trueParticlePosCount).intensity - ...
                    app.data.file(fileId).time(tId).particle(trueParticlePosCount).background;
                if app.data.file(fileId).time(tId).particle(trueParticlePosCount).intensity <= 0
                    app.data.file(fileId).time(tId).particle(trueParticlePosCount).state = 'rejected';
                end
            end
        end

        % removing particles that do not have maximum intensity within stack
        for particleId = 1 : trueParticlePosCount
            if strcmp(app.data.file(fileId).time(tId).particle(particleId).state,'accepted')
                maxParticleIntensityZId = 1;
                maxParticleIntensity = 0;
                for zId = 1 : size(app.data.file(fileId).image,3)
                    particleIntensityZId = mean(app.data.file(fileId).image(tId,1,zId,...
                        app.data.file(fileId).time(tId).particle(particleId).centroid.y - roiRadius : ...
                        app.data.file(fileId).time(tId).particle(particleId).centroid.y + roiRadius,...
                        app.data.file(fileId).time(tId).particle(particleId).centroid.x - roiRadius : ...
                        app.data.file(fileId).time(tId).particle(particleId).centroid.x + roiRadius),'all');
                    if particleIntensityZId > maxParticleIntensity
                        maxParticleIntensityZId = zId;
                        maxParticleIntensity = particleIntensityZId;
                    end
                end
                if ~(maxParticleIntensityZId > 1 && maxParticleIntensityZId < size(app.data.file(fileId).image,3)) && ...
                        size(app.data.file(fileId).image,3) ~= 1
                    app.data.file(fileId).time(tId).particle(particleId).state = 'rejected';
                end
            end
        end
    end
end
end

%%====================SpecificExport=====================%%
function SpecificExport(app)

% extracting ROI radius
roiRadius = app.param.detection.roiRadius;

% extracting number of files
numFiles = length(app.data.file);

for fileId = 1 : numFiles

    % creating new folder
    mkdir(fullfile(app.param.paths.calibrationAndUnknownData,app.data.file(fileId).type,'detected'));

    for tId = 1 : size(app.data.file(fileId).image,1)
        try

            % extracting number of particles
            numParticles = length(app.data.file(fileId).time(tId).particle);

            % extracting image
            image = uint16(app.data.file(fileId).time(tId).maxImage);

            % calculating ROI intensity
            roiInt = max(image(:));

            % looping through particles
            for particleId = 1 : numParticles

                if strcmp(app.data.file(fileId).time(tId).particle(particleId).state,'accepted')

                    % extracting particle centroid
                    centroid.x = app.data.file(fileId).time(tId).particle(particleId).centroid.x;
                    centroid.y = app.data.file(fileId).time(tId).particle(particleId).centroid.y;

                    % assigning ROI pixels
                    for x = round(centroid.x) - roiRadius : round(centroid.x) + roiRadius
                        for y = round(centroid.y) - roiRadius : round(centroid.y) + roiRadius
                            if (x - round(centroid.x)) ^ 2 + (y - round(centroid.y)) ^ 2 > (roiRadius - 0.5) ^ 2 && ...
                                    (x - round(centroid.x)) ^ 2 + (y - round(centroid.y)) ^ 2 < (roiRadius + 0.5) ^ 2
                                image(y,x) = roiInt;
                            end
                        end
                    end
                end
            end

            % exporting annotated images
            if tId == 1
                imwrite(image,...
                    fullfile(app.param.paths.calibrationAndUnknownData,app.data.file(fileId).type,'detected',[erase(app.data.file(fileId).name,'.lsm') '.tif']),'WriteMode','overwrite');
            else
                imwrite(image,...
                    fullfile(app.param.paths.calibrationAndUnknownData,app.data.file(fileId).type,'detected',[erase(app.data.file(fileId).name,'.lsm') '.tif']),'WriteMode','append');
            end
        catch
        end
    end
end
end