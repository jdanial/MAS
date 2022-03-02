function particleLocalizerAuxFn(app,fileId)
% particleLocalizerAuxFn - (Auxillary function)
% localizes particles.
%
% Syntax -
% particleLocalizerAuxFn(app,fileId,particleId)
%
% Parameters -
% - app: MAS UI class.
% - fileId: file #.
% - particleId: particle #.

%% setting up particle image grid for bkgd calculation
roiRadius = app.param.detection.roiRadius;
roiWidth = roiRadius * 2 + 1;

%% looping through time
for tId = 1 : size(app.data.file(fileId).image,1)
    try

        %% extracting image
        image = double(squeeze(app.data.file(fileId).image(tId,1,app.data.file(fileId).maxZId.time(tId),:,:)));

        %% looping across accepted particles
        indices = [];
        centroids = [];
        numParticles = length(app.data.file(fileId).time(tId).particle);
        for particleId = 1 : numParticles
            if strcmp(app.data.file(fileId).time(tId).particle(particleId).state,'accepted')
                indices = [indices particleId];
                centroids(particleId,:) = [app.data.file(fileId).time(tId).particle(particleId).centroid.y ...
                    app.data.file(fileId).time(tId).particle(particleId).centroid.x];
            end
        end

        %% finding centroids
        subImage = [];
        fitData = single([]);
        fitInitParams = single([]);
        for particleId = 1 : length(indices)
            subImage(:,:,particleId) = image(round(centroids(indices(particleId),1)) + roiRadius : -1 :...
                round(centroids(indices(particleId),1)) - roiRadius,...
                round(centroids(indices(particleId),2)) - roiRadius:...
                round(centroids(indices(particleId),2)) + roiRadius);
        end
        maxSubImage = max(subImage(:,:,1 : length(indices)),[],[1,2]);
        minSubImage = min(subImage(:,:,1 : length(indices)),[],[1,2]);
        fitData(:,1 : length(indices)) = ...
            single(reshape(subImage(:,:,1 : length(indices)),[roiWidth * roiWidth,length(indices)]));
        fitInitParams(:,1 : length(indices)) = ...
            single([reshape(maxSubImage(1,1,1:length(indices)),1,[]) - reshape(minSubImage(1,1,1:length(indices)),1,[]) ; ...
            4 * ones(1, length(indices)); ...
            4 * ones(1, length(indices)) ; ...
            1 * ones(1, length(indices)) ; ...
            reshape(minSubImage(1,1,1:length(indices)),1,[])]);

        %% fitting particles
        [parameters, states, chX, ~, ~] = ...
            cpufit(fitData, ...
            [], ...
            ModelID.GAUSS_2D, ...
            fitInitParams, ...
            [], ...
            [], ...
            [], ...
            []);

        %% transferring data to particleData struct

        for particleId = 1 : size(fitData,2)
            if states(particleId) == 0

                % updating sigma
                app.data.file(fileId).time(tId).particle(indices(particleId)).sigma = parameters(4,particleId);
                app.data.file(fileId).time(tId).particle(indices(particleId)).intensity = ...
                    (parameters(1,particleId) .* parameters(4,particleId) .* sqrt(2 * pi)) - ...
                    parameters(5,particleId);
            else
                app.data.file(fileId).time(tId).particle(indices(particleId)).state = 'rejected';
            end
        end

        %% rejecting particles with chisquared value
        medianChX = median(chX);
        for particleId = 1 : size(fitData,2)
            if strcmp(app.data.file(fileId).time(tId).particle(indices(particleId)).state,'accepted')
                if chX(particleId) > medianChX || app.data.file(fileId).time(tId).particle(indices(particleId)).intensity < 0 || parameters(4,particleId) < 0 || parameters(1,particleId) > max(image,[],'all')
                    app.data.file(fileId).time(tId).particle(indices(particleId)).state = 'rejected';
                end
            end
        end
    catch
    end
end
end