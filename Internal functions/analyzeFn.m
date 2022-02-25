function returnFlag = analyzeFn(app)
% analyzeFn() -
% analyzes data processed by BAS.
%
% Syntax -
% analyzeFn(app).
%
% Parameters -
% - app: MAS UI class

%% initializing returnFlag
returnFlag = false;

%% displaying SAS progress
app.msgBox.Value = sprintf('%s','Data analysis started.');
drawnow;

%% calling ParticleDetector
app.msgBox.Value = sprintf('%s','Calculating molecularity.');
drawnow;
MolecularityCalculator(app);

%% exporting files relevant to this function
app.msgBox.Value = sprintf('%s','Exporting analysis data.');
drawnow;
SpecificExport(app);

app.msgBox.Value = sprintf('%s','Done.');
drawnow;
end

%%====================KernelGenerator=====================%%
function MolecularityCalculator(app)

% extracting number of subunits per complex used for calibration
numSubUnitsPerCalibComplex = app.param.analysis.numSubUnitsPerCalibComplex;
binSize = app.param.analysis.binSize;

% extracting number of files
numFiles = length(app.data.file);

% initializing arrays
intCalibration = [];
intUnknown = struct();

% looping over files
for fileId = 1 : numFiles

    % checking if file is a calibration file
    if strcmp(app.data.file(fileId).type,'Calibration')

        % looping through time
        for tId = 1 : size(app.data.file(fileId).image,1)

            % extracting number of particles
            numParticles = length(app.data.file(fileId).time(tId).particle);

            % looping over particles
            for particleId = 1 : numParticles

                % checking if particle is accepted
                if strcmp(app.data.file(fileId).time(tId).particle(particleId).state,'accepted')

                    % calculating the intensity
                    intCalibration = [intCalibration ...
                        app.data.file(fileId).time(tId).particle(particleId).intensity];
                end
            end
        end
    else

        % looping through time
        for tId = 1 : length(app.data.file(fileId).time)
            intUnknownTemp = [];
            
            % extracting number of particles
            numParticles = length(app.data.file(fileId).time(tId).particle);

            % looping over particles
            for particleId = 1 : numParticles

                % checking if particle is accepted
                if strcmp(app.data.file(fileId).time(tId).particle(particleId).state,'accepted')

                    % calculating the intensity
                    intUnknownTemp = [intUnknownTemp ...
                        app.data.file(fileId).time(tId).particle(particleId).intensity];
                end
            end
            try
                intUnknown.time(tId).intensity = [intUnknown.time(tId).intensity intUnknownTemp];
            catch
                intUnknown.time(tId).intensity = intUnknownTemp;
            end
        end
    end
end

% calculating the mean and standard deviation of the calibration molecularity
app.data.intCalibration = intCalibration;
app.data.meanIntCalibration = median(app.data.intCalibration);
app.data.stdIntCalibration = std(app.data.intCalibration);

% calculating molecularity
for tId = 1 : length(intUnknown.time)
    app.data.molecularity.time(tId).index = round((intUnknown.time(tId).intensity ./ ...
        app.data.meanIntCalibration) .* ...
        numSubUnitsPerCalibComplex);
end

% calculating molecularity count
for tId = 1 : length(intUnknown.time)
    for moleculeId = 1 : max(app.data.molecularity.time(tId).index)
        app.data.molecularity.count.time(tId).molecule(moleculeId) = 0;
        for indexId = 1 : length(app.data.molecularity.time(tId).index)
            if app.data.molecularity.time(tId).index(indexId) == moleculeId
                app.data.molecularity.count.time(tId).molecule(moleculeId) = ...
                    app.data.molecularity.count.time(tId).molecule(moleculeId) + 1;
            end
        end
    end
end

% calculating molecularity kernel
for tId = 1 : length(intUnknown.time)
    try
        x_unknown = min(app.data.molecularity.time(tId).index(:)) : 1 : max(app.data.molecularity.time(tId).index(:));
        pd = fitdist(intUnknown.time(tId).intensity','Kernel','BandWidth',binSize);
        y_unknown = pdf(pd,x_unknown);
        app.data.curve.time(tId).y_unknown = y_unknown;
        app.data.curve.time(tId).x_unknown = x_unknown;
    catch
    end
end

% calculating average molecularily 
for tId = 1 : length(intUnknown.time)
    try
        app.data.molecularity.time(tId).mean = mean(app.data.molecularity.time(tId).index);
        app.data.molecularity.time(tId).sd = std(app.data.molecularity.time(tId).index);
    catch
    end
end
end

%%====================SpecificExport=====================%%
function SpecificExport(app)

% creating new folder
mkdir(fullfile(app.param.paths.calibrationAndUnknownData,'analysis'));

% calculating file handle and opening file
fileHandle = fopen(fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Calibration_intensity.txt'),'w');

% writing values
fprintf(fileHandle,'%s\n','Intensity');
for unitId = 1 : length(app.data.intCalibration)
    fprintf(fileHandle,'%d\t%d\n',unitId,app.data.intCalibration(unitId));
end

% closing file
fclose(fileHandle);

% creating figure (box plot of calibration intensities)
figHandle = figure('visible','off');

% drawing plot
boxplot(app.data.intCalibration);
xlabel('Calibration set');
ylabel('Intensity (A.U.)');

% saving figure
saveas(figHandle,fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Calibration_boxplot.png'));

% creating figure (PDF of molecularity)
figHandle = figure('visible','off');

% looping through time
for tId = 1 : length(app.data.molecularity.count.time)

    try

        % drawing kernel density estimator
        plot(app.data.curve.time(tId).x_unknown,app.data.curve.time(tId).y_unknown,...
            'DisplayName',[num2str(app.param.analysis.timeSlice * (tId - 1)) ' min']);
        hold on;
    catch
    end

end
hold off;
legend
xlabel('Molecularity');
ylabel('Frequency');

% saving figure
saveas(figHandle,fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Unknown_molecularity_PDF_combined.png'));

% creating figure (PDF of average molecularity)
figHandle = figure('visible','off');

% calculating file handle and opening file
fileHandle = fopen(fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Unknown_molecularity_average.txt'),'w');

% writing values
fprintf(fileHandle,'%s\t%s\t%s\n','Time (min)','Mean molecularity','Std molecularity');
for tId = 1 : length(app.data.molecularity.count.time)
    try
        fprintf(fileHandle,'%d\t%d\t%d\n',(tId - 1) * app.param.analysis.timeSlice,...
            app.data.molecularity.time(tId).mean,...
            app.data.molecularity.time(tId).sd);
    catch
    end
end

% closing file
fclose(fileHandle);

% looping through time
xConf = [];
yConf = [];
x = [];
y = [];
for tId = 1 : length(app.data.molecularity.count.time)
    if ~isnan(app.data.molecularity.time(tId).mean)
        y = [y app.data.molecularity.time(tId).mean];
        yConf = [yConf app.data.molecularity.time(tId).mean + app.data.molecularity.time(tId).sd];
        x = [x (tId - 1) * app.param.analysis.timeSlice];
        xConf = [xConf (tId - 1) * app.param.analysis.timeSlice];
    end
end
for tId = length(app.data.molecularity.count.time) : - 1 : 1 
    if ~isnan(app.data.molecularity.time(tId).mean)
        yConf = [yConf app.data.molecularity.time(tId).mean - app.data.molecularity.time(tId).sd];
        xConf = [xConf (tId - 1) * app.param.analysis.timeSlice];
    end
end
sdInt = fill(xConf,yConf,'red');
sdInt.FaceColor = [1 0.8 0.8];      
sdInt.EdgeColor = 'none';
hold on;
plot(x,y,'r-');
hold off;
xlabel('Time (min)');
ylabel('Average molecularity');

% saving figure
saveas(figHandle,fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Unknown_molecularity_average.png'));

% looping through time
for tId = 1 : length(app.data.molecularity.count.time)
    try

        % calculating file handle and opening file
        fileHandle = fopen(fullfile(app.param.paths.calibrationAndUnknownData,...
            'analysis',...
            ['Unknown_molecularity_' num2str(app.param.analysis.timeSlice * (tId - 1)) 'min.txt']),'w');

        % writing proportions
        fprintf(fileHandle,'%s\n','Molecularity');
        for unitId = 1 : length(app.data.molecularity.time(tId).index)
            fprintf(fileHandle,'%d\t%d\n',unitId,app.data.molecularity.time(tId).index(unitId));
        end

        % closing file
        fclose(fileHandle);

        % calculating file handle and opening file
        fileHandle = fopen(fullfile(app.param.paths.calibrationAndUnknownData,...
            'analysis',...
            ['Unknown_molecularity_count_' num2str(app.param.analysis.timeSlice * (tId - 1)) 'min.txt']),'w');

        % writing proportions
        fprintf(fileHandle,'%s\n','Molecularity');
        for molecularityId = 1 : length(app.data.molecularity.count.time(tId).molecule)
            fprintf(fileHandle,'%d\t%d\n',molecularityId,app.data.molecularity.count.time(tId).molecule(molecularityId));
        end

        % closing file
        fclose(fileHandle);

        % saving PDF of unknown data in text file
        % calculating file handle and opening file
        fileHandle = fopen(fullfile(app.param.paths.calibrationAndUnknownData,...
            'analysis',...
            ['Unknown_molecularity_PDF_' num2str(app.param.analysis.timeSlice * (tId - 1)) 'min.txt']),'w');

        % creating figure (PDF of molecularity)
        figHandle = figure('visible','off');

        % drawing kernel density estimator
        plot(app.data.curve.time(tId).x_unknown,app.data.curve.time(tId).y_unknown);

        % recording kernel density estimator in text file
        fprintf(fileHandle,'%s\n','Kernel Density Estimator');
        fprintf(fileHandle,'%s\t%s\n','Molecularity','Frequency');
        xCurve = app.data.curve.time(tId).x_unknown;
        yCurve = app.data.curve.time(tId).y_unknown;
        for indexId = 1 : length(app.data.curve.time(tId).x_unknown)
            fprintf(fileHandle,'%f\t%f\n',...
                xCurve(indexId),...
                yCurve(indexId));
        end
        xlabel('Molecularity');
        ylabel('Frequency')

        % saving figure
        saveas(figHandle,fullfile(app.param.paths.calibrationAndUnknownData,...
            'analysis',...
            ['Unknown_molecularity_PDF_' num2str(app.param.analysis.timeSlice * (tId - 1)) 'min.png']));

        % closing file
        fclose(fileHandle);
    catch
    end
end

% saving coordinates and intensities of calibration and unknown data in text file
% calculating file handle and opening file
fileHandle = fopen(fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Calibration_raw.txt'),'w');
fprintf(fileHandle,'%s\t%s\t%s\t%s\n','Intensity','x','y','t');
numFiles = length(app.data.file);
for fileId = 1 : numFiles

    % checking if file is a calibration file
    if strcmp(app.data.file(fileId).type,'Calibration')

        % looping through time
        for tId = 1 : length(app.data.file(fileId).time)

            try

                % extracting number of particles
                numParticles = length(app.data.file(fileId).time(tId).particle);

                % looping over particles
                for particleId = 1 : numParticles

                    % checking if particle is accepted
                    if strcmp(app.data.file(fileId).time(tId).particle(particleId).state,'accepted')

                        if app.data.file(fileId).time(tId).particle(particleId).monomeric

                            % printing to file
                            fprintf(fileHandle,'%f\t%f\t%f\n',...
                                app.data.file(fileId).time(tId).particle(particleId).intensity,...
                                app.data.file(fileId).time(tId).particle(particleId).centroid.x,...
                                app.data.file(fileId).time(tId).particle(particleId).centroid.y);
                        end
                    end
                end
            catch
            end
        end
    end
end

% closing file
fclose(fileHandle);

fileHandle = fopen(fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Unknown_raw.txt'),'w');
fprintf(fileHandle,'%s\t%s\t%s\t%s\n','Intensity','x','y','t');
numFiles = length(app.data.file);
for fileId = 1 : numFiles

    % checking if file is an unknown file
    if strcmp(app.data.file(fileId).type,'Unknown')

        % looping through time
        for tId = 1 : length(app.data.file(fileId).time)

            try

                % extracting number of particles
                numParticles = length(app.data.file(fileId).time(tId).particle);

                % looping over particles
                for particleId = 1 : numParticles

                    % checking if particle is accepted
                    if strcmp(app.data.file(fileId).time(tId).particle(particleId).state,'accepted')

                        % printing to file
                        fprintf(fileHandle,'%f\t%f\t%f\n',...
                            app.data.file(fileId).time(tId).particle(particleId).intensity,...
                            app.data.file(fileId).time(tId).particle(particleId).centroid.x,...
                            app.data.file(fileId).time(tId).particle(particleId).centroid.y);
                    end
                end
            catch
            end
        end
    end
end

% closing file
fclose(fileHandle);

% saving all data
data = app.data;
save(fullfile(app.param.paths.calibrationAndUnknownData,'All data.mat'),'data');
end