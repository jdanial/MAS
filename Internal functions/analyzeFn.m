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
        tId = 1;

        % extracting number of particles
        numParticles = length(app.data.file(fileId).time(tId).particle);

        % looping over particles
        for particleId = 1 : numParticles

            % checking if particle is accepted
            if strcmp(app.data.file(fileId).time(tId).particle(particleId).state,'accepted')

                % calculating the intensity
                intCalibration = [intCalibration ...
                    double(app.data.file(fileId).time(tId).particle(particleId).intensity)];
            end
        end
    else

        % looping through time
        for tId = 1 : app.param.detection.timeLength
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
                intUnknown.file(fileId).time(tId).intensity = [intUnknown.file(fileId).time(tId).intensity double(intUnknownTemp)];
            catch
                intUnknown.file(fileId).time(tId).intensity = double(intUnknownTemp);
            end
        end
    end
end

% calculating the mean and standard deviation of the calibration molecularity
app.data.intCalibration = intCalibration;
probDist = fitdist(app.data.intCalibration','Normal');
app.data.meanIntCalibration = probDist.mu;

% calculating molecularity
for fileId = 1 : numFiles
    if strcmp(app.data.file(fileId).type,'Unknown')
        for tId = 1 : app.param.detection.timeLength
            app.data.file(fileId).molecularity.time(tId).index = round((intUnknown.file(fileId).time(tId).intensity ./ ...
                app.data.meanIntCalibration) .* ...
                numSubUnitsPerCalibComplex);
        end
    end
end

% calculating molecularity kernel
for fileId = 1 : numFiles
    if strcmp(app.data.file(fileId).type,'Unknown')
        for tId = 1 : app.param.detection.timeLength
            try
                app.data.molecularity.time(tId).index = [app.data.molecularity.time(tId).index app.data.file(fileId).molecularity.time(tId).index];
            catch
                app.data.molecularity.time(tId).index = app.data.file(fileId).molecularity.time(tId).index;
            end
        end
    end
end
for tId = 1 : app.param.detection.timeLength
    x_unknown = min(app.data.molecularity.time(tId).index(:)) : 1 : max(app.data.molecularity.time(tId).index(:));
    pd = fitdist(app.data.molecularity.time(tId).index','Kernel','BandWidth',binSize);
    y_unknown = pdf(pd,x_unknown);
    app.data.curve.time(tId).y_unknown = y_unknown;
    app.data.curve.time(tId).x_unknown = x_unknown;
end

% calculating average and std molecularity 
for fileId = 1 : numFiles
    if strcmp(app.data.file(fileId).type,'Unknown')
        for tId = 1 : app.param.detection.timeLength
            [weights,vals] = hist(app.data.file(fileId).molecularity.time(tId).index,unique(app.data.file(fileId).molecularity.time(tId).index));
            [freqs,freqVals] = ecdf(app.data.file(fileId).molecularity.time(tId).index);
            fitFunc = fittype('-exp(-x/a) + 1');
            fitParam = fit(freqVals,freqs,fitFunc,'StartPoint',100);
            app.data.file(fileId).molecularity.time(tId).mean = mean(app.data.file(fileId).molecularity.time(tId).index);
            app.data.file(fileId).molecularity.time(tId).median = median(app.data.file(fileId).molecularity.time(tId).index);
            app.data.file(fileId).molecularity.time(tId).wMean = sum((weights ./ sum(weights)) .* vals);
            app.data.file(fileId).molecularity.time(tId).eMean = fitParam.a;
            try
                app.data.molecularity.time(tId).meanTemp = [app.data.molecularity.time(tId).meanTemp app.data.file(fileId).molecularity.time(tId).mean];
                app.data.molecularity.time(tId).medianTemp = [app.data.molecularity.time(tId).medianTemp app.data.file(fileId).molecularity.time(tId).median];
                app.data.molecularity.time(tId).wMeanTemp = [app.data.molecularity.time(tId).wMeanTemp app.data.file(fileId).molecularity.time(tId).wMean];
                app.data.molecularity.time(tId).eMeanTemp = [app.data.molecularity.time(tId).eMeanTemp app.data.file(fileId).molecularity.time(tId).eMean];
            catch
                app.data.molecularity.time(tId).meanTemp = app.data.file(fileId).molecularity.time(tId).mean;
                app.data.molecularity.time(tId).medianTemp = app.data.file(fileId).molecularity.time(tId).median;
                app.data.molecularity.time(tId).wMeanTemp = app.data.file(fileId).molecularity.time(tId).wMean;
                app.data.molecularity.time(tId).eMeanTemp = app.data.file(fileId).molecularity.time(tId).eMean;
            end
        end
    end
end
for tId = 1 : app.param.detection.timeLength
    app.data.molecularity.time(tId).mean = mean(app.data.molecularity.time(tId).meanTemp);
    app.data.molecularity.time(tId).median = mean(app.data.molecularity.time(tId).medianTemp);
    app.data.molecularity.time(tId).wMean = mean(app.data.molecularity.time(tId).wMeanTemp);
    app.data.molecularity.time(tId).eMean = mean(app.data.molecularity.time(tId).eMeanTemp);
    app.data.molecularity.time(tId).sd = std(app.data.molecularity.time(tId).meanTemp);
    app.data.molecularity.time(tId).wSd = std(app.data.molecularity.time(tId).wMeanTemp);
    app.data.molecularity.time(tId).eSd = std(app.data.molecularity.time(tId).eMeanTemp);
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
for tId = 1 : app.param.detection.timeLength

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
    'Unknown_molecularity_mean.txt'),'w');

% writing values
fprintf(fileHandle,'%s\t%s\t%s\n','Time (min)','Mean molecularity','Std molecularity');
for tId = 1 : app.param.detection.timeLength
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
for tId = 2 : app.param.detection.timeLength
    if ~isnan(app.data.molecularity.time(tId).mean)
        y = [y app.data.molecularity.time(tId).mean];
        yConf = [yConf app.data.molecularity.time(tId).mean + app.data.molecularity.time(tId).sd];
        x = [x (tId - 1) * app.param.analysis.timeSlice];
        xConf = [xConf (tId - 1) * app.param.analysis.timeSlice];
    end
end
for tId = app.param.detection.timeLength : - 1 : 2 
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
ylabel('Mean molecularity');
pbaspect([1 1 1]);
box on;

% saving figure
saveas(figHandle,fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Unknown_molecularity_mean.png'));



% calculating file handle and opening file
fileHandle = fopen(fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Unknown_molecularity_median.txt'),'w');

% writing values
fprintf(fileHandle,'%s\t%s\t%s\n','Time (min)','Median molecularity','Std molecularity');
for tId = 1 : app.param.detection.timeLength
    try
        fprintf(fileHandle,'%d\t%d\t%d\n',(tId - 1) * app.param.analysis.timeSlice,...
            app.data.molecularity.time(tId).median,...
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
for tId = 2 : app.param.detection.timeLength
    if ~isnan(app.data.molecularity.time(tId).median)
        y = [y app.data.molecularity.time(tId).median];
        yConf = [yConf app.data.molecularity.time(tId).median + app.data.molecularity.time(tId).sd];
        x = [x (tId - 1) * app.param.analysis.timeSlice];
        xConf = [xConf (tId - 1) * app.param.analysis.timeSlice];
    end
end
for tId = app.param.detection.timeLength : - 1 : 2 
    if ~isnan(app.data.molecularity.time(tId).median)
        yConf = [yConf app.data.molecularity.time(tId).median - app.data.molecularity.time(tId).sd];
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
ylabel('Median molecularity');
pbaspect([1 1 1]);
box on;

% saving figure
saveas(figHandle,fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Unknown_molecularity_median.png'));



% calculating file handle and opening file
fileHandle = fopen(fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Unknown_molecularity_weighted_mean.txt'),'w');

% writing values
fprintf(fileHandle,'%s\t%s\t%s\n','Time (min)','Weighted mean molecularity','Std molecularity');
for tId = 1 : app.param.detection.timeLength
    try
        fprintf(fileHandle,'%d\t%d\t%d\n',(tId - 1) * app.param.analysis.timeSlice,...
            app.data.molecularity.time(tId).wMedian,...
            app.data.molecularity.time(tId).wSd);
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
for tId = 2 : app.param.detection.timeLength
    if ~isnan(app.data.molecularity.time(tId).wMean)
        y = [y app.data.molecularity.time(tId).wMean];
        yConf = [yConf app.data.molecularity.time(tId).wMean + app.data.molecularity.time(tId).wSd];
        x = [x (tId - 1) * app.param.analysis.timeSlice];
        xConf = [xConf (tId - 1) * app.param.analysis.timeSlice];
    end
end
for tId = app.param.detection.timeLength : - 1 : 2 
    if ~isnan(app.data.molecularity.time(tId).wMean)
        yConf = [yConf app.data.molecularity.time(tId).wMean - app.data.molecularity.time(tId).wSd];
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
ylabel('Weighted mean molecularity');
pbaspect([1 1 1]);
box on;

% saving figure
saveas(figHandle,fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Unknown_molecularity_weighted_mean.png'));



% calculating file handle and opening file
fileHandle = fopen(fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Unknown_molecularity_emperical_mean.txt'),'w');

% writing values
fprintf(fileHandle,'%s\t%s\t%s\n','Time (min)','Emperical mean molecularity','Std molecularity');
for tId = 1 : app.param.detection.timeLength
    try
        fprintf(fileHandle,'%d\t%d\t%d\n',(tId - 1) * app.param.analysis.timeSlice,...
            app.data.molecularity.time(tId).eMedian,...
            app.data.molecularity.time(tId).eSd);
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
for tId = 2 : app.param.detection.timeLength
    if ~isnan(app.data.molecularity.time(tId).eMean)
        y = [y app.data.molecularity.time(tId).eMean];
        yConf = [yConf app.data.molecularity.time(tId).eMean + app.data.molecularity.time(tId).eSd];
        x = [x (tId - 1) * app.param.analysis.timeSlice];
        xConf = [xConf (tId - 1) * app.param.analysis.timeSlice];
    end
end
for tId = app.param.detection.timeLength : - 1 : 2 
    if ~isnan(app.data.molecularity.time(tId).eMean)
        yConf = [yConf app.data.molecularity.time(tId).eMean - app.data.molecularity.time(tId).eSd];
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
ylabel('Emperical mean molecularity');
pbaspect([1 1 1]);
box on;

% saving figure
saveas(figHandle,fullfile(app.param.paths.calibrationAndUnknownData,...
    'analysis',...
    'Unknown_molecularity_emperical_mean.png'));



% looping through time
for tId = 1 : app.param.detection.timeLength
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
        for tId = 1 : app.param.detection.timeLength

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
        for tId = 1 : app.param.detection.timeLength

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