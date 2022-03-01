function saveParamAuxFn(app)
% saveParamAuxFn() -
% saves BAS parameters.
%
% Syntax -
% exportFn(app).
%
% Parameters -
% - app: BAS UI class

% registering parameters
app.param.detection.detect = app.DetectSwitch.Value;
app.param.detection.maxSigma = app.MaximumSigmaEditField.Value;
app.param.detection.roiRadius = app.ROIRadiusEditField.Value;
app.param.analysis.analyze = app.AnalyzeSwitch.Value;
app.param.analysis.numSubUnitsPerCalibComplex = app.NumSubUnitsPerCalibComplexEditField.Value;
app.param.analysis.timeSlice = app.TimeSliceEditField.Value;
app.param.analysis.binSize = app.BinSizeEditField.Value;

% saving parameters
param = app.param;
exportPath = uigetdir;
save(fullfile(exportPath,'defparam.sp'),'param');
end