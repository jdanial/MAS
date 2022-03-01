function main(app)
% main entry point
% calls all MAS functions to detect and analyze single punctas of 
% nanoscopic biological assemblies in diffraction-limited microscopy (i.e.
% confocal and widefield).
%
% Copyright 2022 John S H Danial
% Department of Chemistry, Univerity of Cambridge

% setting folder paths - EDITABLE
app.param.paths.calibrationAndUnknownData = app.CalibrationunknowndataEditField.Value;

% detection parameters - EDITABLE
app.param.detection.detect = app.DetectSwitch.Value == 'Y';
app.param.detection.maxSigma = app.MaximumSigmaEditField.Value;
app.param.detection.roiRadius = app.ROIRadiusEditField.Value + 1;

% analysis parameters - EDITABLE
app.param.analysis.analyse = app.AnalyzeSwitch.Value == 'Y';
app.param.analysis.numSubUnitsPerCalibComplex = app.NumSubUnitsPerCalibComplexEditField.Value;
app.param.analysis.timeSlice = app.TimeSliceEditField.Value;
app.param.analysis.binSize = app.BinSizeEditField.Value;

% setting returnFlag
returnFlag = false;

% running modules - DO NOT EDIT
if app.param.detection.detect && ~returnFlag
    returnFlag = detectFn(app);
end
if app.param.analysis.analyse && ~returnFlag
    extractFn(app);
    analyzeFn(app);
end

% clearing RAM
app.data = struct();
end