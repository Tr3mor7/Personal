function runConfiguration = tape_motion_batch()
%TAPE_MOTION_BATCH Batch-measure the motion of a coloured tape marker in videos.
%
%   tape_motion_batch
%   runs = tape_motion_batch   % also returns each run's settings and status
%
%   For every video, the tape's centroid displacement and its change of
%   angle are measured every N frames and written to Excel.
%
%   GLOBAL SETUP (asked once, used for every run)
%     1. Number of video runs
%     2. RGB tolerance (0-255)
%     3. Units: pixels, or a calibration (pixels per mm, cm, ...)
%     4. Frame interval N between samples (e.g. 30 or 90)
%     5. Output: one combined workbook (one sheet per run) or one workbook per run
%
%   RUN SETUP (asked for each run, all before any analysis starts)
%     - the video
%     - a reference screenshot/photo showing the tape
%     - an output folder (for the marked "_tape" photo, and the workbook in
%       one-workbook-per-run mode)
%     - a rectangle drawn inside the tape on the reference photo; its mean
%       colour is the colour that is tracked
%
%   ANALYSIS then runs unattended. Errors offer retry / reconfigure / skip
%   instead of ending the batch.
%
%   OUTPUT columns (N = frame interval, <unit> = chosen unit)
%     Time_seconds
%     DisplacementMagnitude_EveryNFrames_<unit>
%     AbsoluteAngleChange_EveryNFrames_degrees
%
%   The tracking itself is in track_tape_motion.m.
%
%   Requires MATLAB R2020a or newer (exportgraphics) and the Image Processing Toolbox.

    close all;

    %% Fixed tracking settings (see track_tape_motion for their meaning)
    trackingOptions = struct( ...
        'MaxJumpPixels', 600, ...
        'MinObjectArea', 20, ...
        'ClosingRadius', 3);
    maximumAutomaticAnalysisAttempts = 2;
    maximumFrameReadAttempts = 3;

    %% Global settings
    numberOfRuns = requestPositiveInteger( ...
        'Number of video runs to configure:', 'Batch Tape Analysis', 1);

    trackingOptions.RGBTolerance = requestNumberInRange( ...
        ['RGB tolerance for all runs (0-255):' newline ...
         'Larger values allow more color variation.'], ...
        'Global RGB Tolerance', 35, 0, 255);

    [trackingOptions.PixelsPerUnit, displacementUnit] = requestGlobalCalibration();

    trackingOptions.FrameInterval = requestPositiveInteger( ...
        ['Frames between samples for all runs:' newline ...
         '(e.g. 30 = once per second for a 30 fps video)'], ...
        'Sampling Interval', 30);

    combinedWorkbook = strcmp(questdlg( ...
        ['How should the results be saved?' newline newline ...
         'Combined: one workbook with one sheet per run.' newline ...
         'Per run: one workbook in each run''s output folder.'], ...
        'Excel Output', 'Combined workbook', 'One workbook per run', ...
        'Combined workbook'), 'Combined workbook');

    combinedExcelPath = "";
    if combinedWorkbook
        workbookFolder = selectFolderWithRetry(pwd, ...
            'Select folder for the combined Excel workbook');
        combinedExcelPath = string(fullfile(workbookFolder, ...
            [getSecondLastDirectoryName(workbookFolder), '.xlsx']));
    end

    fprintf('\nGlobal settings:\n');
    fprintf('  Number of runs: %d\n', numberOfRuns);
    fprintf('  RGB tolerance: %.2f\n', trackingOptions.RGBTolerance);
    fprintf('  Displacement unit: %s\n', displacementUnit);
    fprintf('  Frame interval: %d\n', trackingOptions.FrameInterval);
    if combinedWorkbook
        fprintf('  Combined workbook: %s\n\n', combinedExcelPath);
    else
        fprintf('  One workbook per run\n\n');
    end

    %% Configure every run before any analysis
    runConfiguration = repmat(createEmptyRun(), numberOfRuns, 1);
    usedSheetNames = strings(0, 1);

    for runIndex = 1:numberOfRuns
        runConfigured = false;
        while ~runConfigured
            try
                runConfiguration(runIndex) = configureSingleRun( ...
                    runIndex, numberOfRuns, runConfiguration(runIndex));
                runConfigured = true;
            catch setupError
                close all;
                choice = questdlg( ...
                    sprintf(['Run %d setup encountered an error:\n\n%s\n\n' ...
                             'Retry this run, restart it from the beginning, ' ...
                             'or cancel the entire batch?'], ...
                            runIndex, setupError.message), ...
                    sprintf('Run %d Setup Error', runIndex), ...
                    'Retry Run', 'Restart Run', 'Cancel Batch', 'Retry Run');
                if strcmp(choice, 'Cancel Batch') || isempty(choice)
                    error('Batch setup cancelled after run %d error: %s', ...
                        runIndex, setupError.message);
                elseif strcmp(choice, 'Restart Run')
                    runConfiguration(runIndex) = createEmptyRun();
                end
                % "Retry Run" keeps the fields that were already completed.
            end
        end

        sheetName = createUniqueSheetName(runIndex, ...
            char(runConfiguration(runIndex).videoFile), usedSheetNames);
        usedSheetNames(end + 1, 1) = string(sheetName);
        runConfiguration(runIndex).sheetName = string(sheetName);
    end

    close all;
    fprintf('\nAll %d runs are configured. Beginning unattended analysis...\n\n', ...
        numberOfRuns);

    %% Analyse each run, with automatic and user-directed retries
    runResults = cell(numberOfRuns, 1);

    for runIndex = 1:numberOfRuns
        runFinished = false;
        automaticAttempt = 0;

        while ~runFinished
            automaticAttempt = automaticAttempt + 1;
            fprintf('Processing run %d of %d, attempt %d...\n', ...
                runIndex, numberOfRuns, automaticAttempt);

            try
                resultsTable = analyzeConfiguredRun(runConfiguration(runIndex), ...
                    trackingOptions, displacementUnit, maximumFrameReadAttempts);

                if ~combinedWorkbook
                    writeResultsWithRetry(resultsTable, ...
                        char(runConfiguration(runIndex).excelPath), 'Results');
                end

                runResults{runIndex} = resultsTable;
                runConfiguration(runIndex).status = "Completed";
                runConfiguration(runIndex).errorMessage = "";
                runFinished = true;
                fprintf('Completed run %d.\n\n', runIndex);

            catch analysisError
                runConfiguration(runIndex).errorMessage = string(analysisError.message);
                fprintf(2, 'Run %d attempt %d failed: %s\n', ...
                    runIndex, automaticAttempt, analysisError.message);

                if automaticAttempt < maximumAutomaticAnalysisAttempts
                    fprintf('Retrying run %d automatically...\n\n', runIndex);
                    pause(0.5);
                    continue;
                end

                choice = questdlg( ...
                    sprintf(['Run %d could not be completed after %d attempts:\n\n' ...
                             '%s\n\nChoose how to continue.'], ...
                            runIndex, automaticAttempt, analysisError.message), ...
                    sprintf('Run %d Analysis Error', runIndex), ...
                    'Retry Analysis', 'Reconfigure Run', 'Skip Run', 'Retry Analysis');

                if strcmp(choice, 'Retry Analysis')
                    automaticAttempt = 0;
                elseif strcmp(choice, 'Reconfigure Run')
                    [runConfiguration(runIndex), skipped] = reconfigureRun( ...
                        runConfiguration(runIndex), runIndex, numberOfRuns);
                    automaticAttempt = 0;
                    runFinished = skipped;
                else
                    runConfiguration(runIndex).status = "Skipped";
                    runFinished = true;
                    fprintf('Skipped run %d.\n\n', runIndex);
                end
            end
        end
    end

    %% Combined workbook: one sheet per run, including runs that did not finish
    if combinedWorkbook
        writeCombinedWorkbook(char(combinedExcelPath), runConfiguration, runResults, ...
            trackingOptions.FrameInterval, displacementUnit);
    end

    %% Summary
    fprintf('\nBatch analysis finished.\n');
    if combinedWorkbook
        fprintf('Combined workbook: %s (%d sheets)\n', combinedExcelPath, numberOfRuns);
    end
    for runIndex = 1:numberOfRuns
        run = runConfiguration(runIndex);
        fprintf('Run %d | %s | %s\n', runIndex, run.videoFile, run.status);
        if run.status == "Completed"
            if ~combinedWorkbook
                fprintf('  Excel: %s\n', run.excelPath);
            end
            fprintf('  Tape photo: %s\n', run.annotatedPhotoPath);
        else
            fprintf('  Last message: %s\n', run.errorMessage);
        end
    end
end


%% ------------------------------------------------------------------------
%% Run configuration
%% ------------------------------------------------------------------------

function runInfo = createEmptyRun()
    runInfo = struct( ...
        'videoPath', "", 'videoFolder', "", 'videoFile', "", ...
        'referenceImagePath', "", 'photoFile', "", ...
        'outputFolder', "", 'excelPath', "", 'annotatedPhotoPath', "", ...
        'meanReferenceRGB', [NaN, NaN, NaN], ...
        'referenceRectangle', [NaN, NaN, NaN, NaN], ...
        'sheetName', "", 'status', "Not processed", 'errorMessage', "");
end


function [runInfo, skipped] = reconfigureRun(runInfo, runIndex, numberOfRuns)
    skipped = false;
    sheetName = runInfo.sheetName;
    while true
        try
            runInfo = configureSingleRun(runIndex, numberOfRuns, runInfo);
            runInfo.sheetName = sheetName;
            return;
        catch reconfigurationError
            choice = questdlg( ...
                sprintf('Reconfiguration error:\n\n%s\n\nRetry or skip this run?', ...
                        reconfigurationError.message), ...
                'Reconfiguration Error', 'Retry Reconfiguration', 'Skip Run', ...
                'Retry Reconfiguration');
            if strcmp(choice, 'Skip Run') || isempty(choice)
                runInfo.status = "Skipped";
                runInfo.errorMessage = string(reconfigurationError.message);
                skipped = true;
                return;
            end
        end
    end
end


function runInfo = configureSingleRun(runIndex, numberOfRuns, existingRunInfo)
    runInfo = existingRunInfo;
    runLabel = sprintf('Run %d of %d', runIndex, numberOfRuns);
    fprintf('\nConfiguring %s...\n', runLabel);

    %% Video
    videoValid = false;
    while ~videoValid
        if strlength(runInfo.videoPath) == 0 || ~isfile(runInfo.videoPath)
            [videoFile, videoFolder] = selectFileWithRetry( ...
                {'*.mp4;*.avi;*.mov;*.m4v;*.mpg;*.mpeg', 'Video files'; '*.*', 'All files'}, ...
                [runLabel ': Select video'], 'video');
            runInfo.videoFile = string(videoFile);
            runInfo.videoFolder = string(videoFolder);
            runInfo.videoPath = string(fullfile(videoFolder, videoFile));
        end
        try
            videoObject = VideoReader(char(runInfo.videoPath));
            if ~hasFrame(videoObject)
                error('The selected video contains no readable frames.');
            end
            testFrame = readFrame(videoObject);
            if ndims(testFrame) ~= 3 || size(testFrame, 3) ~= 3
                error('The selected video does not contain RGB frames.');
            end
            videoValid = true;
        catch videoError
            choice = questdlg( ...
                sprintf('The selected video could not be opened:\n\n%s\n\nSelect another video or retry this file?', ...
                        videoError.message), ...
                [runLabel ': Video Error'], 'Select Another Video', 'Retry This File', ...
                'Select Another Video');
            if strcmp(choice, 'Select Another Video') || isempty(choice)
                runInfo.videoPath = "";
                runInfo.videoFile = "";
                runInfo.videoFolder = "";
            end
        end
    end
    drawnow;
    pause(0.2);

    %% Reference image
    imageValid = false;
    referenceImage = [];
    while ~imageValid
        if strlength(runInfo.referenceImagePath) == 0 || ~isfile(runInfo.referenceImagePath)
            [photoFile, photoFolder] = selectFileWithRetry( ...
                {'*.png;*.jpg;*.jpeg;*.tif;*.tiff;*.bmp', 'Image files'; '*.*', 'All files'}, ...
                [runLabel ': Select reference screenshot/photo'], 'reference image');
            runInfo.photoFile = string(photoFile);
            runInfo.referenceImagePath = string(fullfile(photoFolder, photoFile));
        end
        try
            referenceImage = imread(char(runInfo.referenceImagePath));
            if ndims(referenceImage) ~= 3 || size(referenceImage, 3) ~= 3
                error('The reference image must be an RGB color image.');
            end
            imageValid = true;
        catch imageError
            choice = questdlg( ...
                sprintf('The reference image could not be read:\n\n%s\n\nSelect another image or retry this file?', ...
                        imageError.message), ...
                [runLabel ': Image Error'], 'Select Another Image', 'Retry This File', ...
                'Select Another Image');
            if strcmp(choice, 'Select Another Image') || isempty(choice)
                runInfo.referenceImagePath = "";
                runInfo.photoFile = "";
            end
        end
    end
    drawnow;
    pause(0.2);

    %% Output folder
    folderValid = false;
    while ~folderValid
        if strlength(runInfo.outputFolder) == 0 || ~isfolder(runInfo.outputFolder)
            startFolder = char(runInfo.videoFolder);
            if isempty(startFolder) || ~isfolder(startFolder)
                startFolder = pwd;
            end
            runInfo.outputFolder = string(selectFolderWithRetry(startFolder, ...
                [runLabel ': Select output directory']));
        end
        try
            assertFolderWritable(char(runInfo.outputFolder));
            folderValid = true;
        catch folderError
            choice = questdlg( ...
                sprintf('The output folder cannot be used:\n\n%s\n\nSelect another folder or retry this folder?', ...
                        folderError.message), ...
                [runLabel ': Output Folder Error'], 'Select Another Folder', ...
                'Retry This Folder', 'Select Another Folder');
            if strcmp(choice, 'Select Another Folder') || isempty(choice)
                runInfo.outputFolder = "";
            end
        end
    end

    %% Tape reference area
    roiValid = false;
    while ~roiValid
        try
            [meanReferenceRGB, referenceRectangle] = ...
                selectReferenceAreaWithRetry(referenceImage, runLabel);
            runInfo.meanReferenceRGB = meanReferenceRGB;
            runInfo.referenceRectangle = referenceRectangle;
            roiValid = true;
        catch roiError
            choice = questdlg( ...
                sprintf('The tape area could not be saved:\n\n%s\n\nRetry tape-area selection or select another image?', ...
                        roiError.message), ...
                [runLabel ': Tape Area Error'], 'Retry Tape Area', 'Select Another Image', ...
                'Retry Tape Area');
            if strcmp(choice, 'Select Another Image')
                runInfo.referenceImagePath = "";
                runInfo.photoFile = "";
                runInfo = configureSingleRun(runIndex, numberOfRuns, runInfo);
                return;
            end
        end
    end

    %% Output names (only once all inputs are valid)
    [~, originalPhotoName, ~] = fileparts(char(runInfo.photoFile));
    runInfo.annotatedPhotoPath = string(fullfile(char(runInfo.outputFolder), ...
        [originalPhotoName, '_tape.png']));
    runInfo.excelPath = string(fullfile(char(runInfo.outputFolder), ...
        [getSecondLastDirectoryName(char(runInfo.outputFolder)), '.xlsx']));

    %% Marked reference photo
    photoSaved = false;
    while ~photoSaved
        try
            saveAnnotatedReferencePhoto(referenceImage, runInfo.referenceRectangle, ...
                char(runInfo.annotatedPhotoPath), runLabel);
            photoSaved = true;
        catch photoError
            choice = questdlg( ...
                sprintf('The marked reference photo could not be saved:\n\n%s\n\nRetry saving or select another output folder?', ...
                        photoError.message), ...
                [runLabel ': Photo Save Error'], 'Retry Save', 'Select Another Folder', ...
                'Retry Save');
            if strcmp(choice, 'Select Another Folder')
                runInfo.outputFolder = "";
                runInfo.excelPath = "";
                runInfo.annotatedPhotoPath = "";
                runInfo = configureSingleRun(runIndex, numberOfRuns, runInfo);
                return;
            end
        end
    end

    runInfo.status = "Configured";
    runInfo.errorMessage = "";
    fprintf('Configured run %d of %d: %s\n', runIndex, numberOfRuns, runInfo.videoFile);
end


%% ------------------------------------------------------------------------
%% Analysis and output
%% ------------------------------------------------------------------------

function resultsTable = analyzeConfiguredRun(runInfo, trackingOptions, ...
        displacementUnit, maximumFrameReadAttempts)

    if strlength(runInfo.videoPath) == 0 || ~isfile(runInfo.videoPath)
        error('The configured video path is missing or invalid.');
    end
    if numel(runInfo.meanReferenceRGB) ~= 3 || any(~isfinite(runInfo.meanReferenceRGB))
        error('The reference RGB values are missing or invalid.');
    end

    videoObject = VideoReader(char(runInfo.videoPath));
    getFrame = @(t) readFrameWithRetry(videoObject, t, maximumFrameReadAttempts);

    results = track_tape_motion(getFrame, videoObject.Duration, videoObject.FrameRate, ...
        runInfo.meanReferenceRGB, trackingOptions);

    resultsTable = table(results.time, results.displacement, results.angleChange, ...
        'VariableNames', resultColumnNames(trackingOptions.FrameInterval, displacementUnit));
end


function names = resultColumnNames(frameInterval, displacementUnit)
    unitText = matlab.lang.makeValidName(char(displacementUnit));
    names = {'Time_seconds', ...
             sprintf('DisplacementMagnitude_Every%dFrames_%s', frameInterval, unitText), ...
             sprintf('AbsoluteAngleChange_Every%dFrames_degrees', frameInterval)};
end


function frame = readFrameWithRetry(videoObject, targetTime, maximumAttempts)
    lastError = [];
    for attempt = 1:maximumAttempts
        try
            videoObject.CurrentTime = max(0, min(targetTime, ...
                max(0, videoObject.Duration - 1e-6)));
            if ~hasFrame(videoObject)
                error('No readable frame exists at this time.');
            end
            frame = readFrame(videoObject);
            if isempty(frame)
                error('MATLAB returned an empty frame.');
            end
            return;
        catch frameError
            lastError = frameError;
            pause(0.15 * attempt);
        end
    end
    rethrow(lastError);
end


function writeResultsWithRetry(resultsTable, excelPath, sheetName)
    while true
        try
            assertFolderWritable(fileparts(excelPath));
            if isfile(excelPath)
                delete(excelPath);
            end
            writetable(resultsTable, excelPath, 'Sheet', sheetName);
            if ~isfile(excelPath)
                error('MATLAB did not create the Excel workbook.');
            end
            return;
        catch writeError
            choice = questdlg( ...
                sprintf('The Excel file could not be written:\n\n%s\n\nClose the workbook if it is open, then retry.', ...
                        writeError.message), ...
                'Excel Write Error', 'Retry Write', 'Stop This Run', 'Retry Write');
            if strcmp(choice, 'Stop This Run') || isempty(choice)
                rethrow(writeError);
            end
        end
    end
end


function writeCombinedWorkbook(workbookPath, runConfiguration, runResults, ...
        frameInterval, displacementUnit)
%WRITECOMBINEDWORKBOOK One sheet per run. Runs that did not finish get an
%   empty sheet with the same columns, so sheet count = run count.
    columns = resultColumnNames(frameInterval, displacementUnit);
    emptyTable = table(zeros(0, 1), zeros(0, 1), zeros(0, 1), 'VariableNames', columns);

    while true
        try
            assertFolderWritable(fileparts(workbookPath));
            if isfile(workbookPath)
                delete(workbookPath);
            end
            for runIndex = 1:numel(runConfiguration)
                resultsTable = runResults{runIndex};
                if isempty(resultsTable)
                    resultsTable = emptyTable;
                end
                writetable(resultsTable, workbookPath, ...
                    'Sheet', char(runConfiguration(runIndex).sheetName), 'Range', 'A1');
            end
            return;
        catch writeError
            choice = questdlg( ...
                sprintf('The combined workbook could not be written:\n\n%s\n\nClose the workbook if it is open, then retry.', ...
                        writeError.message), ...
                'Excel Write Error', 'Retry Write', 'Give Up', 'Retry Write');
            if strcmp(choice, 'Give Up') || isempty(choice)
                rethrow(writeError);
            end
        end
    end
end


function sheetName = createUniqueSheetName(runIndex, videoFile, usedSheetNames)
%CREATEUNIQUESHEETNAME "RunNN_<video name>", valid for Excel (<= 31 chars,
%   no : \ / ? * [ ]) and unique within the workbook.
    [~, videoBaseName, ~] = fileparts(videoFile);
    cleaned = strtrim(regexprep(videoBaseName, '[:\\/?*\[\]]', '_'));
    if isempty(cleaned)
        cleaned = 'Video';
    end
    prefix = sprintf('Run%02d_', runIndex);
    cleaned = cleaned(1:min(end, 31 - numel(prefix)));
    sheetName = [prefix, cleaned];

    duplicateNumber = 2;
    while any(strcmpi(string(sheetName), usedSheetNames))
        suffix = sprintf('_%d', duplicateNumber);
        base = sheetName(1:min(end, 31 - numel(suffix)));
        sheetName = [base, suffix];
        duplicateNumber = duplicateNumber + 1;
    end
end


%% ------------------------------------------------------------------------
%% Dialog helpers
%% ------------------------------------------------------------------------

function value = requestPositiveInteger(promptText, dialogTitle, defaultValue)
    value = NaN;
    while ~isfinite(value) || value < 1 || value ~= round(value)
        answer = inputdlg({promptText}, dialogTitle, [1, 48], {num2str(defaultValue)});
        if isempty(answer)
            requireOrCancel();
            continue;
        end
        value = str2double(answer{1});
        if ~isfinite(value) || value < 1 || value ~= round(value)
            uiwait(errordlg('Enter a positive whole number.', 'Invalid Value', 'modal'));
            value = NaN;
        end
    end
    value = round(value);
end


function value = requestNumberInRange(promptText, dialogTitle, defaultValue, ...
        minimumValue, maximumValue)
    value = NaN;
    while ~isfinite(value) || value < minimumValue || value > maximumValue
        answer = inputdlg({promptText}, dialogTitle, [2, 58], {num2str(defaultValue)});
        if isempty(answer)
            requireOrCancel();
            continue;
        end
        value = str2double(answer{1});
        if ~isfinite(value) || value < minimumValue || value > maximumValue
            uiwait(errordlg(sprintf('Enter a numeric value from %.2f to %.2f.', ...
                minimumValue, maximumValue), 'Invalid Value', 'modal'));
            value = NaN;
        end
    end
end


function requireOrCancel()
    choice = questdlg('This value is required. Retry or cancel the batch?', ...
        'Input Required', 'Retry', 'Cancel Batch', 'Retry');
    if strcmp(choice, 'Cancel Batch') || isempty(choice)
        error('Batch setup was cancelled.');
    end
end


function [pixelsPerUnit, displacementUnit] = requestGlobalCalibration()
    while true
        choice = questdlg(['Choose displacement units for all runs:' newline ...
                           'Use pixels or enter one physical calibration?'], ...
            'Global Displacement Units', 'Pixels', 'Use Calibration', 'Pixels');
        if isempty(choice)
            requireOrCancel();
            continue;
        end
        if strcmp(choice, 'Pixels')
            pixelsPerUnit = 1;
            displacementUnit = "pixels";
            return;
        end

        calibrationInput = inputdlg( ...
            {'Pixels per physical unit for all runs:', 'Physical unit name, such as mm or cm:'}, ...
            'Global Spatial Calibration', [1, 48], {'10', 'mm'});
        if isempty(calibrationInput)
            continue;
        end
        candidatePixelsPerUnit = str2double(calibrationInput{1});
        candidateUnit = string(strtrim(calibrationInput{2}));
        if ~isfinite(candidatePixelsPerUnit) || candidatePixelsPerUnit <= 0
            uiwait(errordlg('Pixels per physical unit must be a positive number.', ...
                'Invalid Calibration', 'modal'));
            continue;
        end
        if strlength(candidateUnit) == 0
            uiwait(errordlg('Enter a physical unit name.', 'Missing Unit', 'modal'));
            continue;
        end
        pixelsPerUnit = candidatePixelsPerUnit;
        displacementUnit = candidateUnit;
        return;
    end
end


function [selectedFile, selectedFolder] = selectFileWithRetry(fileFilter, dialogTitle, fileDescription)
    selectedFile = 0;
    selectedFolder = 0;
    while isequal(selectedFile, 0)
        drawnow;
        pause(0.35);
        try
            [selectedFile, selectedFolder] = uigetfile(fileFilter, dialogTitle, 'MultiSelect', 'off');
        catch selectorError
            selectedFile = 0;
            warning('%s selector error: %s', fileDescription, selectorError.message);
        end
        if isequal(selectedFile, 0)
            choice = questdlg(sprintf(['The %s selection closed or failed.\n' ...
                                       'Retry the selector or cancel the batch?'], fileDescription), ...
                'Selection Not Completed', 'Retry', 'Cancel Batch', 'Retry');
            if strcmp(choice, 'Cancel Batch') || isempty(choice)
                error('Batch cancelled during %s selection.', fileDescription);
            end
        end
    end
end


function selectedFolder = selectFolderWithRetry(startFolder, dialogTitle)
    selectedFolder = 0;
    while isequal(selectedFolder, 0)
        drawnow;
        pause(0.35);
        try
            selectedFolder = uigetdir(startFolder, dialogTitle);
        catch selectorError
            selectedFolder = 0;
            warning('Folder selector error: %s', selectorError.message);
        end
        if isequal(selectedFolder, 0)
            choice = questdlg(['The folder selection closed or failed.' newline ...
                               'Retry the selector or cancel the batch?'], ...
                'Folder Not Selected', 'Retry', 'Cancel Batch', 'Retry');
            if strcmp(choice, 'Cancel Batch') || isempty(choice)
                error('Batch cancelled during folder selection.');
            end
        end
    end
end


function assertFolderWritable(folderPath)
    if ~isfolder(folderPath)
        error('The selected path is not an existing folder.');
    end
    testPath = fullfile(folderPath, sprintf('.tape_write_test_%d.tmp', randi(1e9)));
    fileIdentifier = fopen(testPath, 'w');
    if fileIdentifier == -1
        error('MATLAB could not create a file in this folder.');
    end
    fclose(fileIdentifier);
    delete(testPath);
end


function [meanReferenceRGB, rectangleCoordinates] = selectReferenceAreaWithRetry(referenceImage, runLabel)
%SELECTREFERENCEAREAWITHRETRY Zoom/pan, then draw a rectangle inside the tape.
    while true
        selectionFigure = [];
        try
            selectionFigure = figure('Name', [runLabel ': Select Tape Reference'], ...
                'NumberTitle', 'off', 'Toolbar', 'figure');
            selectionAxes = axes(selectionFigure);
            imshow(referenceImage, 'Parent', selectionAxes);
            title(selectionAxes, {[runLabel ': Step 1 - Zoom or pan to the tape.'], ...
                                  'Click "Continue to ROI Selection" when ready.'});
            try
                axtoolbar(selectionAxes, {'zoomin', 'zoomout', 'pan', 'restoreview'});
            catch
                % the standard figure toolbar is still available
            end
            zoom(selectionFigure, 'on');
            continueButton = uicontrol(selectionFigure, 'Style', 'pushbutton', ...
                'String', 'Continue to ROI Selection', 'Units', 'normalized', ...
                'Position', [0.36, 0.015, 0.28, 0.055], 'FontWeight', 'bold', ...
                'Callback', @(~, ~) uiresume(selectionFigure));
            uiwait(selectionFigure);
            if ~isgraphics(selectionFigure)
                error('The reference-selection window was closed.');
            end
            zoom(selectionFigure, 'off');
            pan(selectionFigure, 'off');
            if isgraphics(continueButton)
                delete(continueButton);
            end

            title(selectionAxes, {[runLabel ': Step 2 - Draw inside the tape.'], ...
                                  'Double-click the rectangle or press Enter when finished.'});
            referenceROI = drawrectangle(selectionAxes, 'Color', 'y', 'LineWidth', 1.5);
            wait(referenceROI);
            if ~isvalid(referenceROI)
                error('The tape rectangle was not completed.');
            end

            roiPosition = round(referenceROI.Position);
            x1 = max(1, roiPosition(1));
            y1 = max(1, roiPosition(2));
            x2 = min(size(referenceImage, 2), x1 + max(1, roiPosition(3)) - 1);
            y2 = min(size(referenceImage, 1), y1 + max(1, roiPosition(4)) - 1);
            if x2 <= x1 || y2 <= y1
                error('The selected tape area is too small.');
            end

            referencePixels = reshape(double(referenceImage(y1:y2, x1:x2, :)), [], 3);
            candidateMeanRGB = mean(referencePixels, 1);
            if any(~isfinite(candidateMeanRGB))
                error('The selected area produced invalid RGB values.');
            end

            rectangleCoordinates = [x1, y1, x2, y2];
            meanReferenceRGB = candidateMeanRGB;
            close(selectionFigure);
            return;

        catch selectionError
            if ~isempty(selectionFigure) && isgraphics(selectionFigure)
                close(selectionFigure);
            end
            choice = questdlg(sprintf('Tape-area selection was not completed:\n\n%s\n\nRetry the tape-area selection?', ...
                                      selectionError.message), ...
                [runLabel ': Retry Tape Area'], 'Retry', 'Return to Run Setup', 'Retry');
            if strcmp(choice, 'Return to Run Setup') || isempty(choice)
                rethrow(selectionError);
            end
        end
    end
end


function saveAnnotatedReferencePhoto(referenceImage, rectangleCoordinates, outputPath, runLabel)
    annotationFigure = figure('Visible', 'off', 'Color', 'w', 'Position', [100, 100, 1000, 700]);
    cleanup = onCleanup(@() close(annotationFigure));
    annotationAxes = axes(annotationFigure);
    imshow(referenceImage, 'Parent', annotationAxes);
    hold(annotationAxes, 'on');
    r = rectangleCoordinates;
    rectangle(annotationAxes, 'Position', [r(1), r(2), r(3) - r(1) + 1, r(4) - r(2) + 1], ...
        'EdgeColor', 'y', 'LineWidth', 3);
    title(annotationAxes, [runLabel ': Selected Tape Reference Area'], 'Interpreter', 'none');
    exportgraphics(annotationAxes, outputPath, 'Resolution', 200);
    if ~isfile(outputPath)
        error('MATLAB did not create the marked reference photo.');
    end
end


function directoryName = getSecondLastDirectoryName(folderPath)
%GETSECONDLASTDIRECTORYNAME Name of the parent of folderPath, used to name the
%   workbook (e.g. .../Trial_03/Output -> Trial_03.xlsx).
    folderPath = regexprep(char(folderPath), '[\\/]+$', '');
    [parentFolder, lastDirectory] = fileparts(folderPath);
    [~, secondLastDirectory] = fileparts(parentFolder);
    if isempty(secondLastDirectory)
        secondLastDirectory = lastDirectory;
    end
    directoryName = regexprep(secondLastDirectory, '[<>:"/\\|?*]', '_');
    if isempty(directoryName)
        directoryName = 'Tape_Motion_Results';
    end
end
