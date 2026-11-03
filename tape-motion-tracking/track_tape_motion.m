function results = track_tape_motion(getFrame, duration, frameRate, referenceRGB, options)
%TRACK_TAPE_MOTION Track a coloured tape marker through a video.
%
%   results = track_tape_motion(getFrame, duration, frameRate, referenceRGB, options)
%
%   This is the analysis core of tape_motion_batch. It has no dialogs and no
%   file output, so it can be scripted and tested on its own.
%
%   Inputs
%     getFrame      function handle, frame = getFrame(timeSeconds), returning
%                   an RGB frame (H x W x 3). It may return [] or throw an
%                   error when a frame cannot be read; that sample is then
%                   recorded as missing (NaN).
%     duration      video duration in seconds
%     frameRate     frames per second
%     referenceRGB  1x3 mean colour of the tape (0-255), e.g. from a
%                   rectangle drawn inside the tape on a reference image
%     options       struct; any missing field takes its default:
%       FrameInterval   frames between samples                    (30)
%       RGBTolerance    max per-channel difference, 0-255          (35)
%       PixelsPerUnit   spatial calibration; 1 = report pixels     (1)
%       MaxJumpPixels   max centroid jump between samples          (600)
%       MinObjectArea   smallest blob kept, pixels                 (20)
%       ClosingRadius   morphological closing radius, pixels       (3)
%
%   Output struct (column vectors, one row per sample after the first frame)
%     time            sample time, seconds
%     displacement    centroid displacement since the previous sample,
%                     in calibrated units (NaN if either sample is missing)
%     angleChange     absolute change of the tape's long-axis angle since
%                     the previous sample, degrees, 0-90 (NaN if missing)
%     detected        true where the tape was found in that sample
%
%   The tape is the connected blob whose colour is within RGBTolerance of
%   referenceRGB on every channel. In the first frame the largest such blob
%   is taken; afterwards, the blob closest to the previous position whose
%   area is within a factor of 5 of the previous area.
%
%   Requires the Image Processing Toolbox (MATLAB) or the image package (Octave).

    if nargin < 5
        options = struct();
    end
    options = apply_defaults(options);

    if ~isfinite(frameRate) || frameRate <= 0
        error('track_tape_motion:frameRate', 'The video frame rate is invalid.');
    end
    if numel(referenceRGB) ~= 3 || any(~isfinite(referenceRGB))
        error('track_tape_motion:referenceRGB', 'referenceRGB must be three finite values.');
    end

    firstFrame = getFrame(0);
    if isempty(firstFrame) || ndims(firstFrame) ~= 3 || size(firstFrame, 3) ~= 3
        error('track_tape_motion:firstFrame', 'The first video frame is missing or not RGB.');
    end

    [previousCentroid, previousArea, previousAngle, found] = choose_tape_component( ...
        tape_mask(firstFrame, referenceRGB, options), [NaN, NaN], NaN, options.MaxJumpPixels);
    if ~found
        error('track_tape_motion:notDetected', ...
            ['Tape was not detected in the first frame. Reselect the reference ', ...
             'area or increase the RGB tolerance.']);
    end

    frameCount = max(1, floor(duration * frameRate));
    sampleFrames = (1 + options.FrameInterval):options.FrameInterval:frameCount;
    sampleTimes = (sampleFrames(:) - 1) / frameRate;
    sampleTimes = sampleTimes(sampleTimes < duration);

    n = numel(sampleTimes);
    results.time = sampleTimes;
    results.displacement = nan(n, 1);
    results.angleChange = nan(n, 1);
    results.detected = false(n, 1);

    previousDetected = true;
    for k = 1:n
        frame = [];
        try
            frame = getFrame(sampleTimes(k));
        catch frameError
            warning('track_tape_motion:frameRead', ...
                'Could not read the frame at %.3f s: %s', sampleTimes(k), frameError.message);
        end

        found = false;
        if ~isempty(frame)
            [centroid, area, angle, found] = choose_tape_component( ...
                tape_mask(frame, referenceRGB, options), ...
                previousCentroid, previousArea, options.MaxJumpPixels);
        end
        results.detected(k) = found;

        if found && previousDetected
            results.displacement(k) = hypot(centroid(1) - previousCentroid(1), ...
                                            centroid(2) - previousCentroid(2)) / options.PixelsPerUnit;
            results.angleChange(k) = abs(line_angle_difference(angle, previousAngle));
        end

        if found
            previousCentroid = centroid;
            previousArea = area;
            previousAngle = angle;
        end
        previousDetected = found;
    end
end


function options = apply_defaults(options)
    defaults = struct('FrameInterval', 30, 'RGBTolerance', 35, 'PixelsPerUnit', 1, ...
                      'MaxJumpPixels', 600, 'MinObjectArea', 20, 'ClosingRadius', 3);
    names = fieldnames(defaults);
    for k = 1:numel(names)
        if ~isfield(options, names{k}) || isempty(options.(names{k}))
            options.(names{k}) = defaults.(names{k});
        end
    end
    if options.FrameInterval < 1 || options.FrameInterval ~= round(options.FrameInterval)
        error('track_tape_motion:options', 'FrameInterval must be a positive whole number.');
    end
    if options.PixelsPerUnit <= 0
        error('track_tape_motion:options', 'PixelsPerUnit must be positive.');
    end
end


function mask = tape_mask(rgbFrame, referenceRGB, options)
%TAPE_MASK Pixels within the RGB tolerance on all three channels, cleaned up.
    frame = double(rgbFrame);
    tol = options.RGBTolerance;
    mask = abs(frame(:, :, 1) - referenceRGB(1)) <= tol & ...
           abs(frame(:, :, 2) - referenceRGB(2)) <= tol & ...
           abs(frame(:, :, 3) - referenceRGB(3)) <= tol;
    if options.ClosingRadius > 0
        if exist('OCTAVE_VERSION', 'builtin')
            disk = strel('disk', options.ClosingRadius, 0);   % Octave: exact disk only
        else
            disk = strel('disk', options.ClosingRadius);      % MATLAB default (as in v2/v3)
        end
        mask = imclose(mask, disk);
    end
    mask = bwareaopen(mask, options.MinObjectArea);
end


function [centroid, area, orientation, found] = choose_tape_component( ...
        mask, previousCentroid, previousArea, maxJumpPixels)
%CHOOSE_TAPE_COMPONENT Pick the blob that is most likely the tape.
    centroid = [NaN, NaN];
    area = NaN;
    orientation = NaN;
    found = false;

    if ~any(mask(:))
        return;
    end
    components = bwconncomp(mask);
    if components.NumObjects == 0
        return;
    end
    props = regionprops(components, 'Centroid', 'Area', 'Orientation');

    hasPrevious = all(isfinite(previousCentroid));
    scores = inf(numel(props), 1);
    for k = 1:numel(props)
        if ~hasPrevious
            scores(k) = -props(k).Area;            % first frame: largest blob
            continue;
        end
        jump = norm(props(k).Centroid - previousCentroid);
        if jump > maxJumpPixels
            continue;
        end
        penalty = 0;
        if isfinite(previousArea) && previousArea > 0
            ratio = props(k).Area / previousArea;
            if ratio < 0.2 || ratio > 5
                continue;
            end
            penalty = 25 * abs(log(ratio));
        end
        scores(k) = jump + penalty;               % nearest blob of similar size
    end

    [best, index] = min(scores);
    if ~isfinite(best)
        return;
    end
    centroid = props(index).Centroid;
    area = props(index).Area;
    orientation = props(index).Orientation;
    found = true;
end


function d = line_angle_difference(current, previous)
%LINE_ANGLE_DIFFERENCE Difference between two line angles (degrees), in [-90, 90).
%   A line at 89 deg and one at -89 deg differ by 2 deg, not 178 deg.
    d = mod(current - previous + 90, 180) - 90;
end
