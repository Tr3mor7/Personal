function test_track_tape_motion()
%TEST_TRACK_TAPE_MOTION Check the tracker on synthetic frames with known motion.
%
%   Run from the tape-motion-tracking folder:
%       addpath('tests'); test_track_tape_motion
%   Works in MATLAB (Image Processing Toolbox) and GNU Octave (pkg load image).
%
%   A 60 x 14 px orange strip moves 40 px right and 10 px down and rotates
%   12 degrees every second. A small orange square elsewhere must be ignored,
%   and an unreadable frame at t = 3 s must give NaN without stopping the run.

    if exist('OCTAVE_VERSION', 'builtin')
        pkg('load', 'image');
    end
    here = fileparts(mfilename('fullpath'));
    addpath(fileparts(here));

    fps = 30;
    duration = 7;
    tapeRGB = [230, 120, 40];
    getFrame = @(t) synthetic_frame(t, tapeRGB);

    options = struct('FrameInterval', 30, 'RGBTolerance', 35, 'PixelsPerUnit', 2);
    r = track_tape_motion(getFrame, duration, fps, tapeRGB, options);

    expectedDisplacement = hypot(40, 10) / options.PixelsPerUnit;   % calibrated units
    assert_close(r.time, (1:6)', 1e-9, 'sample times');
    assert_close(r.displacement([1 2 5 6]), expectedDisplacement * ones(4, 1), 0.75, 'displacement');
    assert_close(r.angleChange([1 2 5 6]), 12 * ones(4, 1), 1.5, 'angle change');

    % one unreadable frame: that sample and the next one have no value
    getFrameWithGap = @(t) frame_or_error(t, tapeRGB, 3);
    warningState = warning('off', 'track_tape_motion:frameRead');
    r2 = track_tape_motion(getFrameWithGap, duration, fps, tapeRGB, options);
    warning(warningState);
    assert(isequal(r2.detected', [true true false true true true]), 'detected flags with a gap');
    assert(all(isnan(r2.displacement([3 4]))) && all(isnan(r2.angleChange([3 4]))), 'NaN after a gap');
    assert_close(r2.displacement(5), expectedDisplacement, 0.75, 'recovery after a gap');

    % line angles wrap: a strip turning from +85 to -85 deg has turned 10 deg, not 170
    r3 = track_tape_motion(@(t) strip_frame([100, 100], 85 - 170 * (t > 0.5), tapeRGB), ...
                           2, fps, tapeRGB, options);
    % (a thin strip near vertical is pixelated: 85 deg measures ~83.9 deg, so allow 3 deg)
    assert_close(r3.angleChange(1), 10, 3, 'angle wrap-around');

    % helpful error when the tape colour is not in the first frame
    try
        track_tape_motion(getFrame, duration, fps, [0, 255, 0], options);
        error('test:noError', 'expected an error for a colour that is not present');
    catch err
        assert(~isempty(strfind(err.message, 'not detected')), 'error message for missing tape');
    end

    fprintf('test_track_tape_motion: all checks passed\n');
end


function frame = synthetic_frame(t, tapeRGB)
    frame = strip_frame([60 + 40 * t, 100 + 10 * t], 5 + 12 * t, tapeRGB);
    frame(200:207, 30:37, 1) = tapeRGB(1);   % small orange distractor (64 px)
    frame(200:207, 30:37, 2) = tapeRGB(2);
    frame(200:207, 30:37, 3) = tapeRGB(3);
end


function frame = frame_or_error(t, tapeRGB, badTime)
    if abs(t - badTime) < 1e-9
        error('test:unreadable', 'simulated unreadable frame');
    end
    frame = synthetic_frame(t, tapeRGB);
end


function frame = strip_frame(centre, angleDeg, tapeRGB)
%STRIP_FRAME 240 x 400 textured background with a rotated 60 x 14 strip.
    [x, y] = meshgrid(1:400, 1:240);
    texture = 8 * sin(x / 7) .* cos(y / 5);
    background = cat(3, 60 + texture, 70 + texture, 90 + texture);
    a = angleDeg * pi / 180;
    u = [cos(a), -sin(a)];            % long axis (image y points down)
    v = [sin(a), cos(a)];
    corners = [centre + 30 * u + 7 * v; centre - 30 * u + 7 * v; ...
               centre - 30 * u - 7 * v; centre + 30 * u - 7 * v];
    inside = poly2mask(corners(:, 1), corners(:, 2), 240, 400);
    frame = background;
    for c = 1:3
        channel = frame(:, :, c);
        channel(inside) = tapeRGB(c) + 0.6 * texture(inside);
        frame(:, :, c) = channel;
    end
    frame = uint8(frame);
end


function assert_close(actual, expected, tolerance, label)
    if any(abs(actual(:) - expected(:)) > tolerance) || numel(actual) ~= numel(expected)
        error('test:mismatch', '%s: got [%s], expected [%s] (tolerance %g)', label, ...
              num2str(actual(:)', '%.3f '), num2str(expected(:)', '%.3f '), tolerance);
    end
end
