# tape-motion-tracking: batch video tracking of a coloured tape marker (MATLAB)

Put a strip of brightly coloured tape (e.g. orange) on the object you're
filming. For each video, this tool measures how far the tape moves and how
much it rotates every *N* frames, and writes the results to Excel. It handles
many videos in one batch: you set everything up first, then it runs
unattended.

## How it works

1. **Colour reference.** On a screenshot or photo from each video, you draw a
   rectangle inside the tape. Its mean RGB colour is what gets tracked.
2. **Segmentation.** In every sampled frame, pixels within the RGB tolerance
   of that colour on all three channels are kept. Small gaps are closed, and
   specks smaller than 20 px are removed.
3. **Picking the tape.** In the first frame, the largest matching blob is the
   tape. After that, it's the blob closest to the previous position, with an
   area within 5× of the previous area. This ignores other orange things in
   the scene.
4. **Measurements, every N frames:**
   * **displacement:** distance the blob's centroid moved since the previous
     sample, in pixels or calibrated units;
   * **angle change:** absolute change of the blob's long-axis angle, in
     degrees. A line has no direction, so +89° → −89° counts as 2°, not 178°.
5. **Missing frames.** If a frame can't be read, or the tape isn't found,
   that sample is recorded as empty (`NaN`) and the run carries on.

## Requirements

* **MATLAB R2020a or newer** with the **Image Processing Toolbox**.
* The core tracker (`track_tape_motion.m`) and its test also run in
  **GNU Octave** with the `image` package. The batch front end needs MATLAB,
  for its dialogs, `VideoReader` and Excel output.

## How to run it

```matlab
cd tape-motion-tracking
tape_motion_batch
```

**Global setup**, asked once and used for every run:

| Setting | Notes |
|---|---|
| number of runs | one run = one video |
| RGB tolerance (0–255) | default 35; raise it if the tape isn't detected, lower it if the background is picked up |
| units | pixels, or a calibration such as `10` pixels per `mm` |
| frame interval *N* | e.g. 30 = once per second for a 30 fps video; the column names include *N* |
| Excel output | **combined**: one workbook with a sheet per run; or **per run**: one workbook in each run's output folder |

**Per-run setup:** the video, a reference screenshot/photo, an output
folder, and the rectangle drawn inside the tape.

Analysis then runs without prompts. If something goes wrong (an unreadable
video, an open workbook, tape not found), you're asked whether to retry,
reconfigure or skip that run, instead of losing the whole batch.

## Inputs and outputs

| | What | Notes |
|---|---|---|
| **In** | video per run | `.mp4`, `.avi`, `.mov`, `.m4v`, `.mpg` (anything `VideoReader` reads); colour |
| **In** | reference image per run | a frame or photo in which the tape is visible |
| **Out** | Excel results | columns `Time_seconds`, `DisplacementMagnitude_EveryNFrames_<unit>`, `AbsoluteAngleChange_EveryNFrames_degrees` |
| **Out** | `<photo>_tape.png` | the reference image with the selected rectangle marked, as a record of what was tracked |

**Workbook names:** a workbook is named after the parent of the folder it's
written to (`…/Trial_03/Output/` → `Trial_03.xlsx`). In combined mode, the
sheets are `Run01_<video name>`, `Run02_…`. A run that fails or is skipped
still gets an empty sheet, so there's always one sheet per run.

## Scripting without dialogs

`track_tape_motion` is the analysis on its own. For example, in a script:

```matlab
v = VideoReader('my_clip.mp4');
getFrame = @(t) read_at(v, t);          % function that returns the frame at time t
r = track_tape_motion(getFrame, v.Duration, v.FrameRate, [230 120 40], ...
                      struct('FrameInterval', 30, 'RGBTolerance', 35));
plot(r.time, r.displacement)

function f = read_at(v, t)
    v.CurrentTime = t;
    f = readFrame(v);
end
```

## Tests

```matlab
cd tape-motion-tracking
addpath('tests'); test_track_tape_motion
```

The test builds synthetic frames in memory, so no video files are needed: an
orange strip moves and rotates by a known amount every second, next to a
small orange distractor. It checks:

* the displacement and angle change, including the calibration;
* that the distractor is ignored;
* that an unreadable frame gives `NaN` and the run recovers;
* angle wrap-around;
* the error message when the colour isn't in the first frame.

It passes in GNU Octave 8.4.

## History

This version merges three earlier scripts (`batch_orange_tape_motion_analysis_v2.m`,
`…_v2_90frames.m` and `…_v3.m`):

* from **v3**: input validation, and the retry/reconfigure/skip handling;
* from **v2_90frames**: the combined workbook with one sheet per run;
* **new**:
  * the frame interval is a setting, instead of 30 or 90 fixed in the code;
  * the analysis was moved into `track_tape_motion.m`, so it can be tested;
  * it's a function, so it no longer clears your workspace.

It also fixes a labelling bug in v2, which sampled every **90** frames but
named its columns `…Every30Frames…`.

## AI disclosure

This code was written with the support of Claude Code (Anthropic's AI coding
assistant), starting from scripts by Gurnoor Hayer. The project and its design
were directed and reviewed by Gurnoor Hayer.
