function bw = phaseCongSkeletonize(im, params)
%PHASECONGSKELETONIZE  Phase congruency ridge detection → binary mask.
%
%   bw = phaseCongSkeletonize(im, params)
%
% Applies Kovesi's phase congruency (phasecong3) to the input image and
% thresholds the maximum-moment response to produce a binary ridge mask.
% Phase congruency detects features based on the consistency of Fourier
% phase across scales, making it completely invariant to absolute intensity
% and robust to intensity gradients across the image.
%
% This is particularly effective for ER network images where fluorophore
% density varies significantly across the field.
%
% DEPENDENCY
%   Requires phasecong3.m (Kovesi, 1996-2012) on the MATLAB path.
%   Typically located in Common_sandbox/Kovesi phase congruency/.
%
% INPUTS
%   im     – 2-D single, normalised [0, 1].
%   params – struct with optional fields:
%              .nscale    – number of log-Gabor wavelet scales.  Default 4.
%              .norient   – number of filter orientations.  Default 6.
%              .minWL     – minimum filter wavelength (pixels).
%                           Controls the finest scale.  Default 3.
%              .mult      – scale factor between successive filter bands.
%                           Default 2.1.
%              .sigmaOnf  – log-Gabor filter bandwidth (ratio of σ to
%                           centre frequency).  Default 0.55.
%              .k         – noise sensitivity: noise threshold is set at
%                           k standard deviations above the mean noise
%                           energy.  Higher k = less sensitive.  Default 2.
%              .threshold – final binarisation level applied to the phase
%                           congruency maximum moment M.  0 = use the
%                           data-driven noise threshold T returned by
%                           phasecong3 (recommended).  Default 0.
%
% OUTPUT
%   bw – logical binary mask, same size as im.
%
% REFERENCE
%   Kovesi, P. (1999). Image features from phase congruency.
%   Videre: Journal of Computer Vision Research, 1(3).

nscale    = getf(params, 'nscale',    4);
norient   = getf(params, 'norient',   6);
minWL     = getf(params, 'minWL',     3);
mult      = getf(params, 'mult',      2.1);
sigmaOnf  = getf(params, 'sigmaOnf',  0.55);
k         = getf(params, 'k',         2.0);
threshold = getf(params, 'threshold', 0);

% phasecong3 requires double input
imD = double(im);

% Call phasecong3: [M m or ft pc EO T]
% M = max moment (edge/ridge strength); T = data-driven noise threshold
[M, ~, ~, ~, ~, ~, T] = phasecong3(imD, nscale, norient, minWL, mult, sigmaOnf, k);

% Select threshold
if threshold <= 0
    applyThresh = T;
else
    applyThresh = threshold * max(M(:) + eps);
end

bw = M > applyThresh;
end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
