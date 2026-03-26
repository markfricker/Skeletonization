function bw = adaptiveHysteresisSkeletonize(im, params)
%ADAPTIVEHYSTERESISSKELETONIZE  Hysteresis threshold with automatic level setting.
%
%   bw = adaptiveHysteresisSkeletonize(im, params)
%
% Extends hysteresisSkeletonize by deriving the high threshold automatically
% from the image using Otsu's method, eliminating the main user-tuning
% burden.  The low threshold is set as a fixed fraction (ratio) of the
% high threshold.  This is analogous to how Canny's edge detector auto-sets
% its two thresholds.
%
% INPUTS
%   im     – 2-D single, normalised [0, 1].
%   params – struct with optional fields:
%              .ratio     – low / high threshold ratio in (0, 1].
%                           lower ratio = more hysteresis connectivity.
%                           Default 0.4.
%              .kSigma    – noise sensitivity: high threshold is raised by
%                           kSigma * estimated noise sigma above Otsu level.
%                           Set to 0 to use raw Otsu.  Default 0.
%
% OUTPUT
%   bw – logical binary mask, same size as im.

ratio  = getf(params, 'ratio',  0.4);
kSigma = getf(params, 'kSigma', 0.0);

% --- Auto high threshold via Otsu ----------------------------------------
threshHigh = graythresh(im);   % Otsu on [0,1] image

% Optional noise-aware upward shift
if kSigma > 0
    % Immerkær noise estimate on the ridge image
    noiseLevel = estimateNoiseLevel(im);   % expects [0,1] image
    threshHigh = min(1, threshHigh + kSigma * noiseLevel);
end

% --- Low threshold -------------------------------------------------------
threshLow = ratio * threshHigh;

% --- Hysteresis binarization ---------------------------------------------
bwHigh = im >= threshHigh;
bwLow  = im >= threshLow;
bw     = imreconstruct(bwHigh, bwLow);
end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
