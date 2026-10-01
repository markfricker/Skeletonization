function bw = adaptiveHysteresisSkeletonize(im, params)
%ADAPTIVEHYSTERESISSKELETONIZE  Hysteresis threshold with automatic level setting.
%
%   bw = adaptiveHysteresisSkeletonize(im, params)
%
% Extends hysteresisSkeletonize by deriving the high threshold automatically
% from the image instead of requiring a fixed, hand-tuned value, eliminating
% the main user-tuning burden.  The low threshold is set as a fixed fraction
% (ratio) of the high threshold.  This is analogous to how Canny's edge
% detector auto-sets its two thresholds.
%
% INPUTS
%   im     – 2-D single, normalised [0, 1].  Works with any upstream
%            enhance method (phase congruency, vesselness/Frangi, log,
%            granulometry, ...) — the threshold is derived from this
%            image's own histogram, not tied to a particular filter.
%   params – struct with optional fields:
%              .threshMethod – 'otsu' (default) | 'triangle' | 'triangleOtsu'.
%                           'otsu' assumes two comparable-variance classes
%                           and tends to bias high / drop dim structure on
%                           the sparse, long-tailed histograms typical of
%                           ridge-filtered images (small bright foreground
%                           fraction on a large dark background) — exactly
%                           the case 'triangle' (Zack et al., 1977) is
%                           designed for. 'triangleOtsu' runs both triangle
%                           and Otsu in log10 space on the nonzero pixels and
%                           takes the minimum -- the same Frangi-threshold
%                           recipe used by Nellie (Lefebvre et al., Nat
%                           Methods 2025); see globalThresholdFast.m.
%              .ratio     – low / high threshold ratio in (0, 1].
%                           lower ratio = more hysteresis connectivity.
%                           Default 0.4.
%              .kSigma    – noise sensitivity: high threshold is raised by
%                           kSigma * estimated noise sigma above the auto
%                           level.  Set to 0 to use the raw auto level.
%                           Default 0.
%
% OUTPUT
%   bw – logical binary mask, same size as im.
%
% DEPENDENCIES
%   globalThresholdFast (Segmentation_sandbox/src/) — only when
%   params.threshMethod = 'triangle' or 'triangleOtsu'.
%   estimateNoiseLevel (DenoiseFilters_sandbox/src/) — only when
%   params.kSigma > 0.

method = lower(string(getf(params, 'threshMethod', 'otsu')));
ratio  = getf(params, 'ratio',  0.4);
kSigma = getf(params, 'kSigma', 0.0);

% A blank image (e.g. a Z-section outside the cell after a background-gated
% Hessian filter) has no foreground. Without this the auto threshold comes
% out as 0, every pixel passes im >= 0, and bwskel turns the whole frame
% into a large X-shaped medial axis.
if ~any(im(:) > 0)
    bw = false(size(im));
    return
end

% --- Auto high threshold ---------------------------------------------------
switch method
    case "triangle"
        [~, threshHigh] = globalThresholdFast(im, 'method', 'triangle');
    case "triangleotsu"
        [~, threshHigh] = globalThresholdFast(im, 'method', 'triangleOtsu');
    case "otsu"
        threshHigh = graythresh(im);   % Otsu on [0,1] image
    otherwise
        error('adaptiveHysteresisSkeletonize:unknownMethod', ...
              'Unknown method "%s". Use ''otsu'', ''triangle'' or ''triangleOtsu''.', method);
end

% Optional noise-aware upward shift
if kSigma > 0
    % Immerkær noise estimate on the ridge image
    noiseLevel = estimateNoiseLevel(im);   % expects [0,1] image
    threshHigh = min(1, threshHigh + kSigma * noiseLevel);
end

% Same failure as the blank image: a zero threshold keeps every pixel.
if threshHigh <= 0
    bw = false(size(im));
    return
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
