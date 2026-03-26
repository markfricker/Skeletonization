function bw = hysteresisSkeletonize(im, params)
%HYSTERESISSKELETONIZE  Two-level hysteresis threshold → binary foreground mask.
%
%   bw = hysteresisSkeletonize(im, params)
%
% Applies Canny-style hysteresis connectivity analysis to the ridge/enhanced
% image: pixels above threshHigh are confirmed foreground; pixels between
% threshLow and threshHigh are included if they are connected (8-connected)
% to a confirmed pixel.  The result is a binary mask ready for
% skeletonPostProcess.
%
% INPUTS
%   im     – 2-D single, normalised [0, 1].  Typically the output of an
%            enhancement step (ridge filter, phase congruency, etc.).
%   params – struct with fields:
%              .threshHigh  – high (confirmation) threshold in [0,1].
%                             Pixels above this are unconditionally foreground.
%                             Default 0.5.
%              .threshLow   – low (connectivity) threshold in [0,1].
%                             Must be <= threshHigh. Default 0.2.
%
% OUTPUT
%   bw – logical binary mask, same size as im.
%
% ALGORITHM
%   imreconstruct(bwHigh, bwLow) performs geodesic dilation of the high-
%   confidence seed (bwHigh) within the permissive mask (bwLow), which is
%   the exact binary equivalent of hysteresis thresholding.

threshHigh = getf(params, 'threshHigh', 0.5);
threshLow  = getf(params, 'threshLow',  0.2);

% Clamp threshLow so it never exceeds threshHigh
threshLow = min(threshLow, threshHigh);

bwHigh = im >= threshHigh;
bwLow  = im >= threshLow;

% Hysteresis: grow high-confidence seeds into low-confidence regions
bw = imreconstruct(bwHigh, bwLow);
end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
