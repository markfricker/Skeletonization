function bw = watershedSkeletonize(im, params)
%WATERSHEDSKELETONIZE  Marker-controlled watershed → binary ridge mask.
%
%   bw = watershedSkeletonize(im, params)
%
% Uses a marker-controlled watershed on the distance transform of an
% initial threshold mask to separate touching ridges.  Unlike the full
% watershedSegment (which returns labelled objects), this function returns
% a binary mask — the ridge foreground with watershed-derived separation
% lines removed — suitable for subsequent skeletonisation.
%
% PIPELINE
%   1. Optional Gaussian smoothing of im.
%   2. Hysteresis threshold → binary foreground mask (bwFg).
%   3. Distance transform of bwFg → local maxima are ridge centres.
%   4. h-minima suppressed watershed on −distance → separation lines.
%   5. Return bwFg with separation lines removed (watershed(L)==0 pixels
%      that lie inside bwFg are suppressed).
%
% INPUTS
%   im     – 2-D single, normalised [0, 1].
%   params – struct with optional fields:
%              .threshold   – foreground binarisation level [0,1].
%                             0 = Otsu automatic.  Default 0.3.
%              .threshLow   – lower hysteresis threshold.  0 = disabled.
%                             Default 0.
%              .smoothSigma – Gaussian pre-smoothing sigma (pixels).
%                             0 = no smoothing.  Default 1.0.
%              .hMinima     – h-minima suppression depth for the distance
%                             transform.  Merges basins shallower than
%                             this value.  0 = no suppression.  Default 2.
%
% OUTPUT
%   bw – logical binary mask, same size as im.

threshold   = getf(params, 'threshold',   0.3);
threshLow   = getf(params, 'threshLow',   0);
smoothSigma = getf(params, 'smoothSigma', 1.0);
hMinima     = getf(params, 'hMinima',     2);

% --- Optional smoothing --------------------------------------------------
if smoothSigma > 0
    im = imgaussfilt(im, double(smoothSigma));
end

% --- Foreground mask (with optional hysteresis) -------------------------
if threshold <= 0
    threshold = graythresh(im);
end
bwHigh = im >= threshold;
if threshLow > 0 && threshLow < threshold
    bwLow = im >= threshLow;
    bwFg  = imreconstruct(bwHigh, bwLow);
else
    bwFg = bwHigh;
end

if ~any(bwFg(:))
    bw = false(size(im));
    return
end

% --- Distance transform and watershed ------------------------------------
D = bwdist(~bwFg);     % distance to nearest background pixel

% Negate so ridges = minima (watershed finds basins = minima)
Dneg = -D;

% h-minima suppression: merges basins shallower than hMinima
if hMinima > 0
    Dneg = imhmin(Dneg, hMinima);
end

% Force background to Inf so watershed stays within bwFg
Dneg(~bwFg) = Inf;

L  = watershed(Dneg);

% Remove watershed ridge lines (L==0) from the foreground mask
bw = bwFg & (L > 0);
end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
