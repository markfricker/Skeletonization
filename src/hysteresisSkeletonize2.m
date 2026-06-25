function bw = hysteresisSkeletonize2(im, params)
%HYSTERESISSKELETONIZE2  NMS + Otsu-on-local-maxima + hysteresis binarisation.
%
%   bw = hysteresisSkeletonize2(im, params)
%
% Implements the binarisation pipeline from:
%   Sten et al. (2024) "A Ridge-Based Detection Algorithm with Filament
%   Overlap Identification for 2D Mycelium Network Analysis."
%   DOI: 10.1016/j.ecoinf.2024.102670
%
% This is a more principled alternative to plain hysteresisSkeletonize:
%   1. Bradley adaptive threshold on the normalised response map (BW1).
%   2. Non-maximum suppression (NMS) in the max-response orientation to find
%      ridge centreline pixels (lmax).
%   3. Otsu threshold on the distribution of locally-maximal response values
%      only → seed map (BW2).  This is more robust than a fixed high threshold
%      because it adapts to the image-wide response distribution.
%   4. 8-connected hysteresis: grow BW2 seeds into BW1 (imfill with seeds).
%
% Works best when im is the output of steerGaussEnhance (or any ridge filter
% that also returns a per-pixel orientation map in params.dirMap).  If
% params.dirMap is absent, NMS is skipped and step 3 uses a global Otsu seed.
%
% INPUTS
%   im     - 2-D single, normalised [0,1].  Typically ridge filter output.
%   params - struct with optional fields:
%              .dirMap      - orientation map (degrees) from the ridge filter,
%                             same size as im.  If absent, NMS is skipped.
%              .sensitivity - Bradley adaptive threshold sensitivity [0,1].
%                             Default 0.5.
%
% OUTPUT
%   bw - logical binary mask, same size as im.
%
% Dependencies
%   nonMaxSuppression  (MycNetAnalysis/src/fcn/external/ must be on path)
%
% See also: hysteresisSkeletonize, steerGaussEnhance, nonMaxSuppression

sensitivity = getf(params, 'sensitivity', 0.5);
dirMap      = getf(params, 'dirMap',      []);

% Step 1: normalise and adaptive threshold
responseMap = mat2gray(double(im));
BW1 = imbinarize(responseMap, 'adaptive', 'Sensitivity', sensitivity);

% Step 2: NMS to find ridge centreline pixels
if ~isempty(dirMap)
    lmax = nonMaxSuppression(double(im), double(dirMap));
else
    % No orientation map — use a local max filter as fallback
    lmax = imregionalmax(im);
end

% Step 3: Otsu threshold on locally-maximal response values only
lmaxRespMap          = responseMap;
lmaxRespMap(~lmax)   = 0;
lmaxRespList         = lmaxRespMap(:);
lmaxRespList(lmaxRespList == 0) = [];

if ~isempty(lmaxRespList)
    T    = otsuthresh(imhist(lmaxRespList));
    BW2  = lmaxRespMap > T;
else
    BW2  = lmax;
end

% Step 4: 8-connected hysteresis — grow BW2 seeds into BW1
seedIdx = find(BW2);
if ~isempty(seedIdx)
    hys = imfill(~BW1, seedIdx, 8);
    bw  = hys & BW1;
else
    bw = BW1;
end

end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
