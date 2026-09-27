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
%   2. Local maxima of the response (lmax) -- Sten et al. use NMS in the
%      max-response orientation; here imregionalmax by default (see below).
%   3. Otsu threshold on the distribution of locally-maximal response values
%      only → seed map (BW2).  This is more robust than a fixed high threshold
%      because it adapts to the image-wide response distribution.
%   4. 8-connected hysteresis: grow BW2 seeds into BW1 (imfill with seeds).
%
% Step 2 (seed candidates) defaults to imregionalmax, which is what
% AnalyzER has always run: the MycNetAnalysis nonMaxSuppression this
% originally called was never reached, because nothing supplied
% params.dirMap. Orientation-guided NMS (Kovesi nonmaxsup, sub-pixel
% interpolated) is available opt-in. On a synthetic ground-truth test
% (2026-09-27) NMS seeds did NOT improve the final mask (F1 0.41 vs 0.45 at
% low noise, ~equal at high noise): NMS keeps whole noise ridges as seed
% candidates, so more noise components survive the hysteresis step.
%
% INPUTS
%   im     - 2-D single, normalised [0,1].  Typically ridge filter output.
%   params - struct with optional fields:
%              .sensitivity - Bradley adaptive threshold sensitivity [0,1].
%                             Default 0.5.
%              .nmsUse      - true: seed candidates from orientation-guided
%                             NMS instead of imregionalmax. Default false.
%              .dirMap      - feature-normal orientation for NMS (degrees,
%                             +ve anticlockwise, taken mod 180), e.g.
%                             steerGaussEnhance's dirMap. Supplying it
%                             implies nmsUse. If NMS is on and dirMap is
%                             absent, orientation is estimated from im with
%                             featureorient (as wsNmsSkeletonize does).
%
% OUTPUT
%   bw - logical binary mask, same size as im.
%
% Dependencies
%   Kovesi phase congruency toolkit: featureorient, nonmaxsup
%   (Common_sandbox/Kovesi phase congruency) -- NMS path only.
%
% See also: hysteresisSkeletonize, wsNmsSkeletonize, steerGaussEnhance,
%           featureorient, nonmaxsup

sensitivity = getf(params, 'sensitivity', 0.5);
dirMap      = getf(params, 'dirMap',      []);
nmsUse      = getf(params, 'nmsUse',      false) || ~isempty(dirMap);

nmsRadius = 1.5;   % Kovesi's suggested 1.2-1.5 (avoids missed maxima at 1)

% Step 1: normalise and adaptive threshold
responseMap = mat2gray(double(im));
BW1 = imbinarize(responseMap, 'adaptive', 'Sensitivity', sensitivity);

% Step 2: seed candidates -- local maxima of the response
if nmsUse
    if isempty(dirMap)
        orient = featureorient(responseMap, 0, 1, 0, 0);   % degrees, 0-180
    else
        orient = mod(double(dirMap), 180);
    end
    lmax = nonmaxsup(responseMap, orient, nmsRadius) > 0;
else
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
