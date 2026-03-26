function bw = hessianCenterlineSkeletonize(im, params)
%HESSIANCENTERLINESKELETONIZE  Sub-pixel Hessian centerline localization.
%
%   bw = hessianCenterlineSkeletonize(im, params)
%
% Computes the 2-D Hessian matrix at each pixel, identifies ridges via the
% sign of the smallest eigenvalue, then applies non-maximum suppression
% along the eigenvector direction (the ridge normal).  This gives the most
% spatially accurate centerline localization of all methods because the
% ridge normal is estimated directly from the image curvature — no
% separate orientation image is required.
%
% For bright tubules on a dark background, the most negative eigenvalue
% λ₂ indicates ridge strength; its eigenvector points across the tubule
% (the normal direction for NMS).
%
% PIPELINE
%   1. Scale-normalised Gaussian second-derivative (Hessian) at scale σ.
%   2. Eigenvalue decomposition: ridge strength = max(−λ₂, 0).
%   3. NMS: retain pixels where strength ≥ interpolated values at ±1 px
%      along the eigenvector of λ₂.
%   4. Threshold on ridge strength.
%
% INPUTS
%   im     – 2-D single, normalised [0, 1].
%   params – struct with optional fields:
%              .sigma      – Gaussian derivative scale (pixels).  Should
%                            match the expected tubule radius.  Default 1.5.
%              .threshold  – ridge strength threshold [0, 1].
%                            0 = Otsu on the strength map.  Default 0.
%
% OUTPUT
%   bw – logical binary mask, same size as im.

sigma     = double(getf(params, 'sigma',     1.5));
threshold =        getf(params, 'threshold', 0);

imD = double(im);
s2  = sigma^2;   % scale normalisation factor

% --- Scale-normalised Hessian via Gaussian smoothing + finite differences
G           = imgaussfilt(imD, sigma);
[Gx,  Gy]  = gradient(G);
[Gxx, Gxy] = gradient(Gx);
[Gyx, Gyy] = gradient(Gy);

% Scale-normalise and symmetrise
Gxx = s2 * Gxx;
Gyy = s2 * Gyy;
Gxy = s2 * (Gxy + Gyx) / 2;

% --- Eigenvalues of 2×2 Hessian (closed form) ---------------------------
tr   = Gxx + Gyy;
disc = sqrt(max((Gxx - Gyy).^2 + 4*Gxy.^2, 0));

lam1 = (tr + disc) / 2;   % larger eigenvalue
lam2 = (tr - disc) / 2;   % smaller eigenvalue  (negative = bright ridge)

% Ridge strength: most-negative eigenvalue, clamped to ≥ 0
strength = single(max(-lam2, 0));

% Normalise strength to [0,1] for threshold consistency
sMax = max(strength(:));
if sMax > 0
    strength = strength / sMax;
end

% --- Eigenvector of λ₂: the ridge normal direction ----------------------
% For the 2×2 symmetric Hessian, the angle of the eigenvector of the
% smaller eigenvalue is:  θ = atan2(Gxy, lam2 - Gyy)
theta = atan2(Gxy, lam2 - Gyy);

% --- Non-maximum suppression along the ridge normal ---------------------
[nY, nX] = size(im);
[xi, yi] = meshgrid(1:nX, 1:nY);

dx = cos(theta);
dy = sin(theta);

sp = interp2(strength, xi + dx, yi + dy, 'linear', 0);
sm = interp2(strength, xi - dx, yi - dy, 'linear', 0);

bw = (strength >= sp) & (strength >= sm) & (strength > 0);

% --- Threshold on ridge strength ----------------------------------------
if threshold <= 0
    threshold = graythresh(strength);
end
bw = bw & (strength >= threshold);
end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
