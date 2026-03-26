function bw = nmsSkeletonize(im, orf, params)
%NMSSKELETONIZE  Non-maximum suppression → inherently thin binary ridge map.
%
%   bw = nmsSkeletonize(im, orf, params)
%   bw = nmsSkeletonize(im, [],  params)   % compute orientation from image
%
% Suppresses ridge-image pixels that are not local maxima along the
% direction perpendicular to the local ridge orientation.  The result is
% inherently single-pixel-wide (no bwskel needed, but skeletonPostProcess
% still runs for spur pruning and area filtering).
%
% The method is the same as the suppression step in Canny's edge detector,
% applied here to ridge images rather than gradient-magnitude images.
%
% INPUTS
%   im     – 2-D single, normalised [0, 1].  Ridge / enhanced image.
%   orf    – 2-D single orientation field in radians, same size as im.
%            Convention: orf encodes the LOCAL RIDGE DIRECTION (along the
%            tubule axis).  The perpendicular (normal) direction used for
%            NMS is therefore orf + pi/2.
%            Pass [] or zeros(size(im)) to compute orientation internally
%            from the Hessian of the ridge image.
%   params – struct with optional fields:
%              .sigma      – Gaussian smoothing applied to im before
%                            computing internal orientation (only used when
%                            orf is empty).  Default 1.5.
%              .threshold  – final intensity threshold applied after NMS.
%                            0 = Otsu automatic.  Default 0.
%
% OUTPUT
%   bw – logical binary mask, same size as im.
%
% DEPENDENCIES
%   Uses interp2 for sub-pixel sampling (MATLAB built-in).

sigma     = getf(params, 'sigma',     1.5);
threshold = getf(params, 'threshold', 0);

[nY, nX] = size(im);
[xi, yi] = meshgrid(1:nX, 1:nY);

% --- Orientation: use supplied orf or derive from Hessian ---------------
if isempty(orf) || ~any(orf(:))
    % Compute ridge normal direction from Hessian eigenvectors
    orf = hessianRidgeNormal(im, sigma);
else
    % orf encodes ridge (along-tubule) direction; normal is perpendicular
    orf = orf + pi/2;
end

% --- NMS: sample ±1 pixel along the ridge normal ----------------------
dx = cos(orf);
dy = sin(orf);

ip = interp2(im, xi + dx, yi + dy, 'linear', 0);
im_ = interp2(im, xi - dx, yi - dy, 'linear', 0);

bw = (im >= ip) & (im >= im_) & (im > 0);

% --- Threshold on ridge strength ----------------------------------------
if threshold <= 0
    threshold = graythresh(im);
end
bw = bw & (im >= threshold);
end

% ---- Hessian-based ridge normal direction --------------------------------
function normalAngle = hessianRidgeNormal(im, sigma)
% Returns per-pixel angle of the Hessian eigenvector corresponding to the
% most negative eigenvalue (the ridge normal direction for bright ridges).
sigma = double(sigma);
G   = imgaussfilt(double(im), sigma);
[Gx,  ~]  = gradient(G);
[Gxx, Gxy] = gradient(Gx);
[~,  Gyy] = gradient(gradient(G));

% Symmetrise cross-derivative
[~, Gyx] = gradient(imgaussfilt(double(im), sigma));
[~, Gyy2] = gradient(Gyx);
Gyy = (Gyy + Gyy2) / 2;

% Angle of the minor eigenvector (ridge normal) for 2×2 Hessian:
%   θ = 0.5 * atan2(2*Gxy, Gxx - Gyy)
normalAngle = single(0.5 * atan2(2*Gxy, Gxx - Gyy));
end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
