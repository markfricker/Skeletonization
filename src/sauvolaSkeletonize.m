function bw = sauvolaSkeletonize(im, params)
%SAUVOLASKELETONIZE  Sauvola local threshold → binary foreground mask.
%
%   bw = sauvolaSkeletonize(im, params)
%
% Applies the Sauvola (1999) adaptive threshold to the enhanced image.
% Sauvola is particularly effective for images with spatially varying
% background because the threshold adapts to the local mean and standard
% deviation.
%
% THRESHOLD FORMULA
%   T(x,y) = μ(x,y) · [ 1 + k · ( σ(x,y)/R − 1 ) ]
%
%   where μ and σ are the local mean and standard deviation computed over
%   a square neighbourhood of side windowSize, R is the dynamic range of σ
%   (nominally 0.5 for images in [0,1]), and k controls sensitivity.
%
%   Behaviour:
%   - High local σ (edge/ridge regions): T < μ  →  easier to exceed threshold
%   - Low local σ (flat background):     T ≈ μ(1-k) →  harder to exceed
%
% INPUTS
%   im     – 2-D single, normalised [0, 1].
%   params – struct with optional fields:
%              .windowSize – neighbourhood size in pixels (odd integer).
%                            Default 31.
%              .k          – sensitivity parameter.  Larger k raises the
%                            threshold in textured regions.  Default 0.2.
%              .r          – dynamic range of σ.  Use 0.5 for [0,1] images.
%                            Default 0.5.
%
% OUTPUT
%   bw – logical binary mask, same size as im.
%
% REFERENCE
%   Sauvola, J. & Pietikainen, M. (2000). Adaptive document image
%   binarization. Pattern Recognition, 33(2), 225-236.

windowSize = getf(params, 'windowSize', 31);
k          = getf(params, 'k',          0.2);
r          = getf(params, 'r',          0.5);

% Ensure odd window size
windowSize = 2*floor(windowSize/2) + 1;

% --- Local mean and variance via integral images (imboxfilt) ------------
mu    = imboxfilt(double(im), windowSize);
mu2   = imboxfilt(double(im).^2, windowSize);
sigma = sqrt(max(mu2 - mu.^2, 0));

% --- Sauvola threshold ---------------------------------------------------
T  = mu .* (1 + k .* (sigma ./ r - 1));
bw = im > single(T);
end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
