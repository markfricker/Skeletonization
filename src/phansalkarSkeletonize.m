function bw = phansalkarSkeletonize(im, params)
%PHANSALKARSKELETONIZE  Phansalkar local threshold → binary foreground mask.
%
%   bw = phansalkarSkeletonize(im, params)
%
% Applies the Phansalkar (2011) adaptive threshold — an extension of
% Sauvola optimised for low-contrast fluorescence microscopy images.  The
% key addition is an exponential term that boosts the threshold in dim
% regions, preventing the collapse of Sauvola's formula when local mean
% intensity is near zero (as is common in sparse ER network images).
%
% THRESHOLD FORMULA
%   T(x,y) = μ(x,y) · [ 1 + p · exp(−q · μ(x,y)) + k · (σ(x,y)/R − 1) ]
%
%   where p and q govern how strongly the threshold is boosted in low-mean
%   regions.  Typical defaults: p=2, q=10.
%
% INPUTS
%   im     – 2-D single, normalised [0, 1].
%   params – struct with optional fields:
%              .windowSize – neighbourhood size in pixels (odd integer).
%                            Default 31.
%              .k          – Sauvola sensitivity weight.  Default 0.25.
%              .r          – dynamic range of σ for [0,1] images.  Default 0.5.
%              .p          – exponential boost amplitude.  Default 2.0.
%              .q          – exponential decay rate.  Default 10.0.
%
% OUTPUT
%   bw – logical binary mask, same size as im.
%
% REFERENCE
%   Phansalkar, N., More, S., Sabale, A., & Joshi, M. (2011). Adaptive
%   local thresholding for detection of nuclei in diversity stained
%   cytology images. In Proc. ICCSP.

windowSize = getf(params, 'windowSize', 31);
k          = getf(params, 'k',          0.25);
r          = getf(params, 'r',          0.5);
p          = getf(params, 'p',          2.0);
q          = getf(params, 'q',          10.0);

% Ensure odd window size
windowSize = 2*floor(windowSize/2) + 1;

% --- Local mean and variance via integral images -------------------------
mu_d  = imboxfilt(double(im), windowSize);
mu2_d = imboxfilt(double(im).^2, windowSize);
sigma = sqrt(max(mu2_d - mu_d.^2, 0));

% --- Phansalkar threshold ------------------------------------------------
%   Extra exponential term raises T when mu is near 0 (dim regions)
T  = mu_d .* (1 + p .* exp(-q .* mu_d) + k .* (sigma ./ r - 1));
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
