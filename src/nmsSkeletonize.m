function bw = nmsSkeletonize(im, orf, params)
%NMSSKELETONIZE  Non-maximum suppression → inherently thin binary ridge map.
%
%   bw = nmsSkeletonize(im, orf, params)
%   bw = nmsSkeletonize(im, [],  params)   % compute orientation from image
%
% Uses Kovesi's featureorient + smoothorient + nonmaxsup pipeline to
% suppress pixels that are not local maxima along the direction
% perpendicular to the local ridge orientation.
%
% INPUTS
%   im     – 2-D single, normalised [0, 1].  Ridge / enhanced image.
%   orf    – 2-D single orientation field in radians, same size as im.
%            Convention: encodes the LOCAL RIDGE DIRECTION (along the
%            tubule axis).  Converted internally to normal direction in
%            degrees for nonmaxsup.
%            Pass [] or zeros(size(im)) to compute orientation via
%            featureorient (same parameters as WS+NMS).
%   params – struct with optional fields:
%              .radius     – sampling radius for nonmaxsup (pixels).
%                            Values 1.2–1.5 avoid broad-peak misses;
%                            larger values (e.g. 3) match WS+NMS behaviour.
%                            Default: 1.5.
%              .threshold  – final ridge-strength threshold after NMS.
%                            0 = Otsu automatic.  Default: 0.
%
% OUTPUT
%   bw – logical binary mask, same size as im.
%
% DEPENDENCIES (must be on MATLAB path)
%   Common_sandbox/Kovesi phase congruency/  – featureorient, smoothorient,
%                                              nonmaxsup

radius    = getf(params, 'radius',    1.5);
threshold = getf(params, 'threshold', 0);

% --- Orientation in degrees [0,180] for nonmaxsup (feature normal) ------
if isempty(orf) || ~any(orf(:))
    % Compute normal-to-ridge orientation directly from image intensity.
    % featureorient returns degrees [0,180] across-ridge (normal direction).
    or = featureorient(double(im), 0, 1, 3, 0);
else
    % orf is ridge direction (along-tubule) in radians.
    % Normal direction = orf + pi/2 → convert to degrees → wrap to [0,180].
    or = mod((orf + pi/2) * (180/pi), 180);
end

% Smooth the orientation field to reduce noise-driven direction flips
or = smoothorient(or, 1.5);

% --- Kovesi NMS ---------------------------------------------------------
nmsOut = nonmaxsup(double(im), or, radius);

% --- Threshold on ridge strength ----------------------------------------
if threshold <= 0
    threshold = graythresh(im);
end
bw = logical(nmsOut) & (im >= single(threshold));
end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
