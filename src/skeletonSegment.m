function [skel, bw] = skeletonSegment(im, params)
%SKELETONSEGMENT  Unified ER skeleton extraction: segment → skeletonise.
%
%   [skel, bw] = skeletonSegment(im, params)
%   [skel, bw] = skeletonSegment(im, params, orf)
%
% Single entry-point for the Skeletonization_sandbox.  Routes the
% enhanced image to the selected segmentation algorithm, then applies
% shared post-processing (area filter → bwskel) via skeletonPostProcess.
%
% INPUTS
%   im     – 2-D single image, normalised [0, 1].  Typically the output of
%            an enhancement step (ridge filter, phase congruency, etc.).
%   params – struct.  Required field:
%              .method  – string, one of:
%                           'none'
%                           'hysteresis'
%                           'adaptiveHysteresis'
%                           'nms'
%                           'sauvola'
%                           'phansalkar'
%                           'watershed'
%                           'phaseCong'
%                           'hessianCenterline'
%            Optional top-level post-processing fields:
%              .pruneLength  – bwskel MinBranchLength (default 10)
%              .minArea      – pre-skeleton area filter (default 10)
%            Method-specific fields (see each algorithm's help):
%              hysteresis:          .threshHigh, .threshLow
%              adaptiveHysteresis:  .ratio, .kSigma
%              nms:                 .sigma, .threshold
%              sauvola:             .windowSize, .k, .r
%              phansalkar:          .windowSize, .k, .r, .p, .q
%              watershed:           .threshold, .threshLow, .smoothSigma, .hMinima
%              phaseCong:           .nscale, .norient, .minWL, .mult, .sigmaOnf, .k, .threshold
%              hessianCenterline:   .sigma, .threshold
%            Optional field used by 'nms':
%              .orf  – 2-D single orientation field (radians), same size as im.
%                      Encodes the local ridge direction (along the tubule).
%                      Pass [] or omit for automatic Hessian-based estimation.
%
% OUTPUTS
%   skel – logical, single-pixel-wide skeleton, same size as im.
%   bw   – logical, pre-skeleton binary mask (before bwskel), same size.
%
% DEPENDENCIES (must be on MATLAB path)
%   skeletonPostProcess        (Skeletonization_sandbox/utils/)
%   hysteresisSkeletonize      (Skeletonization_sandbox/src/)
%   adaptiveHysteresisSkeletonize
%   nmsSkeletonize
%   sauvolaSkeletonize
%   phansalkarSkeletonize
%   watershedSkeletonize
%   phaseCongSkeletonize       (also requires phasecong3 on path)
%   hessianCenterlineSkeletonize

method = getf(params, 'method', 'none');
orf    = getf(params, 'orf',    []);

% --- Route to method -----------------------------------------------------
switch method

    case 'none'
        bw = false(size(im));

    case 'hysteresis'
        bw = hysteresisSkeletonize(im, params);

    case 'adaptiveHysteresis'
        bw = adaptiveHysteresisSkeletonize(im, params);

    case 'nms'
        bw = nmsSkeletonize(im, orf, params);

    case 'sauvola'
        bw = sauvolaSkeletonize(im, params);

    case 'phansalkar'
        bw = phansalkarSkeletonize(im, params);

    case 'watershed'
        bw = watershedSkeletonize(im, params);

    case 'phaseCong'
        bw = phaseCongSkeletonize(im, params);

    case 'hessianCenterline'
        bw = hessianCenterlineSkeletonize(im, params);

    otherwise
        error('skeletonSegment: unknown method ''%s''', method);
end

% --- Shared post-processing: area filter + bwskel -----------------------
skel = skeletonPostProcess(bw, params);
end

% ---- local helper -------------------------------------------------------
function v = getf(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end
