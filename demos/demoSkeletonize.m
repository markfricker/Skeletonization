%DEMOSKELETONIZE  Compare skeleton extraction methods on a synthetic ER network.
%
% Demonstrates the full pipeline:
%   1. Synthetic reticulate network image (makeSkeletonTestImage)
%   2. Phase congruency enhancement (phasecong3 — Kovesi)
%   3. All skeleton methods compared side-by-side
%
% PATH REQUIREMENTS (add to MATLAB path before running):
%   Skeletonization_sandbox/src/
%   Skeletonization_sandbox/utils/
%   Skeletonization_sandbox/demos/
%   Common_sandbox/Kovesi phase congruency/   (for phasecong3)
%
% The demo is self-contained and does not require the AnalyzER GUI.

clear; close all; clc;

% =========================================================================
% 1. Synthetic test image
% =========================================================================
fprintf('Generating reticulate test image...\n');
I = makeSkeletonTestImage(256, 42);

% =========================================================================
% 2. Phase congruency enhancement (default first-pass enhance)
% =========================================================================
fprintf('Running phase congruency enhancement...\n');

% Default phasecong3 parameters tuned for thin tubular structures:
%   nscale=4, norient=6, minWL=3, mult=2.1, sigmaOnf=0.55, k=2
[M, ~, ~, ft, ~, ~, T] = phasecong3(double(I), 4, 6, 3, 2.1, 0.55, 2);

% M = maximum moment (ridge + edge strength)
% ft = local phase angle: pi/2 = bright line, 0 = step edge
% Retain only bright-line features (ft > 0 selects ridges over edges)
ridgeStrength = single(M) .* single(ft > 0);
ridgeStrength = ridgeStrength / (max(ridgeStrength(:)) + eps);

% =========================================================================
% 3. Skeleton methods
% =========================================================================
fprintf('Running skeleton methods...\n');

% Shared post-processing parameters
postP.pruneLength = 8;
postP.minArea     = 15;

% --- hysteresis -----------------------------------------------------------
p_hyst            = postP;
p_hyst.threshHigh = 0.40;
p_hyst.threshLow  = 0.15;
p_hyst.method     = 'hysteresis';
[skel_hyst, bw_hyst] = skeletonSegment(ridgeStrength, p_hyst);

% --- adaptiveHysteresis ---------------------------------------------------
p_ahyst        = postP;
p_ahyst.ratio  = 0.4;
p_ahyst.method = 'adaptiveHysteresis';
[skel_ahyst, ~] = skeletonSegment(ridgeStrength, p_ahyst);

% --- nms (orientation from image) ----------------------------------------
p_nms           = postP;
p_nms.sigma     = 1.5;
p_nms.threshold = 0;    % auto Otsu
p_nms.method    = 'nms';
[skel_nms, ~] = skeletonSegment(ridgeStrength, p_nms);

% --- sauvola --------------------------------------------------------------
p_sauv            = postP;
p_sauv.windowSize = 31;
p_sauv.k          = 0.15;
p_sauv.r          = 0.5;
p_sauv.method     = 'sauvola';
[skel_sauv, ~] = skeletonSegment(ridgeStrength, p_sauv);

% --- phansalkar -----------------------------------------------------------
p_phan            = postP;
p_phan.windowSize = 31;
p_phan.k          = 0.25;
p_phan.r          = 0.5;
p_phan.method     = 'phansalkar';
[skel_phan, ~] = skeletonSegment(ridgeStrength, p_phan);

% --- watershed ------------------------------------------------------------
p_ws              = postP;
p_ws.threshold    = 0.30;
p_ws.threshLow    = 0.10;
p_ws.smoothSigma  = 1.0;
p_ws.hMinima      = 2;
p_ws.method       = 'watershed';
[skel_ws, ~] = skeletonSegment(ridgeStrength, p_ws);

% --- phaseCong ------------------------------------------------------------
% Use the raw image here (phasecong3 is called internally)
p_pc           = postP;
p_pc.nscale    = 4;
p_pc.norient   = 6;
p_pc.minWL     = 3;
p_pc.k         = 2.0;
p_pc.threshold = 0;   % auto data-driven threshold
p_pc.method    = 'phaseCong';
[skel_pc, ~] = skeletonSegment(I, p_pc);

% --- hessianCenterline ----------------------------------------------------
p_hcl           = postP;
p_hcl.sigma     = 1.5;
p_hcl.threshold = 0;   % auto Otsu
p_hcl.method    = 'hessianCenterline';
[skel_hcl, ~] = skeletonSegment(ridgeStrength, p_hcl);

% =========================================================================
% 4. Display
% =========================================================================
fprintf('Plotting results...\n');

skels   = {skel_hyst, skel_ahyst, skel_nms, skel_sauv, ...
           skel_phan, skel_ws, skel_pc, skel_hcl};
labels  = {'hysteresis', 'adaptiveHysteresis', 'nms', 'sauvola', ...
           'phansalkar', 'watershed', 'phaseCong', 'hessianCenterline'};
nMethod = numel(skels);

figure('Name', 'Skeleton method comparison', 'Color', 'w', ...
    'Position', [50, 50, 1400, 700]);

% Row 1: raw image and enhanced ridge map
subplot(3, 5, 1);
imshow(I, []); title('Raw image', 'Interpreter', 'none');

subplot(3, 5, 2);
imshow(ridgeStrength, []); title('Phase cong. ridge map', 'Interpreter', 'none');

subplot(3, 5, 3);
imshow(bw_hyst, []); title('Hysteresis BW mask', 'Interpreter', 'none');

% Rows 2-3: skeletons overlaid on grayscale image
for m = 1:nMethod
    ax = subplot(3, 5, 4 + m);
    imshow(I, []); hold on;
    [r, c] = find(skels{m});
    plot(c, r, '.', 'Color', [1 0.3 0.3], 'MarkerSize', 1);
    title(labels{m}, 'Interpreter', 'none', 'FontSize', 8);
    hold off;
end

sgtitle('Skeleton method comparison — synthetic ER network', ...
    'FontSize', 12, 'FontWeight', 'bold');

% =========================================================================
% 5. Quantitative summary
% =========================================================================
fprintf('\n%-24s  %8s  %8s\n', 'Method', 'Skel px', 'BranchPts');
fprintf('%s\n', repmat('-', 1, 44));
for m = 1:nMethod
    nPx  = nnz(skels{m});
    nBr  = nnz(bwmorph(skels{m}, 'branchpoints'));
    fprintf('%-24s  %8d  %8d\n', labels{m}, nPx, nBr);
end

fprintf('\nDone.\n');
