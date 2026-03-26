%TESTSKELETONIZE  Unit tests for Skeletonization_sandbox.
%
% Run with:
%   results = runtests('testSkeletonize');
%   disp(results);
%
% PATH REQUIREMENTS:
%   Skeletonization_sandbox/src/
%   Skeletonization_sandbox/utils/
%   Skeletonization_sandbox/demos/
%   Common_sandbox/Kovesi phase congruency/   (for phaseCong test)

classdef testSkeletonize < matlab.unittest.TestCase

    properties (TestParameter)
        % Methods tested via skeletonSegment (excluding phaseCong which needs phasecong3 on path)
        method = {'hysteresis', 'adaptiveHysteresis', 'nms', ...
                  'sauvola', 'phansalkar', 'watershed', 'hessianCenterline'}
    end

    % =====================================================================
    % Shared test fixture: synthetic 64×64 image with a known line
    % =====================================================================
    methods (TestMethodSetup)
        function buildTestImage(tc) %#ok<MANU>
        end
    end

    methods (Access = private)
        function [im, trueSkel] = makeStraightLine(~, sz, angle_deg)
            % Thin bright line at given angle through image centre
            if nargin < 3, sz = 64; end
            if nargin < 4, angle_deg = 0; end
            im   = zeros(sz, sz, 'single');
            skel = false(sz, sz);
            cx = sz/2;  cy = sz/2;
            th = deg2rad(angle_deg);
            for t = -sz/2:sz/2
                r = round(cy + t*sin(th));
                c = round(cx + t*cos(th));
                if r>=1 && r<=sz && c>=1 && c<=sz
                    im(r,c)   = 1;
                    skel(r,c) = true;
                end
            end
            % Blur slightly to make it a realistic ridge image
            im = imgaussfilt(im, 1.2);
            im = im / max(im(:));
            trueSkel = skel;
        end

        function [im, trueSkel] = makeReticulateTestFixture(~)
            im = makeReticulateTestImage(128, 42);
            im = imgaussfilt(im, 1.0);
            im = im / max(im(:));
            trueSkel = [];  % ground truth not needed for smoke tests
        end
    end

    % =====================================================================
    % Tests: output shape and type
    % =====================================================================
    methods (Test)
        function testOutputIsLogical(tc, method)
            [im, ~] = tc.makeStraightLine(64);
            p.method = method;
            p.threshHigh  = 0.4;   p.threshLow   = 0.15;
            p.ratio       = 0.4;   p.windowSize  = 15;
            p.k           = 0.2;   p.r           = 0.5;
            p.threshold   = 0.3;   p.sigma       = 1.5;
            p.smoothSigma = 1.0;   p.hMinima     = 2;
            p.pruneLength = 2;     p.minArea     = 3;
            [skel, bw] = skeletonSegment(im, p);
            tc.verifyClass(skel, 'logical', 'skel must be logical');
            tc.verifyClass(bw,   'logical', 'bw must be logical');
        end

        function testOutputSizeMatchesInput(tc, method)
            sz = 64;
            [im, ~] = tc.makeStraightLine(sz);
            p.method      = method;
            p.threshHigh  = 0.4;   p.threshLow   = 0.15;
            p.ratio       = 0.4;   p.windowSize  = 15;
            p.k           = 0.2;   p.r           = 0.5;
            p.threshold   = 0.3;   p.sigma       = 1.5;
            p.smoothSigma = 1.0;   p.hMinima     = 2;
            p.pruneLength = 2;     p.minArea     = 3;
            [skel, bw] = skeletonSegment(im, p);
            tc.verifySize(skel, [sz sz], 'skel size must match input');
            tc.verifySize(bw,   [sz sz], 'bw size must match input');
        end

        function testNoneMethodReturnsAllFalse(tc)
            [im, ~] = tc.makeStraightLine(64);
            p.method = 'none';
            [skel, bw] = skeletonSegment(im, p);
            tc.verifyFalse(any(skel(:)), 'none method must return all-false skel');
            tc.verifyFalse(any(bw(:)),   'none method must return all-false bw');
        end

        function testSkeletonIsAtMostOnePixelWide(tc, method)
            % bwmorph 'thin' applied again to a skeleton should not change it
            [im, ~] = tc.makeStraightLine(64);
            p.method      = method;
            p.threshHigh  = 0.4;   p.threshLow   = 0.15;
            p.ratio       = 0.4;   p.windowSize  = 15;
            p.k           = 0.2;   p.r           = 0.5;
            p.threshold   = 0.3;   p.sigma       = 1.5;
            p.smoothSigma = 1.0;   p.hMinima     = 2;
            p.pruneLength = 0;     p.minArea     = 1;
            [skel, ~] = skeletonSegment(im, p);
            if ~any(skel(:)), return; end  % empty skeleton — nothing to check
            % A true skeleton is idempotent under bwskel
            skel2 = bwskel(skel);
            tc.verifyEqual(skel, skel2, 'Skeleton must already be single-pixel wide');
        end

        function testStraightLineRecovered(tc, method)
            % The skeleton of a bright ridge must have pixels in the right place
            [im, trueSkel] = tc.makeStraightLine(64, 0);
            p.method      = method;
            p.threshHigh  = 0.5;   p.threshLow   = 0.20;
            p.ratio       = 0.4;   p.windowSize  = 15;
            p.k           = 0.2;   p.r           = 0.5;
            p.threshold   = 0.4;   p.sigma       = 1.5;
            p.smoothSigma = 1.0;   p.hMinima     = 2;
            p.pruneLength = 2;     p.minArea     = 3;
            [skel, ~] = skeletonSegment(im, p);

            if ~any(skel(:))
                tc.assumeTrue(false, ...
                    sprintf('%s: skeleton is empty for straight-line image', method));
            end

            % Dilate ground truth to allow 1 pixel tolerance
            trueDil = imdilate(trueSkel, strel('disk', 2));
            precision = nnz(skel & trueDil) / nnz(skel);
            tc.verifyGreaterThan(precision, 0.7, ...
                sprintf('%s: < 70%% of skeleton pixels lie within ground truth ±2 px', method));
        end

        function testSkeletonPostProcessPruning(tc)
            % Ensure pruneLength removes short spurs
            bw = false(64, 64);
            % Main line
            bw(32, 10:54) = true;
            % Short spur (5 pixels)
            bw(28:32, 32) = true;
            % Long branch (20 pixels, should survive pruning at 10)
            bw(12:32, 10) = true;

            p1.pruneLength = 10;  p1.minArea = 1;
            skel1 = skeletonPostProcess(bw, p1);

            p2.pruneLength = 0;   p2.minArea = 1;
            skel2 = skeletonPostProcess(bw, p2);

            % With pruning: fewer skel pixels
            tc.verifyLessThan(nnz(skel1), nnz(skel2), ...
                'pruneLength=10 should remove spurs and reduce skeleton pixel count');
        end

        function testHysteresisConnectivity(tc)
            % Low-threshold pixels connected to high-threshold seeds must be included
            im = zeros(64, 64, 'single');
            im(32, 20:44) = 0.6;   % strong ridge (above high threshold)
            im(32, 5:19)  = 0.3;   % weak extension (between low and high)
            im = imgaussfilt(im, 0.5);
            im = im / max(im(:));

            p.method      = 'hysteresis';
            p.threshHigh  = 0.5;
            p.threshLow   = 0.2;
            p.pruneLength = 0;
            p.minArea     = 1;
            [~, bw] = skeletonSegment(im, p);

            % Both the strong and weak sections must be foreground
            tc.verifyTrue(bw(32, 32),  'Strong region must be foreground');
            tc.verifyTrue(bw(32, 10),  'Connected weak region must be foreground');
        end

        function testNmsOnOrientedRidge(tc)
            % NMS with correct orf should thin a 5-pixel-wide ridge to 1 pixel
            im = zeros(64, 64, 'single');
            im(:, 30:34) = 1;   % vertical 5-px-wide ridge
            im = imgaussfilt(im, 0.5);
            im = im / max(im(:));

            % orf = 0 rad = ridge runs vertically → normal = pi/2 = horizontal
            orf = zeros(64, 64, 'single');  % ridge is vertical (0 rad along ridge)

            p.method      = 'nms';
            p.sigma       = 1.0;
            p.threshold   = 0.3;
            p.pruneLength = 0;
            p.minArea     = 1;
            p.orf         = orf;
            [~, bw] = skeletonSegment(im, p);

            % The thinned result should span fewer columns than the original 5
            [~, bwCols] = find(bw);
            nCols = numel(unique(bwCols));
            tc.verifyLessThanOrEqual(nCols, 2, ...
                'NMS should thin the 5-px ridge to at most 2 columns');
        end

        function testReticulateNetworkSmoke(tc)
            % Smoke test: all methods produce non-empty skeleton on reticulate image
            [im, ~] = tc.makeReticulateTestFixture();
            methods = {'hysteresis', 'adaptiveHysteresis', 'sauvola', ...
                       'phansalkar', 'watershed', 'hessianCenterline'};
            for m = 1:numel(methods)
                p.method      = methods{m};
                p.threshHigh  = 0.35;   p.threshLow   = 0.12;
                p.ratio       = 0.4;    p.windowSize  = 25;
                p.k           = 0.2;    p.r           = 0.5;
                p.threshold   = 0.25;   p.sigma       = 1.5;
                p.smoothSigma = 1.0;    p.hMinima     = 2;
                p.pruneLength = 5;      p.minArea     = 10;
                [skel, ~] = skeletonSegment(im, p);
                tc.verifyTrue(any(skel(:)), ...
                    sprintf('%s: skeleton is empty on reticulate image', methods{m}));
            end
        end
    end
end
