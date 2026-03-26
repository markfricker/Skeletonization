function I = makeReticulateTestImage(sz, seed)
%MAKERETICULATETESTIMAGE  Synthetic ER-like reticulate tubular network image.
%
% USAGE
%   I = makeReticulateTestImage()           % 256×256, fixed seed
%   I = makeReticulateTestImage(sz)         % sz×sz image
%   I = makeReticulateTestImage(sz, seed)
%
% DESCRIPTION
%   Generates a single-precision fluorescence image of a reticulate
%   tubular network resembling the endoplasmic reticulum (ER) in confocal
%   microscopy.  The network topology is derived from a perturbed grid of
%   nodes connected by Delaunay triangulation edges, producing the polygon-
%   mesh appearance characteristic of the ER.
%
%   Network properties (matched to ER confocal data):
%     Tubule sigma    : 1.2 px (≈ 2.8 px FWHM after PSF)
%     Junction nodes  : brighter than tubule midpoints (fluorophore pooling)
%     PSF blurring    : Gaussian sigma 0.8 px (confocal diffraction limit)
%     Shot noise      : Poisson, 60 peak photons
%     Readout noise   : Gaussian σ = 0.015 (normalised)
%
% OUTPUT
%   I  – single-precision image in [0, 1], size [sz × sz]
%
% GROUND TRUTH (accessible via outputs from the nested build)
%   The true network skeleton is a binary image of the Delaunay edges.
%   To obtain it, run makeReticulateTestImage with two outputs:
%     [I, trueSkel] = makeReticulateTestImage(...)    (future extension)

if nargin < 1, sz   = 256; end
if nargin < 2, seed = 42;  end

rng(seed);
canvas = zeros(sz, sz, 'single');

% --- Network topology: perturbed regular grid of nodes -------------------
gridN  = 7;          % number of grid points per axis
margin = round(sz * 0.10);
positions = linspace(margin, sz - margin, gridN);
[gx, gy]  = meshgrid(positions, positions);

% Random perturbation to break regularity (± 40% of grid spacing)
spacing = sz / (gridN + 1);
gx = gx + (rand(size(gx)) - 0.5) * spacing * 0.8;
gy = gy + (rand(size(gy)) - 0.5) * spacing * 0.8;

% Clamp to image bounds with margin
gx = max(margin, min(sz - margin, gx));
gy = max(margin, min(sz - margin, gy));

nodes = [gx(:), gy(:)];   % [nNodes × 2]:  col, row

% Delaunay triangulation gives the reticulate edge set
DT       = delaunayTriangulation(nodes(:,1), nodes(:,2));
edgeList = DT.edges;   % [nEdges × 2] indices into nodes

% --- Tubule rendering parameters -----------------------------------------
tubuleSigma   = 1.2;   % Gaussian cross-section sigma (pixels)
baseIntensity = 0.8;   % peak tubule fluorescence (normalised)

% --- Render each edge as a series of Gaussian spots ---------------------
for e = 1:size(edgeList, 1)
    n1 = nodes(edgeList(e, 1), :);   % [col, row]
    n2 = nodes(edgeList(e, 2), :);

    % Vary intensity slightly along each tubule
    intVar = baseIntensity * (0.6 + 0.4 * rand);

    len  = norm(n2 - n1);
    nPts = ceil(len) + 1;

    for k = 0:nPts
        t  = k / nPts;
        cx = n1(1) + t * (n2(1) - n1(1));   % col (x)
        cy = n1(2) + t * (n2(2) - n1(2));   % row (y)
        canvas = addSpot(canvas, cx, cy, tubuleSigma, intVar, sz);
    end
end

% --- Junction nodes are brighter (fluorophore pooling) -------------------
nodeIntensity = baseIntensity * 1.4;
nodeSigma     = tubuleSigma * 1.3;
for n = 1:size(nodes, 1)
    canvas = addSpot(canvas, nodes(n, 1), nodes(n, 2), nodeSigma, nodeIntensity, sz);
end

% --- Gaussian PSF blur (confocal diffraction limit) ----------------------
canvas = imgaussfilt(canvas, 0.8);

% --- Poisson shot noise --------------------------------------------------
peakPhotons = 60;
canvas = single(poissrnd(double(canvas) * peakPhotons)) / peakPhotons;

% --- Gaussian readout noise ----------------------------------------------
readNoise = 0.015;
canvas = canvas + readNoise * randn(size(canvas), 'single');
canvas = max(0, min(1, canvas));

I = canvas;
end

% ---- local helper: render a Gaussian spot at (cx, cy) ------------------
function canvas = addSpot(canvas, cx, cy, sigma, intensity, sz)
r  = ceil(4 * sigma);
x1 = max(1, floor(cx) - r);  x2 = min(sz, ceil(cx) + r);
y1 = max(1, floor(cy) - r);  y2 = min(sz, ceil(cy) + r);
if x1 > x2 || y1 > y2, return; end

[gx, gy] = meshgrid(x1:x2, y1:y2);
g = intensity * exp(-((gx - cx).^2 + (gy - cy).^2) / (2 * sigma^2));
canvas(y1:y2, x1:x2) = max(canvas(y1:y2, x1:x2), single(g));
end
