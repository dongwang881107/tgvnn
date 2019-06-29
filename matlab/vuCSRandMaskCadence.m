function M = vuCSRandMaskCadence (dims, accel, cadence, w, omitCorners, seed)

%%vuCSRandMaskCadence

%

% M = vuCSRandMaskTwiddle(dims,accel,w,omitCorners) returns a

% binary acquisition mask in Matlab's Fourier ordering.

%

%  dims: [nx ny nz] data size. ny is the 1st PE dimension, nx is the

%  readout dimension, and nz is the 2nd PE dimension.

%

%  accel: integer acceleration factor

%

%  cadence: period of mask repetition, e.g. 4 yields 4 masks that are

%  maximally orthogonal (?).

%

%  w: OPTIONAL radius of central fully sampled disk.  DEFAULT: 10.

%

%  omitCorners: OPTIONAL boolean option to trim off corners from mask, leading

%  to an inscribed circle pattern that could be useful with homodyning in the

%  future.  DEFAULT: false.

 

% dims = [144 144 40];

% accel = 4;

% cadence = 24;

% w = 5;

% seed = 1;

 

nx = dims(1);

ny = dims(2);

accel = round(accel);

 

 

 

if nargin < 3, cadence = 24; end

if nargin < 4, w = 5; end

if nargin < 5, omitCorners = false; end

if nargin < 6, seed = 11235; end

 

rs = RandStream.create('mt19937ar','Seed',seed);

RandStream.setGlobalStream(rs);

 

if numel(dims) == 3

  nz = dims(3);

  M = false(ny,nz);

  

  % add central window and (maybe) zero corners

  [Zfull Yfull] = meshgrid(1:nz,1:ny);

  R2 = (Yfull - ny/2).^2 + (Zfull - nz/2).^2;

  %M(ny/2-w+1:ny/2+w,nz/2-w+1:nz/2+w) = true;

  M(R2 <= w^2) = true;

  if omitCorners

    fprintf('>> corners omitted\n');

    M(R2 >= max(ny,nz)^2/4) = false;

  end

  

  % random assign each remaining point to one of the unique masks in the cadence

  idxsRem = find(~M);

  nsamples = floor(cadence / accel);

  positions = zeros(numel(idxsRem),nsamples);

  for k = 1:numel(idxsRem)

    j = randperm(cadence);

    positions(k,:) =  j(1:nsamples);

  end

  

  %

  M = repmat(M,[1 1 cadence]);

  for k  = 1:numel(idxsRem)

    [row col] = ind2sub([ny nz], idxsRem(k));

    M(row,col,positions(k,:)) = true;

  end

  

else  % 2D mask

  M = false(ny,1);

  M(round(ny/2)-w:round(ny/2)+w) = true;

  

  % random assign each remaining point to one of the unique masks in the cadence

  idxsRem = find(~M);

  nsamples = floor(cadence / accel);

  positions = zeros(numel(idxsRem),nsamples);

  for k = 1:numel(idxsRem)

    j = randperm(cadence);

    positions(k,:) =  j(1:nsamples);

  end

  

  %

  M = repmat(M,[1 cadence]);

  for k  = 1:numel(idxsRem)

    M(idxsRem(k),positions(k,:)) = true;

  end

end

 

%% plot it

% figure(1);

% subplot(121);

% imagesc(mosaic(M));

% axis image;

% subplot(122);

% imagesc(sum(M,3));

% colorbar;

% axis image;

% %colormap(gray);

 

 

% fprintf('Delivered acceleration: %g\n', numel(M) / nnz(M));