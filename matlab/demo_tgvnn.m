%% 3DTGVNN demo on dynamic data

clear;
close all;

%% Load Data
load('../data/pincat.mat');
img = new;
img_name = 'pincat';

% load('../data/perfusion.mat');
% img = x;
% img_name = 'perfusion';

% load('../data/breast.mat');
% img = abs(data(:,1:128,:));
% img_name = 'breast';

img = abs(img);
img = img/max(img(:));
[m,n,d] = size(img);

%% Set up parameters
reduction = 5e-1;
ratio = sqrt(2);
alpha = 4e-3;
beta = 5e-1;
sigma = 0.25;
tau = 0.25;
mu = 1;
iter = 500;

%% Generate samplying Operator
rng(5089);

% Psuedo-radial sampling
numSpokes = 32;
mask = strucrand(m,n,d,numSpokes);

% % Cartesian sampling
% accel = 10;
% width = 8;
% mask_seed = 1000;
% mask = vuCSRandMaskCadence([n,m],accel,d,width,false,mask_seed);
% mask = repmat(mask,[n,1,1]);
% mask = reshape(mask,[m,n,d]);

fprintf('Sampling factor is %.2f%%\n',numel(find(mask==1))/m/n/d*100);
[A,AT] = FFTOperator(mask);

imgb = A(img);
imgz = AT(imgb);

%% Put parameter into input
input.A = A;
input.AT = AT;
input.B = imgb;
input.alpha0 = ratio*alpha;
input.alpha1 = alpha;
input.alpha = alpha;
input.beta = beta;
input.sigma = sigma;
input.tau = tau;
input.mu = mu;
input.img = img;
input.iter = iter;
input.reduction = reduction;

%% Loop through mode
tic;
[imgl,imgs] = tgvnn(input);
imgr = imgl + imgs;
t = toc;
fprintf('Elapsed time: %.2f s\n',t);

% Display Results
ser_imgz = -20*log10(norm((abs(imgz(:)) - img(:)))/norm(img(:)));
ssim_imgz = ssim(abs(img),abs(imgz), ...
            [0.01 0.03], fspecial('gaussian',11,1.5),1);

ser_imgr = -20*log10(norm((abs(imgr(:)) - img(:)))/norm(img(:)));
ssim_imgr = ssim(abs(img),abs(imgr), ...
            [0.01 0.03], fspecial('gaussian',11,1.5),1);
        
fprintf('The SER  of zerofill: %.2f dB\n',ser_imgz);
fprintf('The SER  of recon:    %.2f dB\n',ser_imgr);
fprintf('The SSIM of zerofill: %.2f dB\n',ssim_imgz);
fprintf('The SSIM of recon:    %.2f dB\n',ssim_imgr);
        
% Save results
out_path = strcat('../result/',img_name,'_tgvnn_',num2str(numSpokes));
save(out_path,'imgl','imgs','ser_imgr','ssim_imgr','t');
