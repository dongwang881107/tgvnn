# tgvnn: Dynamic MR Image Reconstruction using TGV and Low Rank Decomposition

tgvnn is an open-source CUDA-based GPU code to perform compressed sensing based dynamic MRI reconstruction as fast as possible with validated accuracy. The propoed model decomposes the dynamic MR image as a sparse component L and a low rank component S, which are contrainted by second order TGV and nuclear norm respectively. The model is 

<a href="https://www.codecogs.com/eqnedit.php?latex=\LARGE&space;\min_{L,S}&space;\frac{1}{2}\|A(L&plus;S)-B\|_{\mathrm{F}}^2&plus;\mathrm{TGV}^2_\alpha(S)&plus;\beta\|L\|_*," target="_blank"><img src="https://latex.codecogs.com/png.latex?\LARGE&space;\min_{L,S}&space;\frac{1}{2}\|A(L&plus;S)-B\|_{\mathrm{F}}^2&plus;\mathrm{TGV}^2_\alpha(S)&plus;\beta\|L\|_*," title="\LARGE \min_{L,S} \frac{1}{2}\|A(L+S)-B\|_{\mathrm{F}}^2+\mathrm{TGV}^2_\alpha(S)+\beta\|L\|_*," /></a>

where A is the sampling matrix and B is the undersampled data. tgvnn implements the proposed model using Primal-Dual algorithm.

## Dependencies

This code has been tested on 
- Ubuntu 18.04 & macOS High Sierra
- CUDA 10.0 & MATLAB 2018a/2017b

## Installing

To install, run the following commands in the `src/` subdirectory:
- `make`
- `sudo make install`


## Usage
```
Usage: runtgv [OPTION] <img.ra> <mask.ra>
	-o, --output <output.ra>	 recon output RA file
	-i, --iter n			 iteration number
	-a, --alpha n			 parameter for TGV
	-b, --beta n			 parameter for nuclear norm
	-s, --sigma n			 dual stepsize
	-t, --tau n			 primal stepsize
	-m, --mu n			 temporal stepsize
	-G, --gridsize n		 set GPU gridsize
	-B, --blocksize n		 set GPU blocksize
	-h				 show this help
```
## Test
```
Reading input ...
Input image: ../data/pincat.ra
Input mask:  ../data/mask.ra
gridsize:  128
blocksize: 128
alpha: 0.0040
beta:  0.5000
sigma: 0.2500
tau:   0.2500
iter:  500
rows = 128, cols = 128, ndyn = 50, N = 819200

Running recon ...
The SER of zerofill: 20.53 dB
The SER of recon:    32.74 dB
Elapsed time: 7.65 s

Saving output ...
Output: ../result/recon.ra
```

## Funding Sources
This work is supported by National Natural Science Foundation of China (NO. 91330101 and NO. 11531005).
