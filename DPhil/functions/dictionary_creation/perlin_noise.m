function [perlinNoise, noises] = perlin_noise(vecLength, varargin)
%PERLIN_NOISE Function to generate 1-Dimensional Perlin noise.
%
% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
%
%
% Uses a uniform random number generator and an interpolator.
% Written using http://freespace.virgin.net/hugo.elias/models/m_perlin.htm
% as a guide.
%
% PERLINNOISE = PERLIN_NOISE(VECLENGTH, VARARGIN)
%
% VECLENGTH
%   scalar to set how long the resulting Perlin noise vector should be.
%
% VARARGIN{1}
%   scalar to set the lower limit when scaling the resulting Perlin noise
%   vector.
%
% VARARGIN{2}
%   scalar to set the upper limit when scaling the resulting Perlin noise
%   vector.
%
% VARARGIN{3}
%   scalar to set the 'persistence' used when adding the noise vectors.
%
%

%%
% set the initial perlin noise vector to a vector of zeros.
perlinNoise = zeros(1,vecLength);
% we will generate twice as many points as we need, to avoid problems with
% differences in the lengths of the interpolated noise vectors.
numGeneratedPts = 2*vecLength;

if nargin>3
    p = varargin{3};
else
    p = 0.5;
end

%set the number of noise vectors to add
if nargin>4
    numOctaves = varargin{4};
else
    numOctaves = 5;
end
 numOctaves = 5;
maxperiod = 2^(numOctaves-1);

noises = zeros(numOctaves,numGeneratedPts*2);

%%
for  n = 0:numOctaves-1
    %generate a different uniformly distributed noise vector for each
    %iteration.
    % values are from 0 to 1.
    rng(n)
    noise = rand(1,numGeneratedPts/maxperiod);
    
    %interpolate the noise vector of this iteration.
    %period of the interpolation increases with n.
    period = 2^n;
    noise = interp1(1:(numGeneratedPts/maxperiod),noise,1:(1/(period*maxperiod)):(numGeneratedPts/(period*maxperiod)),'spline');

    % add together the outputs from different noise functions.
    %amplitude a of the added noise, depends on persitence p.
    a = p^n;
    perlinNoise = perlinNoise + noise(1:vecLength)*a;
    
   
end


%%
if nargin>=2
    lower = varargin{1};
    upper = varargin{2};
else
    lower = 0;
    upper = 1;
end


perlinNoise = (upper - lower)*(perlinNoise/max(perlinNoise)) + lower;
end

