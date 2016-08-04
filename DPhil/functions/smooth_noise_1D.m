function [ smoothedNoise ] = smooth_noise_1D( x )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

smoothedNoise = random(x)/2 + random(x-1)/4 + random(x+1)/4;


end

