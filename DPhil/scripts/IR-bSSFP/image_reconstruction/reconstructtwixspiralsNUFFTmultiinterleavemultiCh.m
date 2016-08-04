%% Reconstruct Raw Data using NUFFT - Spiral Sampled
%
%
clear all
%% 1. Set up paths
%Are you working on jalapeno00 or locally?
% workingdir = '/home/fs0/jallen/Documents/MATLAB/DPhil';
workingdir = '/Users/jallen/Documents/MATLAB/DPhil';
addpath(genpath(workingdir)); % sometimes causes MATLAB to freeze
%
savingdir = '/Users/jallen/Documents/DPhil';
addpath(genpath(savingdir));

% If working on jalapeno00, uncomment the following lines:
% addpath(genpath('/Applications/fsl/'))
% addpath(genpath('/usr/local/fsl/bin'))
% addpath(genpath('/opt/fmrib/fsl/etc/matlab'))
addpath('/Users/jallen/Documents/MATLAB/irt',...
    '/Users/jallen/Documents/MATLAB/vdspiral')
addpath(genpath('/Users/jallen/Documents/MATLAB/DPhil'),...
    genpath('/Users/jallen/Documents/MATLAB/mapVBVD_20150918'))
run('/Users/jallen/Documents/MATLAB/irt/setup')

%% Spiral Reconstruction
disp('Prepare Spiral Trajectories...')
run('make_spirals')
disp('Prepare Spiral Trajectories...Finished')

ADCshifts = 0:10:110;

summed_image = zeros(256);

for ADCshift = ADCshifts(1)
    figure;
    for nChannel = [20]
        
        %% Load raw data
        disp('Load Raw Data...')
        filename = (['/Users/jallen/Documents/DPhil/data/TWIX/20160415/meas_MID7',num2str(28+(ADCshifts(1+(ADCshift*0.1))*0.1)),'_JA_IR_bSSFP_fp_tr10_fa10_noinv_shift',num2str(ADCshifts(1+(ADCshift*0.1))),'_FID',num2str(9497+(ADCshifts(1+(ADCshift*0.1))*0.1))]);
        twix_obj = mapVBVD(filename);
        twix_obj.image.flagRemoveOS = 1; %remove factor of 2 oversampling along each trajectory
        image_data = twix_obj.image{''};
        disp('Load Raw Data...Finished')
        
        %% Prepare Data for NUFFT
        disp('Prepare Data Vectors for NUFFT...')
        numInterleaves = 48;
        for n = 1:numInterleaves
            leaf_values((n-1)*nSamplePts+1:n*nSamplePts) = cast(squeeze(image_data(:,nChannel,n)),'double'); %measured k-space data values
            samplelocations(1,(n-1)*nSamplePts+1:n*nSamplePts) = spiralkspace(n,1:nSamplePts); %desired k-space sample positions
        end
        % magnitude of radii must be < 5
        samplelocations = complex(0.5*real(samplelocations)/max(real(samplelocations)),0.5*imag(samplelocations)/max(imag(samplelocations))) ;

        radii = sqrt((real(samplelocations)).^2 + (imag(samplelocations)).^2);
        %
        MLfile = [real(samplelocations); imag(samplelocations)]';
        MLvec_of_samples = leaf_values;
        fft_values=[MLfile MLvec_of_samples'];
        disp('Prepare Data Vectors for NUFFT...Finished')
        %
        
        disp('Perform NUFFT...')
        
        %% reconstruction, based on Lior's use of Fessler's code
        N = 256; %dimension of image
        
        %%| Create a "image geometry" structure that describes the sampling
        %| characteristics of a single 2d image.
        ig=image_geom('nx',N,'ny',N,'dx',1,'offsets','dsp');
        % define a mask
        mask=ig.circ(1+ig.nx/2, 1+ig.ny/2) > 0;
       % mask= ones(N)>0;
        ig.mask=mask;
        mask_locs=find(mask);
        
        % set NUFFT parameters
        J = [6 6];
        N2=[N N];
        nufft_args = {N2, J, 2*N2, N2/2, 'table', 2^10, 'minmax:kb'};
        %Construct Gmri object for MR image reconstruction, or for other problems
        %| that require reconstruction of a function from nonuniform Fourier samples.
        %| This object can also model relaxation and/or off-resonance effects for MRI.
        Gm = Gmri(fliplr(fft_values(:,1:2)), ig.mask, ...
            'fov', ig.fov, 'basis', {'rect'}, 'nufft', nufft_args);
        
        % Build roughness penalty regularization "object" based on C = Cdiff() object,
        % for regularized solutions to inverse problems.
        beta = 2^-7 * size(fft_values,1);
        Rn = Robject(ig.mask, 'type_denom', 'matlab', 'potential', 'hyper3', 'beta', 2^2*beta, 'delta', 0.3);
        %| penalized weighted least squares (PWLS)
        %| with convex non-quadratic regularization,
        %| minimized via preconditioned conjugate gradient algorithm.
        niter = 40;
        %reconstruct image
        xh = pwls_pcg1(zeros(Gm.dim(2),1), Gm, 1, fft_values(:,3), Rn, 'niter', niter);
        
        reconstructed_image=zeros(N);
        reconstructed_image(mask_locs)=xh;
        
        %%
        disp('Perform NUFFT...Finished')
        dt = datetime('now');
        abs_reconstructed_image = abs(reconstructed_image);
        summed_image = summed_image + abs_reconstructed_image.^2;
        combined_image = sqrt(summed_image);
        % save image
    %    save([savingdir,'/MAT-files/images/reconstructed_imageMI_Ch',num2str(nChannel),'_ADCshift',num2str(ADCshift),'',datestr(dt),'.mat'],'abs_reconstructed_image')
       % subplot(4,8,nChannel)
        imagesc(abs_reconstructed_image);
        title(['ADCshift:',num2str(ADCshift),' Ch',num2str(nChannel)])
     %   save([savingdir,'/MAT-files/images/reconstructed_imageMI_Ch',num2str(nChannel),'_ADCshift',num2str(ADCshift),'',datestr(dt),'.mat'],'abs_reconstructed_image')
    end
    matlab2tikz([savingdir,'/MAT-files/images/reconstructed_imageMI_Ch',num2str(nChannel),'_ADCshift',num2str(ADCshift),'',datestr(dt),'.tex'])
    
end


