% added in -nojvm warning suppression, alexg, oct 2012
warning('off', 'MATLAB:HandleGraphics:noJVM')

if ~exist('T2pair','var') 
    disp('ERROR: T2pair must be specified')
    return
end

if length(T2pair)~=2
    disp('ERROR: T2pair must contain 2 values')
    return
end

figure
set(gcf,'Position',[    50   679   740   409])
TE_vals = linspace(0,150,1000);
sig = zeros(length(TE_vals),2);
for iT2 = 1:2
    sig(:,iT2) = (exp(-TE_vals/T2pair(iT2)));
end
plot(TE_vals,sig(:,1),'linewidth',2); hold all
plot(TE_vals,sig(:,2),'linewidth',2); 
h1 = plot(TE_vals,abs(sig(:,1)-sig(:,2)));
set(h1,'linewidth',2)
xlabel('TE (ms)')
ylabel('Relative signal')
grid on
legend(['T2 = ' num2str(T2pair(1)) 'ms'], ['T2 = ' num2str(T2pair(2)) 'ms'], 'Contrast');
fontScale(1.4)

% added in -nojvm warning reactived, alexg, oct 2012
warning('on', 'MATLAB:HandleGraphics:noJVM')
