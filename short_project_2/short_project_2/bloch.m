% Bloch Equation Simulation, Excercise C-3a
% -----------------------------------------
% 

%
%	This code will plot contrast as a function
%	of TE and TR.   You may want to come back to 
%	this code soon!!

% 	Ranges of TE and TR to use.

TE = [0:20:1000];
TR = [10:100:5000];

T11 = 600;
T21 = 100;
T12 = 1000;
T22 = 150;

nte = length(TE);
ntr = length(TR);

C = zeros(nte,ntr);
%x = ones(size(TE'))*TR;
%y = TE'*ones(size(TR));



for k=1:nte
	for p=1:ntr
		if (TE(k)>TR(p))
			C(k,p)=0;
		else
			% Can modify this equation to suit! 

			C(k,p)= abs( abs(sesignal(T11,T21,TE(k),TR(p),0))-abs(sesignal(T12,T22,TE(k),TR(p),0)) ) / sqrt(TR(p));
			%C(k,p)=TE(k)+TR(p);
		end;
	end;
	tt = sprintf('%d %% complete.',round(100*k/nte));
	disp(tt);
end;

Cmx = max(size(colormap));
% Normalize C to [0,Cmx]
Cp = C-min(C(:));
Cp = Cmx*Cp/max(Cp(:));

% Make Color Scale bar on right:
Cs = [1:nte]' * Cmx/nte * ones(1,round(ntr/20));

Cp = [Cp Cs];


image(TR,TE,Cp);
xlabel('TR');
ylabel('TE');
title('Contrast-efficiency between tissue A and tissue B');
