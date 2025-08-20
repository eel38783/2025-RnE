function psi0 = V_RSU(d) 
% V2I Mode Channel SNR

Pt = 20; % (dBm)
PL = 20*log10((4*pi*d)/lambda); % (dB)
g = Pt/PL;

B = 2.5*10^11; %300(GHz)
lambda = 10^(-2); %10 (mm)
N0= 20; 

psi0 = g/(B*N0);
end

