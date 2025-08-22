function psi0 = V_RSU(d0)
% V2I Mode Channel SNR

    Pt = 20;  % (dBm)
    lambda = 0.01;  %10 (mm)
    PL = 20 * log10((4*pi*d0)/lambda);
    g = Pt/PL;

    B = 2.5*10^(11);  %300(GHz)
    N0 = 20;

    psi0 = g/(B*N0);
end
