
function [psi] = Rate(veh_pos)
    psi = zeros(15, 15, 2);
    
    for i=1:15 
        d_i0 = norm(veh_pos(i,1:3) - [0,0,10]); %RSU Height <-- 10m
        psi(i,1,1) = log10(1 + V_RSU(d_i0)); 
        psi(i,1,2) = log10(1 + V_RSU(d_i0)); 
        for j=2:15
            psi2 = V_RIS(4,4,theta,UAV,veh_pos);
            psi(i,j,1) = log10(1 + V2V_SNR(veh_pos,i,j));
            psi(i,j,2) = log10(1 + psi2(i,j));
        end
    end 
                  
end 

function y = distance(m, n, veh_pos)
    d_mn = veh_pos(m, 1:2) - veh_pos(n, 1:2);
    y = sqrt((d_mn(1)^2 + (d_mn(2))^2);
end