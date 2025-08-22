
function T = Opt_func(K, rho, Ds, C, psi, t_tole, t_hold)


    cvx_begin
        variable epsilon(15,15,2)
        Time = 0;
        for i=1:(K(2)+K(3))
            for j=1:(K(2)+K(3)+K(4))
                if j==i; continue; end
                    for k=1:2 
                        Time = Time + epsilon(i,j,k)*(C(i)/rho(j) + (2*Ds(i))/(psi(i,j,k))) + epsilon(i,i,1)*(C(i)/rho(j));
                    end
            end
        end

        for i=1:K(1) 
            Time = Time + C(i)/rho(i);
        end

        minimize real(Time)
        subject to
            epsilon >= 0
            epsilon <= 1
            %for i=1:(K(2)+K(3)) 
                %for j=1:(K(2)+K(3)+K(4))
                    %for k=1:2 
                        %epsilon(i,j,k) >= 0;
                        %epsilon(i,j,k) <= 1;
                    %end
                %end
            %end
            % relaxation factor
            for i=1:K(1)
                epsilon(i,1,1) * (C(i)/rho(i)+(2*Ds(i)/psi(i,1,1))) ...
                     <= min(t_tole(i),t_hold(i,1,1));
            end
        
            for i=1:(K(2)+K(3)) 
                for j=2:(K(2)+K(3)+K(4)) 
                    real(epsilon(i,j,1) * (C(i)/rho(j)+(2*Ds(i)/psi(i,j,1)))) ...
                        <= min(t_tole(i), t_hold(i,j,1));
                end
            end

            for i=1:(K(2)+K(3))
                for j=2:(K(2)+K(3)+K(4)) 
                    epsilon(i,j,2)*(C(i)/rho(j)+(2*Ds(i)/psi(i,j,2))) ...
                        <= min(t_tole(i), t_hold(i,j,2));
                end
            end
            % time constraints
            e = 0;
            for j=1:(K(2)+K(3)+K(4)) 
                for k=1:2 
                    for i=(K(1)+K(2)+K(3)) 
                        e = e + epsilon(i,j,k);
                    end
                    e == 1;
                end
            end
            % to ensure 1:1 match

            for i=1:(K(2)+K(3)) 
                for j=2:(K(2)+K(3)+K(4)) 
                    epsilon(i,j,1) + eps(i,j,2) == 1; 
                end 
            end
        % to ensure V2V Mode is either on the ground or through UAV-RIS
    cvx_end

    T = Time;
end