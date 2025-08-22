function psi1 = V2V_SNR(veh_pos,K)
    psi1 = zeros(length(K));
    gamma_vals = linspace(-70, 70, 2001);


    for i=1:15
        for j=1:15
            d_tr = norm(veh_pos(i) - veh_pos(j));
            pdf_vals = V_2V(d_tr);
            cdf_vals = cumtrapz(gamma_vals, pdf_vals);
            cdf_vals = cdf_vals / cdf_vals(end);

              % 2. Uniform(0,1) 확률 하나 생성
            u = rand();

            % 3. u에 가장 가까운 cdf 값의 인덱스 찾기
            [~, idx] = min(abs(cdf_vals - u));

            % 4. 해당 인덱스의 gamma 값을 반환
            psi1(i,j) = gamma_vals(idx);
        end
    end
end
