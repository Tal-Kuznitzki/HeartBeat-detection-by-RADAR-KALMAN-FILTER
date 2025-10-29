function [radar_i_fitted, radar_q_fitted] = ellipseCreator(radar_i, radar_q,plotBool)
    
        radar_i_centered = radar_i - mean(radar_i);
        radar_q_centered = radar_q - mean(radar_q);
        data = [radar_i_centered radar_q_centered];
        C = cov(data);
        [V, D] = eig(C);
       % [eigvals, idx] = sort(diag(D), 'ascend');
      %  D = diag(eigvals);
      %  V = V(:, idx);
        %normalizing 
        D_inv_sq = diag(1 ./ sqrt(diag(D)));
        % by using A = VDV' 
        transformed_C = (D_inv_sq * V' * data')' ;
        radar_i_fitted = transformed_C(:,1);
        radar_q_fitted = transformed_C(:,2);

        if plotBool
            EllipsePlotter(radar_i, radar_q,radar_i_centered,radar_q_centered,radar_i_fitted,radar_q_fitted)
        end
        
end