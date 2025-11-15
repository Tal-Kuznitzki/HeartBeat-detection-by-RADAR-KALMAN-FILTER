function [radar_i_fitted, radar_q_fitted] = ellipseCreator(radar_i, radar_q,plotBool)
    
        radar_i_centered = radar_i - (max((radar_i))+min((radar_i)))/2;%mean(radar_i);
        radar_i_centered = radar_i_centered(:);
        radar_q_centered = radar_q - (max((radar_q))+min((radar_q)))/2;%mean(radar_q);
        radar_q_centered = radar_q_centered(:);
        data = [radar_i_centered radar_q_centered];
       
        radVec= vecnorm(data,2,2);
        meanR= mean(radVec);
        data2= data./meanR;
        C = cov(data);
        [V, D] = eig(C);
        [eigvals, idx] = sort(diag(D), 'ascend');
        D = diag(eigvals);
        V = V(:, idx);
        % Ensure first eigenvector points "right" (positive x component)
        if V(1,1) < 0
            V(:,1) = -V(:,1);
        end
        % Ensure right-handed basis (no mirror flip)
        if det(V) < 0
            V(:,2) = -V(:,2);
        end
        %normalizing 

        D_inv_sq = diag(1 ./ sqrt(diag(D)));
        % by using A = VDV' 
        transformed_C = (D_inv_sq * V' * data')' ;
        radar_i_fitted =  transformed_C(:,1);
        radar_q_fitted =  transformed_C(:,2);

        data_fit = [radar_i_fitted radar_q_fitted];

        % Option A: scale so the AVERAGE radius is 1
        r = sqrt(sum(data_fit.^2, 2));   % N×1
        s = 1 / mean(r);
        data_fit = s * data_fit;
        
        radar_i_fitted = data2(:,1); %data_fit(:,1);
        radar_q_fitted = data2(:,2); %data_fit(:,2);
  
        % we skip the whole covariance and mormalize by the mean of the
        % radiuses of the i-q vectors.



        if plotBool
            EllipsePlotter(radar_i, radar_q,radar_i_centered,radar_q_centered,radar_i_fitted,radar_q_fitted)
        end
        
end
