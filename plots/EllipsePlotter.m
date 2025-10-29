function [] = EllipsePlotter(radar_i, radar_q,radar_i_centered,radar_q_centered,radar_i_fitted,radar_q_fitted)
        figure('Name', 'Ellipse fitting me');
        subplot(3, 1, 1);
        scatter(radar_i, radar_q,'filled'); 
        title(['radar_i vs. radar_q - regular ']);
        xlabel('radar_i');
        ylabel('radar_q');
        grid on;
        subplot(3,1,2);
        scatter(radar_i_centered, radar_q_centered,'filled'); % 5 is the marker size
        title(['radar_i vs. radar_q - centered ']);
        xlabel('radar_i');
        ylabel('radar_q');
        grid on;
        subplot(3,1,3);
        scatter(radar_i_fitted, radar_q_fitted,'filled'); % 5 is the marker size
        title(['radar_i vs. radar_q - after transformation ']);
        xlabel('radar_i');
        ylabel('radar_q');
        grid on;
   end