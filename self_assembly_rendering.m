clear
clc
set(0,'DefaultFigureWindowStyle','normal')
addpath('./Matlab');
addpath('./Matlab/Image_Analysis_v0626');
figure(10)
set(gcf,'position',[10 10 900 900])
movie = false;
data = load('traj/xyz0.dat');
pnum = max(data(:,1))+1;
a = 0.8;
frame = size(data,1)/pnum;
if movie
    v = VideoWriter('traj_demo.avi');
    v.FrameRate = 2;
    open(v);
end
cluster = [0,0; 0,-1];
for f = 1:10:frame
    curr_frame = data((f-1)*pnum+1:f*pnum,:);
    
    % plot background
    clf
    hold on
    plot_electrodes(0);
    
    % separately handle each cluster
    for i = 1:size(cluster,1)
        % pick indices of particles within the same cluster
        idx = round(curr_frame(:,2)*a/60) == cluster(i,1) & ...
              round(curr_frame(:,3)*a/60) == cluster(i,2);
        cluster_frame = curr_frame(idx,:);
        
        % image analysis on single cluster
        [P_Type, G_PSI6, G_C6, L_PSI6, L_C6] = IMAGE_ANALYSIS(cluster_frame);
        cluster_frame = cluster_frame*a;
        
        % morphology parameters of single cluster
        norm_cluster_frame = cluster_frame(:,2:3) - mean(cluster_frame(:,2:3));
        phi = fminsearch(@(t)minGry(norm_cluster_frame,t),22.5);
        R = [cos(phi*pi/180), -sin(phi*pi/180);...
             sin(phi*pi/180), cos(phi*pi/180)];
        lambda = diag((norm_cluster_frame*R)'*(norm_cluster_frame*R)/pnum);
        if lambda(1) < lambda(2)
            phi = 90+ phi;
            lambda = flipud(lambda);
        end
        c = abs(diff(sqrt(lambda)))/sqrt(sum(lambda));
        
        % plot all particles
        plot(cluster_frame(:,2),cluster_frame(:,3),'ko',...
            'markersize', 6,...
            'MarkerFaceColor','k');
        
        % plot grain boudary
        if G_C6>=0.8 && G_PSI6 <= 0.93
            gb_ptr = P_Type == -3;
            % grain boundary positions
            x_gb = cluster_frame(gb_ptr,2);
            y_gb = cluster_frame(gb_ptr,3);
            % grain boundary center
            x_gb_mean = mean(x_gb);
            y_gb_mean = mean(y_gb);
            % linear regression slope and intercept
            slope_gb = sum((x_gb-x_gb_mean).*(y_gb-y_gb_mean))...
                /sum((x_gb-x_gb_mean).^2);
            intercept_gb = y_gb_mean - slope_gb*x_gb_mean;
            % grain boundary line
            alpha = atan(slope_gb);
            norm = [cos(alpha), sin(alpha)];
            % positive and negative range of grain boundary
            max_gb = max((x_gb-x_gb_mean).*norm(1) + (y_gb-y_gb_mean).*norm(2));
            min_gb = min((x_gb-x_gb_mean).*norm(1) + (y_gb-y_gb_mean).*norm(2));
            % plot grain boundary particles
            plot(cluster_frame(gb_ptr,2), cluster_frame(gb_ptr,3),...
                'o',...
                'markersize', 6,...
                'MarkerFaceColor','y',...
                'MarkerEdgeColor','y');
            % plot grain boundary line
            plot((min_gb:max_gb)*cos(alpha)+x_gb_mean,...
                (min_gb:max_gb)*sin(alpha)+y_gb_mean,'-',...
                'color',[1, .8, 0],...
                'linewidth',4)
        end
        
        % plot morphology
        % cluster center position
        x_morph_mean = mean(cluster_frame(:,2));
        y_morph_mean = mean(cluster_frame(:,3));
        % scaled particle positions
        x_morph = cluster_frame(:,2) - x_morph_mean;
        y_morph = cluster_frame(:,3) - y_morph_mean;
        % morphology long-axis direction
        norm = [cos(phi*pi/180), sin(phi*pi/180)];
        % positive and negative range along long-axis
        max_morph = max(x_morph.*norm(1) + y_morph.*norm(2));
        min_morph = min(x_morph.*norm(1) + y_morph.*norm(2));
        % end points of position range
        morph1 = [max_morph*norm(1)+x_morph_mean,...
            max_morph*norm(2)+y_morph_mean];
        morph2 = [min_morph*norm(1)+x_morph_mean,...
            min_morph*norm(2)+y_morph_mean];
        dp = morph2 - morph1;
        % quiver connect end points in both directions
        quiver(morph1(1),morph1(2),dp(1),dp(2),0,...
            'c','linewidth',3,'MaxHeadSize',.4)
        quiver(morph2(1),morph2(2),-dp(1),-dp(2),0,...
            'c','linewidth',3,'MaxHeadSize',.4)
    end
    
    % figure customization
    xticks('')
    yticks('')
    axis([-90 91 -90 90])
    pbaspect([1,1,1])
    sb = subplot(1, 1, 1);
    set(sb,'position',[0 0 1 1])
    box on
    hold off
    annotation('textbox', ...
        [.9 .0 .1 .1], ...
        'String', num2str(f/10), ...
        'EdgeColor', 'none', ...
        'BackgroundColor','w',...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom',...
        'FontSize', 16);
    drawnow;
    if movie
        this_frame = getframe(gcf);
        writeVideo(v,this_frame.cdata);
    end
end
if movie
    close(v)
end


%         gyr = cluster_frame'*cluster_frame/pnum;
%         S = eig(gyr);
%         c = abs(sqrt(S(1)) - sqrt(S(2)))/sqrt(S(1)) + sqrt(S(2));