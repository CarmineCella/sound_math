
% rawScatter3.m
%
function rawScatter3 (x_file, y_file, z_file, cluster_file, labels_file, centroids)
    x = rawread (x_file, 'float');
    y = rawread (y_file, 'float');
    z = rawread (z_file, 'float');
    
    labels = rawread (labels_file, 'int');
    %dispersions = rawread ('dispersions.raw', 'double');
    
    [cx, cy, cz] = textread (cluster_file, '%f %f %f', centroids);
   % [mx, my, mz] = textread (modelings_file, '%f %f %f', modelings);

    colors = labels;
    figure
    scatter3 (x, y, z, 30, colors, 'filled'); %, 'MarkerEdgeColor','k');    
    hold on
    sz = 50; %1 + (dispersions .* 100);
       
   % scatter3 (cx, cy, cz, sz, 'k', 'x', 'LineWidth', 3, 'MarkerEdgeColor','k')    
    
      for k=1:centroids
         scatter3 (cx(k), cy(k), cz (k), 'LineWidth',1, 'MarkerEdgeColor','k', 'MarkerFaceColor',[0 1 0], 'SizeData', 72)
        % [xe, ye, ze] = ellipsoid (cx(k), cy(k), cz(k), dispersions(k) * .001, dispersions(k) * .01, dispersions(k) * .001);
        % surfl (xe, ye, ze)
      end

    xlabel (x_file)
    ylabel (y_file)
    zlabel (z_file)
   
    figure
    subplot (3, 1, 1)
    plot (x)
    title (x_file)
    subplot (3, 1, 2)
    plot (y)
    title (y_file)
    subplot (3, 1, 3)
    plot (z)
    title (z_file)    
   
    %figure
    %patch (x, y, z)
    %title ('Patch representation of the features');
    
    %figure
    %subplot (2, 1, 1)
    %rose (labels)
    %title ('Labels');
    %subplot (2, 1, 2)
    %plot (dispersions)
    %title ('Dispersions')
    %hold on
    %dm = ones (length (dispersions), 1);
    %plot (dm .* mean (dispersions), 'r')    
end
    
% eof
