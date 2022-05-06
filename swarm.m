close all;

f = figure();
N = 50;
x = 5*randn(1,N);
y = 5*randn(1,N);
V = 5; % velocity magnitude (all particles have the same velocity magnitude)
v = 5*randn(N,2);
for n=1:N
    v(n,:) = v(n,:)./norm(v(n,:));
end
vnew = v;
rneighbor = 0.5; % search radius for neighbors
rsearch = 3;
% weighting of velocity/cm favoritism based on neighbor location relative to
% velocity vector direction
vw_forward = 0.99;
vw_90degrees = 0.1;
cmw_forward = 0.9;
cmw_90degrees = 0.7;
vw_prod = vw_forward*vw_90degrees;
cmw_prod = cmw_forward*cmw_90degrees;
vweight_a = log((vw_forward-vw_prod)/(vw_90degrees-vw_prod));
vweight_b = log((1-vw_90degrees)/vw_90degrees);
cmweight_a = log((cmw_forward-cmw_prod)/(cmw_90degrees-cmw_prod));
cmweight_b = log((1-cmw_90degrees)/cmw_90degrees);

p = quiver(x,y,0.1*v(:,1)',0.1*v(:,2)',0,'LineWidth',2.5);
xrange = 20;
yrange = 20;
xlim([-xrange/2 xrange/2]);
ylim([-yrange/2 yrange/2]);
dt = 0.01;

% parameters for storing particle locations
search_grid_step = floor(sqrt(rneighbor*rsearch)); % don't make too small for memory reasons
x_grid = int32(xrange/search_grid_step);
y_grid = int32(yrange/search_grid_step);
search_grid_contains_particle = cell(x_grid,y_grid); % list of particles in grid location
for xg=1:x_grid
    for yg=1:y_grid
        search_grid_contains_particle{xg,yg} = [];
    end
end
search_grid_particle_locs = zeros(x_grid,y_grid,N,2); % x and y coordinate of particle

% wrap initial locations
for n=1:N
    if x(n) >= xrange/2
        x(n) = x(n) - xrange;
    elseif x(n) < -xrange/2
        x(n) = x(n) + xrange;
    end
    if y(n) >= yrange/2
        y(n) = y(n) - yrange;
    elseif y(n) < -yrange/2
        y(n) = y(n) + yrange;
    end
end

for t=1:800
    % integrate position
    xnew = x + (dt*v(:,1)');
    ynew = y + (dt*v(:,2)');
    for n=1:N
        if xnew(n) >= xrange/2
            xnew(n) = xnew(n) - xrange;
        elseif xnew(n) < -xrange/2
            xnew(n) = xnew(n) + xrange;
        end
        if ynew(n) >= yrange/2
            ynew(n) = ynew(n) - yrange;
        elseif ynew(n) < -yrange/2
            ynew(n) = ynew(n) + yrange;
        end
    end
    x_cm = sum(xnew)/N;
    y_cm = sum(ynew)/N;
    % update particle locations in search grid
    for n=1:N
        % clear old position
        x_grid_loc = int32(floor((x(n)+xrange/2)/search_grid_step))+1;
        y_grid_loc = int32(floor((y(n)+yrange/2)/search_grid_step))+1;
        search_grid_contains_particle{x_grid_loc,y_grid_loc}(search_grid_contains_particle{x_grid_loc,y_grid_loc} == n) = [];
        % update with new position
        xnew_grid_loc = int32(floor((xnew(n)+xrange/2)/search_grid_step))+1;
        ynew_grid_loc = int32(floor((ynew(n)+yrange/2)/search_grid_step))+1;
        search_grid_contains_particle{xnew_grid_loc,ynew_grid_loc}(end+1) = n;
        search_grid_particle_locs(xnew_grid_loc,ynew_grid_loc,n,:) = [xnew(n) ynew(n)];
    end
    x = xnew;
    y = ynew;
    % select velocity vector
    for n=1:N
        % determine neighbors from grid
        % https://diglib.eg.org/bitstream/handle/10.2312/cgvc20191258/055-063.pdf
        % first get cells to search
        num_neighbor = [0 0];
        v_cart_neighbor = zeros(2,2);
        cm_neighbor = zeros(2,2);
        for ri=1:2
            if ri == 1
                r = rneighbor;
            else
                r = rsearch;
            end
            xgrid_max = int32(floor((x(n)+r+xrange/2)/search_grid_step))+1;
            xgrid_min = int32(floor((x(n)-r+xrange/2)/search_grid_step))+1;
            ygrid_max = int32(floor((y(n)+r+yrange/2)/search_grid_step))+1;
            ygrid_min = int32(floor((y(n)-r+yrange/2)/search_grid_step))+1;
            for xg=unique(xgrid_min:xgrid_max)
                for yg=unique(ygrid_min:ygrid_max)
                    if xg < 1 || xg > x_grid
                        continue
                    end
                    if yg < 1 || yg > y_grid
                        continue
                    end
                    for np=search_grid_contains_particle{xg,yg}
                        if np == n
                            continue
                        end
                        np_loc = squeeze(search_grid_particle_locs(xg,yg,np,:))';
                        d = sum((np_loc - [x(n) y(n)]).^2);
                        if d <= r^2
                            % weight neighbors in front more strongly (i.e. follow the "leaders")
                            v_dot_dx = dot(np_loc - [x(n) y(n)], v(n,:)/V);
                            norm_dot = v_dot_dx/sqrt(d);
                            v_weight = 1/(1+exp(-vweight_a*norm_dot+vweight_b));
                            cm_weight = 1/(1+exp(-cmweight_a*norm_dot+cmweight_b));
                            v_cart_neighbor(ri,:) = v_weight*v_cart_neighbor(ri,:) + v(np,:);
                            cm_neighbor(ri,:) = cm_weight*cm_neighbor(ri,:) + np_loc;
                            num_neighbor(ri) = num_neighbor(ri) + 1;
                        end
                    end
                end
            end
        end
        if num_neighbor(1) > 0
            % nearby neighbors found -> match velocity vector
            d_cm = cm_neighbor(1,:)./num_neighbor(1) - [x(n) y(n)];
            v_neighbor = v_cart_neighbor(1,:)./num_neighbor(1);
            v_urge = 0.1;
            rand_urge = 0.01;
            cm_urge = 0.1;
        elseif num_neighbor(2) > 0
            % further neighbors found -> flock to them
            d_cm = cm_neighbor(2,:)./num_neighbor(2) - [x(n) y(n)];
            v_neighbor = v_cart_neighbor(2,:)./num_neighbor(2);
            v_urge = 0.2;
            rand_urge = 0.05;
            cm_urge = 0.5;
        else
            % just wander around; increase randomness of d_theta and d_vn
            d_cm = [0 0];
            v_neighbor = [0 0];
            rand_urge = 0.1;
            v_urge = 0;
            cm_urge = 0;
        end
        vnew(n,:) = v(n,:) + v_urge*v_neighbor + rand_urge*randn(1,2) + cm_urge*d_cm;
        vnew(n,:) = V.*vnew(n,:)./norm(vnew(n,:));
    end
    v = vnew;
    % update plot
    set(p, 'xdata', x, 'ydata', y, 'udata', 0.8/V*v(:,1)', 'vdata', 0.8/V*v(:,2)');
    refreshdata(f,'caller');
    drawnow limitrate;
    F(t) = getframe;
end

writerObj = VideoWriter('test2.avi');
open(writerObj);
writeVideo(writerObj, F);
close(writerObj);