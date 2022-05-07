close all;

Nframes = 1500;
f = figure();
N = 150;
x = 10*randn(1,N);
y = 10*randn(1,N);
V = 5; % velocity magnitude (all particles have the same velocity magnitude)
v = 5*randn(N,2);
for n=1:N
    v(n,:) = v(n,:)./norm(v(n,:));
end
vnew = v;
rneighbor = 1; % search radius for neighbors
rsearch = 2;
% weighting of velocity/cm favoritism based on neighbor location relative to
% velocity vector direction
[vweight_a, vweight_b] = make_logistic(0.1, 0.99);
[cmweight_a, cmweight_b] = make_logistic(0.7, 0.9);
[vdotvpweight_a, vdotvpweight_b] = make_logistic(0.05, 0.999);
[leaderweight_a, leaderweight_b] = make_logistic(0.5, 0.9);
leader_factor = zeros(1,N);

colors = zeros(4,3*N,'uint8');
colors(4,:) = 255;
cmap = hot(256);

q = quiver(x,y,0.8/V*v(:,1)',0.8/V*v(:,2)',0,'LineWidth',2.5);
xrange = 20;
yrange = 20;
xlim([-xrange/2 xrange/2]);
ylim([-yrange/2 yrange/2]);
dt = 0.01;
txt = {strcat("N = ", num2str(N), ", |v| = ", num2str(V), ", \Deltat = ", num2str(dt)), strcat("R_{neighbor} = ", num2str(rneighbor), ", R_{search} = ", num2str(rsearch))};
text(-xrange/2+0.5,-yrange/2+0.5,txt,'HorizontalAlignment','left','VerticalAlignment','bottom');

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

display(x);
display(y);

% parameters for storing particle locations
search_grid_step = 0.5;%floor(sqrt(rneighbor*rsearch)); % don't make too small for memory reasons
x_grid = int32(xrange/search_grid_step);
y_grid = int32(yrange/search_grid_step);
search_grid_contains_particle = cell(x_grid,y_grid); % list of particles in grid location
for xg=1:x_grid
    for yg=1:y_grid
        search_grid_contains_particle{xg,yg} = [];
    end
end
search_grid_particle_locs = zeros(x_grid,y_grid,N,2); % x and y coordinate of particle

for t=1:Nframes
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
        cm_neighbor_uw = zeros(2,2); % unweighted center of mass
        avg_v_dot_vp = [0 0];
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
                            vdot_vp_weight = 1/(1+exp(vdotvpweight_a*norm_dot+vdotvpweight_b));
                            avg_v_dot_vp(ri) = avg_v_dot_vp(ri) + dot(v(n,:)/V,v(np,:))*vdot_vp_weight; % weight velocity alignment of neighbors behind current particle more strongly
                            v_cart_neighbor(ri,:) = v_weight*v_cart_neighbor(ri,:) + v(np,:);
                            cm_neighbor(ri,:) = cm_weight*cm_neighbor(ri,:) + np_loc;
                            cm_neighbor_uw(ri,:) = cm_neighbor_uw(ri,:) + np_loc;
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
            leader_factor(n) = 1/(1+exp(-leaderweight_a*avg_v_dot_vp(1)/num_neighbor(1)+leaderweight_b));
        elseif num_neighbor(2) > 0
            % further neighbors found -> flock to them
            d_cm = cm_neighbor(2,:)./num_neighbor(2) - [x(n) y(n)];
            v_neighbor = v_cart_neighbor(2,:)./num_neighbor(2);
            v_urge = 0.2;
            rand_urge = 0.05;
            cm_urge = 0.5;
            leader_factor(n) = 1/(1+exp(-leaderweight_a*avg_v_dot_vp(2)/num_neighbor(2)+leaderweight_b));
        else
            % just wander around; increase randomness of d_theta and d_vn
            d_cm = [0 0];
            v_neighbor = [0 0];
            rand_urge = 0.1;
            v_urge = 0;
            cm_urge = 0;
            leader_factor(n) = 0;
        end
        leader_factor(n) = 2*(10^(10*(leader_factor(n)-1))+0.5*leader_factor(n))/3;
        rand_urge = rand_urge + 0.1*leader_factor(n)^5;
        v_urge = (1 - 0.8*leader_factor(n))*v_urge;
        cm_urge = (1 - 0.5*leader_factor(n))*cm_urge;
        vnew(n,:) = v(n,:) + v_urge*v_neighbor + rand_urge*randn(1,2) + cm_urge*d_cm;
        vnew(n,:) = V.*vnew(n,:)./norm(vnew(n,:));
    end
    %leader_idx = uint8(256*leader_factor)+1;
    v = vnew;
    % update plot
    set(q, 'xdata', x, 'ydata', y, 'udata', 0.8/V*v(:,1)', 'vdata', 0.8/V*v(:,2)');
    % hacky color setting, works in R2021b, may break in the future though
    % heads take a 4x3N vector (each 3 consective 4x3 chunks is one color used
    % for coloring the 3 points of the arrow head), tails take a 4x2N vector
    leader_factor = 0.5*leader_factor;
    colors(1:3,1:3:3*N) = uint8(256*cmap(uint8(256*leader_factor)+1,:)');
    colors(1:3,2:3:3*N) = uint8(256*cmap(uint8(256*leader_factor)+1,:)');
    colors(1:3,3:3:3*N) = uint8(256*cmap(uint8(256*leader_factor)+1,:)');
    set(q.Head,'ColorBinding','interpolated','ColorData',colors);
    colors(1:3,1:2:2*N) = uint8(256*cmap(uint8(256*leader_factor)+1,:)');
    colors(1:3,2:2:2*N) = uint8(256*cmap(uint8(256*leader_factor)+1,:)');
    set(q.Tail,'ColorBinding','interpolated','ColorData',colors(:,1:2*N));
    refreshdata(f,'caller');
    drawnow limitrate;
    F(t) = getframe;
end

fname = strcat(datestr(now,'yymmdd_HHMMSS_'),'N',num2str(N),'_v',num2str(V),'_dt',num2str(dt),'_Rn',num2str(rneighbor),'_Rs',num2str(rsearch),'.avi');
writerObj = VideoWriter(fname);
open(writerObj);
writeVideo(writerObj, F);
close(writerObj);