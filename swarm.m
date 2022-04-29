close all;

f = figure();
N = 50;
x = 5*randn(1,N);
y = 5*randn(1,N);
v = 0.5*rand(1,N);
vnew = v;
theta = 2*pi*randn(1,N);
thetanew = theta;
rneighbor = 1; % search radius for neighbors
rsearch = 3;

p = plot(NaN,NaN,'.','MarkerSize',6);
xrange = 20;
yrange = 20;
xlim([-xrange/2 xrange/2]);
ylim([-yrange/2 yrange/2]);
set(p, 'XDataSource', 'x', 'YDataSource', 'y');
dt = 0.01;

% parameters for storing particle locations
search_grid_step = 0.2; % don't make too small for memory reasons
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

for i=1:10000
    % integrate position
    xnew = x + (dt*v).*cos(theta);
    ynew = y + (dt*v).*sin(theta);
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
                            if ri == 1
                                v_cart_neighbor(ri,:) = r^2/d*v_cart_neighbor(ri,:) + [v(np)*cos(theta(np)) v(np)*sin(theta(np))];
                            else
                                v_cart_neighbor(ri,:) = v_cart_neighbor(ri,:) + [v(np)*cos(theta(np)) v(np)*sin(theta(np))];
                            end
                            cm_neighbor(ri,:) = cm_neighbor(ri,:) + np_loc;
                            num_neighbor(ri) = num_neighbor(ri) + 1;
                        end
                    end
                end
            end
        end
        v_cart = [v(n)*cos(theta(n)) v(n)*sin(theta(n))];
        if num_neighbor(1) > 0
            % nearby neighbors found -> match velocity vector
            d_cm = cm_neighbor(1,:)./num_neighbor(1) - [x(n) y(n)];
            v_cart = v_cart_neighbor(1,:)./num_neighbor(1) + 0.01*d_cm;
            vnew(n) = 0.001*v(n) + 0.999*norm(v_cart) + 0.5*randn;
            d_theta = 0.1*randn;
        elseif num_neighbor(2) > 0
            % further neighbors found -> flock to them
            d_cm = cm_neighbor(2,:)./num_neighbor(2) - [x(n) y(n)];
            v_cart = 0.3*v_cart + 0.5*v_cart_neighbor(2,:) + 2*d_cm;
            vnew(n) = norm(v_cart) + 0.5*randn;
            d_theta = 0.5*randn;
        else
            % just wander around; increase randomness of d_theta and d_vn
            d_theta = 0.01*randn;
            vnew(n) = v(n) + 0.02*randn;
        end
        R = [cos(d_theta) -sin(d_theta); sin(d_theta) cos(d_theta)];
        v_cart = (R*v_cart')';
        if abs(vnew(n)) > 10
            vnew(n) = sign(vnew(n))*10;
        end
        thetanew(n) = angle(v_cart(1) + 1j*v_cart(2));
    end
    v = vnew;
    theta = thetanew;
    % update plot
    refreshdata(f,'caller');
    drawnow limitrate;
end