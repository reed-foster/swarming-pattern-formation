close all;

f = figure();
N = 1000;
x = 5*randn(1,N);
y = 5*randn(1,N);
v = 5*randn(1,N);
theta = 2*pi*randn(1,N);
r = 1; % search radius for neighbors

p = plot(NaN,NaN,'.','MarkerSize',6);
xrange = 100;
yrange = 100;
xlim([-xrange/2 xrange/2]);
ylim([-yrange/2 yrange/2]);
set(p, 'XDataSource', 'x', 'YDataSource', 'y');
dt = 0.1;

% parameters for storing particle locations
search_grid_step = 0.5; % don't make too small for memory reasons
x_grid = int32(xrange/search_grid_step);
y_grid = int32(yrange/search_grid_step);
search_grid_contains_particle = cell(x_grid,y_grid); % list of particles in grid location
for xg=1:x_grid
    for yg=1:y_grid
        search_grid_contains_particle{xg,yg} = [];
    end
end
search_grid_particle_locs = zeros(x_grid,y_grid,N,2); % x and y coordinate of particle
outside_contains_particle = [];
outside_particle_locs = zeros(N,2);

for i=1:10000
    % integrate position
    xnew = x + (dt*v).*cos(theta);
    ynew = y + (dt*v).*sin(theta);
    % update particle locations in search grid
    for n=1:N
        if any(outside_contains_particle(:) == n)
            outside_contains_particle(outside_contains_particle == n) = [];
        else
            x_grid_loc = int32(floor((x(n)+xrange/2)/search_grid_step))+1;
            y_grid_loc = int32(floor((y(n)+yrange/2)/search_grid_step))+1;
            search_grid_contains_particle{x_grid_loc,y_grid_loc}(search_grid_contains_particle{x_grid_loc,y_grid_loc} == n) = [];
        end
        if abs(xnew(n)) >= xrange/2 || abs(ynew(n)) >= yrange/2
            outside_contains_particle(end+1) = n;
            outside_particle_locs(n,:) = [xnew(n) ynew(n)];
        else
            xnew_grid_loc = int32(floor((xnew(n)+xrange/2)/search_grid_step))+1;
            ynew_grid_loc = int32(floor((ynew(n)+yrange/2)/search_grid_step))+1;
            search_grid_contains_particle{xnew_grid_loc,ynew_grid_loc}(end+1) = n;
            search_grid_particle_locs(xnew_grid_loc,ynew_grid_loc,n,:) = [xnew(n) ynew(n)];
        end
    end
    x = xnew;
    y = ynew;
    % select velocity vector
    for n=1:N
        % determine neighbors from grid
        % https://diglib.eg.org/bitstream/handle/10.2312/cgvc20191258/055-063.pdf
        % first get cells to search
        xsearch_min = x(n)-r;
        xsearch_max = x(n)+r;
        ysearch_min = y(n)-r;
        ysearch_max = y(n)+r;
        num_neighbor = 0;
        v_neighbor = 0;
        theta_neighbor = 0;
        if xsearch_max >= xrange/2 || xsearch_min <= -xrange/2 || ysearch_max >= yrange/2 || ysearch_min <= -yrange/2
            % search outside
            for np=outside_contains_particle
                if np == n
                    continue
                end
                np_loc = outside_particle_locs(np,:);
                d = sum((np_loc - [x(n) y(n)]).^2);
                if d <= r^2
                    v_neighbor = v_neighbor + v(np);
                    theta_neighbor = theta_neighbor + theta(np);
                    num_neighbor = num_neighbor + 1;
                end
            end
        end
        xgrid_max = int32(floor((xsearch_max+xrange/2)/search_grid_step))+1;
        xgrid_min = int32(floor((xsearch_min+xrange/2)/search_grid_step))+1;
        ygrid_max = int32(floor((ysearch_max+yrange/2)/search_grid_step))+1;
        ygrid_min = int32(floor((ysearch_min+yrange/2)/search_grid_step))+1;
        for xg=unique([xgrid_min xgrid_max])
            for yg=unique([ygrid_min ygrid_max])
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
                    np_loc = squeeze(search_grid_particle_locs(xg,yg,np,:));
                    d = sum((np_loc - [x(n); y(n)]).^2);
                    if d <= r^2
                        v_neighbor = v_neighbor + v(np);
                        theta_neighbor = theta_neighbor + theta(np);
                        num_neighbor = num_neighbor + 1;
                    end
                end
            end
        end
        v_n = v(n);
        t_n = theta(n);
        if num_neighbor > 0
            v_n = v_neighbor/num_neighbor;
            t_n = theta_neighbor/num_neighbor;
        end
        v(n) = 2; %(25/(x(n)^2 + y(n)^2))*0.2*(rand - 0.5) + v_n;
        theta(n) = 0.9*mod(theta(n)*0.8 + (2*pi*0.02)*randn + 0.2*t_n,2*pi)-0.1*atan2(y(n), x(n));
    end
    % update plot
    refreshdata(f,'caller');
    drawnow limitrate;
end