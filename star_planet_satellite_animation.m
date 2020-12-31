function[] = star_planet_satellite_animation()
% star_planet_satellite_animation : function to create an animation
% to model the orbitography and the gravitational fields of a
% star - planet - satellite system.
%
% Author & support nicolas.douillet (at) free.fr, 2007-2020.


% Computational parameters
sz = 51; % size of the space grid, odd number >= 3; default : 51
xy_centre = 1 + floor(0.5*sz); % position of the grid centre; default : 1 + floor(0.5*sz)

% Star position
xs = 0; % default : 0
ys = 0; % default : 0

% Planet initial position
xp = floor(0.225*(sz-1)); % default : floor(0.225*(sz-1))
yp = floor(0.225*(sz-1)); % default : floor(0.225*(sz-1));

resolution = 60; % nb angle steps over one loop; default : 60
Phi_step = pi/resolution; % rotation angle step
m = 0:Phi_step:2*pi-Phi_step;

sat_year_period_nb = 5; % satellite number of revolutions during one year of the planet; default : 5

% Planet elliptic path semi major and minor axes
a = 1; % default : 1
b = 1; % default : 1

planet_path = cat(1,a*cos(m)*xp-b*sin(m)*yp,...
                    a*sin(m)*xp+b*cos(m)*yp);
                
rp = sqrt(sum(planet_path.^2,1)); % planet_orbit_radius
sat_dst = 0.25*mean(rp); % satellite (relative) distance to planet; default : 0.25*mean(rp)

sat_orb = cat(1,cos((1-sat_year_period_nb)*m)*sat_dst - sin((1-sat_year_period_nb)*m)*sat_dst,...
                sin((1-sat_year_period_nb)*m)*sat_dst + cos((1-sat_year_period_nb)*m)*sat_dst);
          
sat_path = planet_path + sat_orb;


% Star & planet size and gravity function radius
star_weight = 24;   % default : 24
star_radius_function = 0.45*(sz-1);   % default : 0.45*(sz-1) 
planet_weight = 3*star_weight/4; % default : 3*star_weight/4
planet_radius_function = 0.45*(sz-1); % default : 0.45*(sz-1); relative alternative : planet_weight*star_radius_function/star_weight


% Display parameters
time_lapse = 6/resolution; % default : 6/resolution
title_text = 'Star - planet - satellite system gravitational fields and orbitography modelling';
title_on = true;
filename = 'star_planet_satellite_system.gif';
az = -22; % default : -22
el = 80;  % default :  80

% Display settings
h = figure;
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);
axis tight manual;
planet_path_on = true;        % default : true
sat_path_on = true;           % default : true
cmap = 'hot';                 % default : 'hot'
star_size = star_weight;      % default : star_weight
planet_size = 0.5*star_size;  % default : 0.5*star_size
sat_size = 0.5*planet_weight; % default : 0.5*planet_weight

% Static star space shape 
Z_star = compute_space_shape(sz,xs,ys,-star_weight,star_radius_function);


for s = 1:numel(m)
    
    Z_planet = compute_space_shape(sz,planet_path(1,s),planet_path(2,s),-planet_weight,planet_radius_function);
    Z = Z_star + Z_planet;
    
    mesh(Z), hold on;
    colormap(cmap);
    alpha 0;        
    
    % --- Star --- %
    plot3(xy_centre,xy_centre,0,'o','Color',[1 1 0],'Markersize',star_size,'Linewidth',star_size), hold on;
    
    % --- Planet --- %
    zp = griddata(1:sz,1:sz,Z,xy_centre+planet_path(1,:),xy_centre+planet_path(2,:));            
    
    if planet_path_on        
        plot3(xy_centre+cat(2,planet_path(1,:),planet_path(1,1)),xy_centre+cat(2,planet_path(2,:),planet_path(2,1)),cat(2,zp(1,:),zp(1,1)),'Color',[1 0 1],'Linewidth',2), hold on;
    end
    
    plot3(xy_centre+planet_path(1,s),xy_centre+planet_path(2,s),zp(1,s),'o','Color',[1 0 1],'Markersize',planet_size,'Linewidth',planet_size), hold on;
    
    % --- Satellite --- %
    zs = griddata(1:sz,1:sz,Z,xy_centre+sat_path(1,:),xy_centre+sat_path(2,:));
    
    if sat_path_on        
        plot3(xy_centre+cat(2,sat_path(1,:),sat_path(1,1)),xy_centre+cat(2,sat_path(2,:),sat_path(2,1)),cat(2,zs(1,:),zs(1,1)),'Color',[0 1 0],'Linewidth',2), hold on;
    end
    
    plot3(xy_centre+planet_path(1,s)+sat_orb(1,s),xy_centre+planet_path(2,s)+sat_orb(2,s),zs(1,s),'o','Color',[0 1 0],'Markersize',sat_size,'Linewidth',sat_size), hold on;        
    
    % --- Display settings --- %
    ax = gca;
    ax.Clipping = 'off';
    set(ax,'Color',[0 0 0]);
    
    axis off;
    
    if title_on
        title(title_text,'FontSize',20,'Color',[1 1 1]), hold on;
    end
    
    view(az,el);
    
    % --- .gif creation --- %
    drawnow;
    
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
        
    if s == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',Inf,'DelayTime',time_lapse);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',time_lapse);
    end
    
    clf;
    
end


end % star_planet_satellite_animation


function [z] = compute_space_shape(sz, xc, yc, height, rmax)


z = zeros(sz);

scale_min = -1;
scale_max = 1;
step = (scale_max-scale_min)/(sz-1);

xc = step*xc;
yc = step*yc;
rmax = step*rmax;

sample_vect = scale_min:step:scale_max;
[i,j] = meshgrid(sample_vect);
idx_vect = round((0.5*(sz-1))*sample_vect+0.5*(sz-1)+1);

r = sqrt((i-xc).^2+(j-yc).^2);
z(idx_vect,idx_vect) = height*((1-r).^4).*(4*r+1);
f = r >= rmax;
z(f) = 0;


end % compute_space_shape