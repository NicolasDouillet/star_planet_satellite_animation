function[] = star_planet_satellite_animation()
% star_planet_satellite_animation : function to create an animation
% to model the orbitography and the gravitational fields of a
% star - planet - satellite system.
%
% Author & support nicolas.douillet (at) free.fr, 2007-2020.


% Computational parameters
sz = 51; % grid size (odd number >= 3)
xy_centre = 0.5*(sz-1); % grid centre position

% Star position
xs = 0;
ys = 0;

% Planet initial position
xp = floor(0.225*(sz-1));
yp = floor(0.225*(sz-1));

sat_dst = 3*sqrt(2); % satellite distance to planet

resolution = 60; % nb angle steps over one loop
Phi_step = pi/resolution; % angular step
m = 0:Phi_step:2*pi-Phi_step;

sat_year_period_nb = 5; % satellite number of revolution during one year of the planet

% Ellipse semi major and minor axes values
a = 1;
b = 1;
planet_path = cat(1,a*cos(m)*xp-b*sin(m)*yp,...
                    a*sin(m)*xp+b*cos(m)*yp);

sat_orb = cat(1,cos((1-sat_year_period_nb)*m)*sat_dst - sin((1-sat_year_period_nb)*m)*sat_dst,...
                sin((1-sat_year_period_nb)*m)*sat_dst + cos((1-sat_year_period_nb)*m)*sat_dst);
          
sat_path = planet_path + sat_orb;

% Display parameters
time_lapse = 0.1;
title_text = 'Star - planet - satellite system gravitational fields and orbitography modelling';
title_on = true;
filename = 'star_planet_satellite_system.gif';
az = -22;
el = 80;

% Display settings
h = figure;
set(h,'Position',get(0,'ScreenSize'));
set(gcf,'Color',[0 0 0]);
axis tight manual;
planet_path_on = true;
sat_path_on = true;


for s = 1:length(m)
    
    Phi_star = CS_RBF(sz,xs,ys,-15,floor(0.45*(sz-1)));
    Phi_planet = CS_RBF(sz,planet_path(1,s),planet_path(2,s),-2,6.5);
    Phi = Phi_star + Phi_planet;
    
    mesh(Phi), hold on;
    colormap([0 1 1]);
    alpha 0;
    
    if planet_path_on
        plot3(xy_centre+cat(2,planet_path(1,:),planet_path(1,1)),xy_centre+cat(2,planet_path(2,:),planet_path(2,1)),zeros(1,length(m)+1),'Color',[1 0 1],'Linewidth',2), hold on;
    end
    
    % Star
    plot3(xy_centre,xy_centre,0,'o','Color',[1 1 0],'Linewidth',40), hold on;
    
    % Planet
    plot3(xy_centre+planet_path(1,s),xy_centre+planet_path(2,s),0,'o','Color',[1 0 1],'Linewidth',15), hold on;
    
    % Satellite
    zs = griddata(1:sz,1:sz,Phi,xy_centre+sat_path(1,:),xy_centre+sat_path(2,:));
    
    if sat_path_on        
        plot3(xy_centre+cat(2,sat_path(1,:),sat_path(1,1)),xy_centre+cat(2,sat_path(2,:),sat_path(2,1)),cat(2,zs(1,:),zs(1,1)),'Color',[0 1 0],'Linewidth',2), hold on;
    end
    
    plot3(xy_centre+planet_path(1,s)+sat_orb(1,s),xy_centre+planet_path(2,s)+sat_orb(2,s),zs(1,s),'o','Color',[0 1 0],'Linewidth',10), hold on;
    
    ax = gca;
    ax.Clipping = 'off';
    set(ax,'Color',[0 0 0]);
    
    axis off;
    
    if title_on
        title(title_text,'FontSize',16,'Color',[1 1 1]), hold on;
    end
    
    view(az,el);
    
    drawnow;
    
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    % Write to the .gif file
    if s == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',Inf,'DelayTime',time_lapse);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',time_lapse);
    end
    
    clf;
    
end


end % star_planet_satellite_animation


function phi = CS_RBF(sz, xc, yc, height, rmax)
%
% Author & support nicolas.douillet (at) free.fr, 2007-2020.


if ~mod(sz,2) % sz = odd number
    sz = sz+1;
end

phi = zeros(sz,sz);
scale_min = -1;
scale_max = 1;
step = (scale_max-scale_min)/(sz-1);
xc = step*xc;
yc = step*yc;
rmax = step*rmax;


for i = scale_min:step:scale_max  % loop on the matrix
    
    u = round((0.5*(sz-1))*i+0.5*(sz-1)+1);  % conversion into indices
    
    for j = scale_min:step:scale_max
        
        v = round((0.5*(sz-1))*j+0.5*(sz-1)+1);
        r = norm([i-yc;j-xc]);  % radius computation
        
        phi(u,v) = ((1-r)^4).*(4*r+1);
        
        if r >= rmax 
            
           phi(u,v) = 0;
           
        end        
        
    end
    
end

phi = height*phi;


end % CS_RBF