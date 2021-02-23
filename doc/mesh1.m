clc
clear

resolution = 5;

R         = 6371229.0;
d2r       = pi/180;

a         = R/sqrt(3);

dlambda   = resolution*d2r;
dtheta    = resolution*d2r;
lambda    = -pi:dlambda:pi;
theta     = -pi/2:dtheta:pi/2;

nx        = size(lambda,2);
ny        = size(theta ,2);

lambda    = repmat(lambda,ny,1);
theta     = repmat(theta ,nx,1)';
[X,Y,Z]   = sph2cart(lambda,theta,R);

x         = a*tan(lambda);
y         = a*tan(theta).*sec(lambda);
a_matrix  = ones(size(x))*a;

point_idx = abs(x)<=a & abs(y)<=a;

% Plot Origin point
plot3(0,0,0,'o','Color','k')
hold on

% Plot Cube
plot3(a_matrix(point_idx),x       (point_idx),y(point_idx)       ,'.','Color','b')
hold on
plot3(x       (point_idx),a_matrix(point_idx),y(point_idx)       ,'.','Color','b')
hold on
plot3(x       (point_idx),y       (point_idx),a_matrix(point_idx),'.','Color','b')
hold on
plot3(-a_matrix(point_idx),x       (point_idx),y(point_idx)       ,'.','Color','b')
hold on
plot3(x       (point_idx),-a_matrix(point_idx),y(point_idx)       ,'.','Color','b')
hold on
plot3(x       (point_idx),y       (point_idx),-a_matrix(point_idx),'.','Color','b')
hold on


% Plot Sphere
plot3(X,Y,Z,'+','Color','r')