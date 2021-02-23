clc
clear

ndx       = 31;
dlambda   = 3;
dtheta    = 3;

projection= 'Equiangular'; % Choose from 'Equidistant' or 'Equiangular'

R         = 6371229.0;
d2r       = pi/180;

dlambda   = dlambda*d2r;
dtheta    = dtheta *d2r;

a         = R/sqrt(3); % Length of the cubic edges

if strcmp(projection,'Equidistant')
    % Equidistant
    dx        = 2 * a / ndx;
    x         = -a:dx:a;
    y         = x;
elseif strcmp(projection,'Equiangular')
    % Equiangular
    lambda      = -pi/4-dlambda:dlambda:pi/4+dlambda;
    theta       = -pi/4-dtheta :dtheta :pi/4+dtheta;
    x           = a*tan(lambda);
    y           = a*tan(theta );
end

nx        = size(x,2);
ny        = size(y,2);
a_matrix  = ones(nx,ny)*a;

x         = repmat(x,ny,1);
y         = repmat(y,nx,1)';

r         = sqrt(a_matrix.^2 + x.^2 + y.^2);

cart_coord(1,:,:) = R./r.*a_matrix;
cart_coord(2,:,:) = R./r.*x;
cart_coord(3,:,:) = R./r.*y;

% face 1
X1         = squeeze(cart_coord(1,:,:));
Y1         = squeeze(cart_coord(2,:,:));
Z1         = squeeze(cart_coord(3,:,:));

xx1      = a*Y1./X1;
yy1      = a*Z1./X1;

ghost_X1l = X1(:,1);
ghost_Y1l = Y1(:,1);
ghost_Z1l = Z1(:,1);
ghost_X1r = X1(:,end);
ghost_Y1r = Y1(:,end);
ghost_Z1r = Z1(:,end);
ghost_X1t = X1(end,:);
ghost_Y1t = Y1(end,:);
ghost_Z1t = Z1(end,:);
ghost_X1b = X1(1,:);
ghost_Y1b = Y1(1,:);
ghost_Z1b = Z1(1,:);

ghost_x14 = -a*ghost_X1l./ghost_Y1l;
ghost_y14 = -a*ghost_Z1l./ghost_Y1l;
ghost_x12 = -a*ghost_X1r./ghost_Y1r;
ghost_y12 = a*ghost_Z1r./ghost_Y1r;
ghost_x15 = a*ghost_Y1t./ghost_Z1t;
ghost_y15 = -a*ghost_X1t./ghost_Z1t;
ghost_x16 = -a*ghost_Y1b./ghost_Z1b;
ghost_y16 = -a*ghost_X1b./ghost_Z1b;

inverse_lambda14 = atan(ghost_x14/a);
inverse_theta14  = atan(ghost_y14/a);

% plot(ghost_x14./xx4(:,end))
% plot(ghost_y14./yy4(:,end))

% plot3(X,Y,Z,'.','Color','b');
surf(X1,Y1,Z1,'EdgeColor','k','FaceColor','c')
hold on

% face 2
X2         = -squeeze(cart_coord(2,:,:));
Y2         = squeeze(cart_coord(1,:,:));
Z2         = squeeze(cart_coord(3,:,:));

xx2      = -a*X2./Y2;
yy2      = a*Z2./Y2;

% plot3(X,Y,Z,'.','Color','r');
surf(X2,Y2,Z2,'EdgeColor','k','FaceColor','w')
hold on

% face 3
X3         = -squeeze(cart_coord(1,:,:));
Y3         = -squeeze(cart_coord(2,:,:));
Z3         = squeeze(cart_coord(3,:,:));

xx3      = a*Y3./X3;
yy3      = -a*Z3./X3;

% plot3(X,Y,Z,'.','Color','k');
surf(X3,Y3,Z3,'EdgeColor','k','FaceColor','g')
hold on

% face 4
X4         = squeeze(cart_coord(2,:,:));
Y4         = -squeeze(cart_coord(1,:,:));
Z4         = squeeze(cart_coord(3,:,:));

xx4      = -a*X4./Y4;
yy4      = -a*Z4./Y4;

% plot3(X,Y,Z,'.','Color','c');
surf(X4,Y4,Z4,'EdgeColor','k','FaceColor','b')
hold on

% face 5
X5         = -squeeze(cart_coord(3,:,:));
Y5         = squeeze(cart_coord(2,:,:));
Z5         = squeeze(cart_coord(1,:,:));

xx5      = a*Y5./Z5;
yy5      = -a*X5./Z5;

% plot3(X,Y,Z,'.','Color','g');
surf(X5,Y5,Z5,'EdgeColor','k','FaceColor','r')
hold on

% face 6
X6         = squeeze(cart_coord(3,:,:));
Y6         = squeeze(cart_coord(2,:,:));
Z6         = -squeeze(cart_coord(1,:,:));

xx6      = -a*Y6./Z6;
yy6      = -a*X6./Z6;

% plot3(X,Y,Z,'.','Color','m');
surf(X6,Y6,Z6,'EdgeColor','k','FaceColor','m')
hold on


print('-dpng','-opengl','-r300',['MCV_grid','.png']);