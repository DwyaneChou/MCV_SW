clc
clear

ndx      = 31;
dalpha   = 3;
dbeta    = 3;

projection= 'Equidistant'; % Choose from 'Equidistant' or 'Equiangular'

R         = 6371229.0;
d2r       = pi/180;

dalpha   = dalpha*d2r;
dbeta    = dbeta *d2r;

a         = R/sqrt(3);

if strcmp(projection,'Equidistant')
    % Equidistant
    dx        = 2 * a / ndx;
    x         = -a:dx:a;
    y         = x;
elseif strcmp(projection,'Equiangular')
    % Equiangular
    alpha      = -pi/4:dalpha:pi/4;
    beta       = -pi/4:dbeta :pi/4;
    x           = a*tan(alpha);
    y           = a*tan(beta );
end

nx        = size(x,2);
ny        = size(y,2);
a_matrix  = ones(nx)*a;

x         = repmat(x,nx,1);
y         = repmat(y,ny,1)';

r         = sqrt(a_matrix.^2+x.^2+y.^2);

cart_coord(1,:,:) = R./r.*a_matrix;
cart_coord(2,:,:) = R./r.*x;
cart_coord(3,:,:) = R./r.*y;

% face 1
X         = squeeze(cart_coord(1,:,:));
Y         = squeeze(cart_coord(2,:,:));
Z         = squeeze(cart_coord(3,:,:));

% plot3(X,Y,Z,'.','Color','b');
surf(X,Y,Z,'EdgeColor','k','FaceColor','c')
hold on

% face 2
X         = -squeeze(cart_coord(2,:,:));
Y         = squeeze(cart_coord(1,:,:));
Z         = squeeze(cart_coord(3,:,:));

% plot3(X,Y,Z,'.','Color','r');
surf(X,Y,Z,'EdgeColor','k','FaceColor','w')
hold on

% face 3
X         = -squeeze(cart_coord(1,:,:));
Y         = -squeeze(cart_coord(2,:,:));
Z         = squeeze(cart_coord(3,:,:));

% plot3(X,Y,Z,'.','Color','k');
surf(X,Y,Z,'EdgeColor','k','FaceColor','g')
hold on

% face 4
X         = squeeze(cart_coord(2,:,:));
Y         = -squeeze(cart_coord(1,:,:));
Z         = squeeze(cart_coord(3,:,:));

% plot3(X,Y,Z,'.','Color','c');
surf(X,Y,Z,'EdgeColor','k','FaceColor','b')
hold on

% face 5
X         = -squeeze(cart_coord(3,:,:));
Y         = squeeze(cart_coord(2,:,:));
Z         = squeeze(cart_coord(1,:,:));

% plot3(X,Y,Z,'.','Color','g');
surf(X,Y,Z,'EdgeColor','k','FaceColor','r')
hold on

% face 6
X         = squeeze(cart_coord(3,:,:));
Y         = squeeze(cart_coord(2,:,:));
Z         = -squeeze(cart_coord(1,:,:));

% plot3(X,Y,Z,'.','Color','m');
surf(X,Y,Z,'EdgeColor','k','FaceColor','m')
hold on


print('-dpng','-opengl','-r300',['MCV_grid','.png']);