close all;
%Dipole field plot
p = 1;
e = 8.85*10^(-12);
x  =linspace(-0.5 , 0.5, 25);
z = linspace(-0.5 , 0.5, 25);
[X, Z ] = meshgrid(x,z);
R=sqrt(X.^2+Z.^2) ; 
EX =( p .* 3 .* X .* Z ./ R.^5 )./ (4.*pi.*e );
EZ = p./( 4 .* pi .* e ) .* ( 3.* Z.^2 ./R.^5 -1./ R.^3);

%// normalize the vectors so the arrows are visible
V = [EX(:) EZ(:)];
% Vn = bsxfun(@rdivide, V, sqrt(sum(V.^2,2)));
Vn=V;
Exn = reshape(Vn(:,1), size(EX));
Ezn = reshape(Vn(:,2), size(EZ));

%Create field that is rotated by 90 degrees
theta = 90; % rotate 90 degrees
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
ROTATED = R*[Exn(:)'; Ezn(:)'];
Exnr=reshape(ROTATED(1,:),size(EX));
Eznr=reshape(ROTATED(2,:),size(EZ));

%Add the rotated fields together:
Exn_tot = (Exnr + Exn)/2;
Ezn_tot = (Eznr + Ezn)/2;
theta = -45; % rotate 45 degrees
R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
ROTATED = R*[Exn_tot(:)'; Ezn_tot(:)'];
Exnr2=reshape(ROTATED(1,:),size(EX));
Eznr2=reshape(ROTATED(2,:),size(EZ));



figure(1);
quiver ( X , Z , Exn , Ezn , 'AutoScaleFactor',0.8,'linewidth',1) ; 
figure(2);
quiver ( X , Z , Exnr , Eznr , 'AutoScaleFactor',0.8,'linewidth',1) ; 
axis equal

%Dipoles can be expressed as sum of rotated dipoles??
figure(3);
quiver(X,Z,Exn_tot,Ezn_tot , 'AutoScaleFactor',0.8,'linewidth',1) ;
axis equal

figure(4);
quiver(X,Z,Exnr2,Eznr2 , 'AutoScaleFactor',0.8,'linewidth',1) ;
axis equal
%%
%Add streamlines:
p = 1;
e = 8.85*10^(-12);
x  =linspace(-0.5 , 0.5, 101);
z = linspace(-0.5 , 0.5, 101);
[X, Z ] = meshgrid(x,z );
R=sqrt(X.^2+Z.^2) ; 
EX =( p .* 3 .* X .* Z ./ R.^5 )./ (4.*pi.*e );
EZ = p./( 4 .* pi .* e ) .* ( 3.* Z.^2 ./R.^5 -1./ R.^3);

%// normalize the vectors so the arrows are visible
V = [EX(:) EZ(:)];
Vn = bsxfun(@rdivide, V, sqrt(sum(V.^2,2)));
Exn = reshape(Vn(:,1), size(EX));
Ezn = reshape(Vn(:,2), size(EZ));
[row, col] = find(isnan(EX));
EX(row,col) = 0;
EZ(row,col) = EZ(row-1,col);
hold on;

for i=[-6,-5,-4,-3,-2,-1,1,2,3,4,5,6]
    streamline1 = stream2(X,Z,EX,EZ,i*0.08,0,[0.01,100000]);
    streamline1 = cropToSingleRevolution(streamline1);
    plot(streamline1{1}(:,1),streamline1{1}(:,2),'r','linewidth',1);  
end

%%
m1 = [0;0;1];
m2 = [0;1;0];
r = [rand(1);rand(1);rand(1)];
r = r/norm(r);

dir1 = dot(m1,r)*r - m1;
dir2 = dot(m2,r)*r - m2;
sum = dir1+dir2;
theta = -45;
R = [1 0 0; 0 cosd(theta) -sind(theta);0 sind(theta) cosd(theta)];
rot = R*sum;
dir1
dir2
rot = rot/norm(rot)
