function plot_MJHP_2theta(filename, MINY, MAXY)

fileID = fopen(filename,'r');
filestring = split(filename,["_", "_", "."]);
compound = filestring(1);
type = filestring(2);
theta = filestring(3);
mjout = filestring(4);
type = strcat("(", type);
type = strcat(type, " type)");
fig_title = strcat(compound, "  ");
fig_title = strcat(fig_title, type);

%Header
str = fscanf(fileID, '%s',1); %comp name
str = fscanf(fileID, '%s',1); %astar
astar = fscanf(fileID,'%f %f %f',[1 3]);
str = fscanf(fileID, '%s',1); %bstar
bstar = fscanf(fileID,'%f %f %f',[1 3]);
str = fscanf(fileID, '%s',1); %cstar
cstar = fscanf(fileID,'%f %f %f',[1 3]);
str = fscanf(fileID, '%s',1); %number of steps
nstep = fscanf(fileID, '%d',1);
str = fscanf(fileID, '%s',1); %starting Energy
Elow = fscanf(fileID, '%f',1);
str = fscanf(fileID, '%s',1); %size of step
Emesh = fscanf(fileID, '%f',1);
str = fscanf(fileID, '%s',1); %Fermi degree
deg_Ef = fscanf(fileID, '%f',1);
str = fscanf(fileID, '%s', 1); %number of reflections
npxrd = fscanf(fileID, '%d', 1); 
%End of header

%Initialize Arrays
twotheta = [];
H = [];
K = [];
L = [];
Estep = [];
rflcE = [];

%Read in pxrdE information
for np = 1:npxrd
    str = fscanf(fileID, '%s', 1); %pxrd#
    twotheta_val = fscanf(fileID, '%f', 1);
    H_val = fscanf(fileID, '%d', 1);
    K_val = fscanf(fileID, '%d', 1);
    L_val = fscanf(fileID, '%d',1);
    twotheta = [twotheta, twotheta_val];
    H = H_val;
    K = K_val;
    L = L_val;

    %pxrd by Estep
    size_rflc = [2 Inf];
    RFLC = fscanf(fileID, '%d %f', size_rflc);
    step = RFLC(1,:);
    rflcE(np,:) = RFLC(2,:);
    Estep = (step + Elow)/Emesh;

end

%initializing data for gaussian
Z = 0;
Ztot=0;
Etot = 0;
miny = min(Estep);
if ((miny < MINY) && (MINY ~= 0))
    miny = MINY
end
maxy = max(Estep);
if ((maxy > MAXY) && (MAXY ~=0))
    maxy = MAXY
end
% miny = -10;

minx = 0;
maxx = 60;
[X,Y] = meshgrid(minx:.05:maxx,miny:.05:maxy);

for n = 1:npxrd
    stmt = ['Calculating for powder reflection#  ', num2str(n)];
    disp(stmt);
    for dE = 1:nstep
        RFpeak = rflcE(n,dE);
        theta = twotheta(n);
        energy = Estep(dE);
        sigma_x = 0.8;
        sigma_y = 0.3;
        Egauss = RFpeak.*exp(-(((X-theta)./sigma_x).^2+((Y-energy)./sigma_y).^2));
        Etot = Etot + Egauss;
    end
    Z = Z+Etot;
end

maxZ = max(Z,[],'all');
minZ = min(Z,[],'all');
abs_minZ = abs(minZ);
if maxZ > abs_minZ
    norm = maxZ;
else
    norm = abs_minZ;
end
Z = Z./norm;

%Plotting MJpxrd surface
figure;
surf(X,Y,Z,'edgecolor','none');
%setting co
red = [1,0,0];
white = [1,1,1];
blue = [0,0,1];
length = 100;
arr_blue = [linspace(blue(1),white(1),length)', linspace(blue(2),white(2),length)', linspace(blue(3),white(3),length)'];
arr_red =  [linspace(white(1),red(1),length)', linspace(white(2),red(2),length)', linspace(white(3),red(3),length)'];   
colors = cat(1, arr_blue, arr_red);
colormap(colors);
% clim([-1 1]);
clim([-.5 .5]);
colorbar;
colorbar('Ticks', ({})); %no ticks on colorbar;
view([0,0,1]);
hold on;

% %plotting Ef
fermi = 0.0;
plot3([0,60],[fermi,fermi],[maxZ, maxZ],':','linewidth',2,'color',[0,0,0]);
plot3([deg_Ef, deg_Ef],[maxy, miny], [maxZ, maxZ],':','linewidth',2.5,'color', [0.6,0.6,0.6]);
hold on;
plot3([0,60],[fermi,fermi],[minZ, minZ],':','linewidth',2,'color',[0,0,0]);
plot3([deg_Ef, deg_Ef],[maxy, miny], [minZ, minZ],':','linewidth',2.5,'color',[0.6, 0.6, 0.6]);
hold on

miny
maxy
%making the plot pretty
axis([0 60 miny maxy]);
set(gca,'Fontsize',13);
title('\fontsize{13}', fig_title);
xlabel('\fontsize{18} 2\theta','Interpreter','tex')
ylabel('\fontsize{18} Energy (eV)');
fontname(gcf,"Arial")
yticks({});
plot([0,0],[maxy,miny],'linewidth',1,'color',[0,0,0]);
plot([0,60],[maxy,maxy],'linewidth',1,'color',[0,0,0]);
plot([0,60],[miny, miny],'linewidth',1,'color',[0,0,0]);
plot([60,60],[maxy,miny],'linewidth',1, 'color',[0,0,0]);
fclose('all');
