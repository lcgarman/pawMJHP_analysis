function plot_MJHP_HKL(filename, MINY, MAXY,path_direc)

 %%Open the *.mjout file
fileID = fopen(filename,'r');
%split file name to become title 
filestring = split(filename,["_", "_", "."]);
compound = filestring(1);
type = filestring(2);
type = strcat("(", type);
type = strcat(type, " type)");
reflection = filestring(3);
reflection = strcat("[", reflection);
reflection = strcat(reflection, "]");
mjout = filestring(4);
fig_title = strcat(compound, "  ");
fig_title = strcat(fig_title, type);
fig_title = strcat(fig_title, "  ");
fig_title = strcat(fig_title, reflection);

%Header
str = fscanf(fileID, '%s',1); %file name
str = fscanf(fileID, '%s',1); %astar
astar = fscanf(fileID,'%f %f %f',[1 3]);
str = fscanf(fileID, '%s',1); %bstar
bstar = fscanf(fileID,'%f %f %f',[1 3]);
str = fscanf(fileID, '%s',1); %cstar
cstar = fscanf(fileID,'%f %f %f',[1 3]);
str = fscanf(fileID, '%s',1); %number of steps
totstep = fscanf(fileID, '%d',1);
str = fscanf(fileID, '%s',1); %starting Energy
Elow = fscanf(fileID, '%f',1);
str = fscanf(fileID, '%s',1); %size of step
Emesh = fscanf(fileID, '%f',1);
str = fscanf(fileID, '%s',1); %HKL
HKL = fscanf(fileID,'%d %d %d',[1 3]);
%End of header

HATOEV = 27.211396;
%Read in reflection Energies
size_rflc = [2 Inf];
rflc = fscanf(fileID,'%d %f',size_rflc);
step = rflc(1,:);
%convert mjhp values from Ha to eV
mjhp_val = rflc(2,:)*HATOEV;
%convert step number into an energy value
Yval = (step + Elow)/Emesh;

%find maximum mjhp values 
mjhp_val = mjhp_val(Yval<MAXY & Yval>MINY);
Yval = Yval(Yval<MAXY & Yval>MINY);

maxMJ = max(mjhp_val);
minMJ = min(mjhp_val);
abs_minMJ = abs(minMJ);
if maxMJ > abs_minMJ
    heightMJ = maxMJ;
else
    heightMJ = abs_minMJ;
end

%initialize varibale for gaussian
nsts = size(Yval);
nsts = nsts(2);
miny=  MINY;
maxy = MAXY;
Eval = miny:0.001:maxy;
nEval = size(Eval);
nEval = nEval(2);

%perform gaussian over mjhp contributions
for i = 1:nEval
    gauss_sum = 0;
    for j = 1:nsts
        sigma = 0.1;%standard deviation; controls width 
        Ei = Eval(i);
        mjhp_j = mjhp_val(j);
        Yval_j = Yval(j);
        gauss_mjhp = mjhp_j*exp(-(Ei-Yval_j)^2/(2*sigma^2));
        gauss_sum = gauss_mjhp + gauss_sum;
    end
    smooth_mjhp(i) = gauss_sum;
end

%find max and min of the smoothed data
max_smooth = max(smooth_mjhp);
min_smooth = min(smooth_mjhp);
abs_min_smooth = abs(min_smooth);

%normalize the smoothed data so that max and min match original max/minMJ
for i = 1:nEval
    if smooth_mjhp(i)>0
        conserve_mjhp(i) = smooth_mjhp(i) * maxMJ / max_smooth;
        mjhp_color(i) = conserve_mjhp(i)/maxMJ;
    else
        conserve_mjhp(i) = smooth_mjhp(i) * abs_minMJ / abs_min_smooth;
        mjhp_color(i) = conserve_mjhp(i)/abs_minMJ;
    end
end

figure('Position', [360,198,300,400]);
%plot reflection energy values
red = [1,0,0];
blue = [0,0,1];
black = [0,0,0];
length = 100;
arr_blue = [linspace(blue(1),black(1),length)', linspace(blue(2),black(2),length)', linspace(blue(3),black(3),length)'];
arr_red = [linspace(black(1),red(1),length)', linspace(black(2),red(2),length)', linspace(black(3),red(3),length)'];
blue_to_red = cat(1, arr_blue, arr_red);
col = mjhp_color;
%plot the smoothed + normalized mjhp
fig = patch([-conserve_mjhp nan],[Eval nan],[col nan],[col nan],'edgecolor', 'interp', 'linewidth',2.8); 
% fig = fill3(-conserve_mjhp, Eval, col , col, 'edgecolor', 'interp', 'linewidth', 2, 'FaceColor','interp','FaceVertexAlphaData', alphadata,'FaceAlpha','interp'); 
colormap(blue_to_red);
% colorbar;
hold on;

% % plot data without gaussian
% % plot(-mjhp_val, Yval);

minx = -heightMJ;
maxx = heightMJ;
fermi = 0.0;

% % % making the plot pretty
hold on;
view([0,90]);
axis([minx maxx miny maxy]);
x_range = maxx-minx;
x_step = x_range/4;
set(gca,'TickDir','out');
% xticks(minx:x_step:maxx);
% xticks([minx, 0, maxx]);

xtickformat('%.2f')
title(['\fontsize{12}' fig_title]);
xlabel('\fontsize{18} -MJHP (eV)');
ylabel('\fontsize{18} Energy (eV)');
yticks({});
fontname(gcf,"Arial");
set(gca,'Fontsize',13);
plot([0,0],[miny, maxy],':', 'linewidth',1,'color',[0,0,0]);
plot([minx,maxx],[fermi,fermi],':', 'linewidth',2,'color',[0,0,0]);
plot([minx, minx],[miny,maxy],'linewidth',1,'color',[0,0,0]);
plot([minx,maxx],[miny,miny],'linewidth',1,'color',[0,0,0]);
plot([minx,maxx],[maxy, maxy],'linewidth',1,'color',[0,0,0]);
plot([maxx,maxx],[miny,maxy],'linewidth',1,'color',[0,0,0]);

%%plot additional lines on the plot
% plot([minx,maxx],[-2.1, -2.1],':','linewidth',2,'color',blue);
% plot([minx,maxx],[2.7, 2.7],':','linewidth',2,'color',red);

%save figure in directory
if path_direc ~= 0
    fileinput = split(filename,["."]);
    filestore = fileinput(1);
    filestore = strcat(filestore, "-MJHP");
    filestore_tif = strcat(filestore, ".tif");
    filestore_fig = strcat(filestore, ".fig");
    saveas(fig,fullfile(path_direc,[filestore_tif]));
    saveas(fig,fullfile(path_direc,[filestore_fig]));
end

fclose('all');


