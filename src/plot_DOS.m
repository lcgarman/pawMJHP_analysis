function plot_DOS(filename, option, MINY, MAXY)
% %option=1->no guass directly plot
% %option2->perform gauss

fileID = fopen(filename,'r');
for k=1:7
    tline = fgets(fileID);
end
str = fscanf(fileID, '%s %s %s %s',4); %Fermi energy preamble
fermi = fscanf(fileID, '%f',1)
for k=9:16
    tline = fgets(fileID);
end

% filestring = split(filename,["-", "."]);
% compound = filestring(1);
% psp = filestring(2);
% dosfile = filestring(3);
% psp = strcat("(", psp);
% psp = strcat(psp, ")");
% fig_title = strcat(psp, "  ");
% fig_title = strcat(fig_title, compound);

%Read in density info
sizedos = [3 Inf];
dos = fscanf(fileID,'%f %f %f',sizedos);
energy_Ha = dos(1,:);
density = dos(2,:);
integrated = dos(3,:);
energy = (energy_Ha-fermi) * 27.21139613;
fclose(fileID);

%plot data (1) NO gaussian
if (option == 1) 
    figure;
    plot3(energy, density, integrated, 'LineWidth', 2, 'Color',[0, 0, 0]);
    hold on;
    maxy = max(density);
end

%plot data (2) with gaussian
sigma = 0.001; %standard deviation; controls width
if (option == 2)
    nsts = size(density);
    nsts = nsts(2);
    maxx = max(energy);
    minx = min(energy);
    Eval = minx:0.001:maxx;
    nEval = size(Eval);
    nEval = nEval(2);
    for i = 1:nEval
        DOS = 0;
        for j = 1:nsts
            Ei = Eval(i);
            density_v = density(j);
            energy_v = energy(j);
            gauss_dos = density_v*exp(-(Ei-energy_v)^2/(2*sigma));
            DOS = gauss_dos + DOS;
        end
        DOSplot(i) = DOS;
    end

    figure('Position', [360,198,300,400]);
    plot(Eval, DOSplot,'linewidth',2.2,'color',[0 0 0]);
    DOS_belowEf = DOSplot(Eval<0);
    E_belowEf = Eval(Eval<0);
    hold on;
    area(E_belowEf, DOS_belowEf, 'FaceColor','black');
    maxy = max(DOSplot);
end

view([90 -90]);
maxy = 600;
miny = 0;
maxx = MAXY;
minx = MINY;

% % Print a line at the ideal electron count
% e_count = 1614;
% e_nearest = abs(e_count - integrated);
% min_nearest = min(e_nearest);
% pseudogap = energy(e_nearest == min_nearest);
% plot([pseudogap, pseudogap],[0,maxy],'--','linewidth',2,'color',[0.5,0.5,0.5]);

%making the plot pretty
set(gca,'Fontsize',13);
fontname(gcf,"Arial")
% title('\fontsize{18}', fig_title);
set(gca,'TickDir','out');
% xlabel('\fontsize{18}\fontname{arial} Energy (eV)');
set(gca,'ytick',[]);
hold on;
axis([minx maxx miny maxy]);
hold on;
fermi = 0.0;
plot([fermi,fermi], [miny,maxy],':','linewidth',2,'color',[0,0,0]);
plot([minx,minx],[miny,maxy],'linewidth',1,'color',[0,0,0]);
plot([minx,maxx],[miny,miny],'linewidth',1,'color',[0,0,0]);
plot([minx,maxx],[maxy, maxy],'linewidth',1,'color',[0,0,0]);
plot([maxx,maxx],[miny,maxy],'linewidth',1,'color',[0,0,0]);
fontname(gcf,"Arial")

fclose('all');