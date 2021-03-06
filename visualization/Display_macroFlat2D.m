%%
%%--- Plot data computed with Vicsek (data must be in 'data')
%%

%% Parameters visualisation
choiceDensity = 1;		% 1: ρ(x),u(x)  2: ρ(x,y),u(x,y) 3: g(θ)
lengthArrow2D = .5;	        % for 2: ρ(x,y),u(x,y)
isTheoricCurveTheta = 0;	% for 3: g(θ)

shouldSave    = 0;
shouldPlotEnd = 1;


%%------------------ 0.1) Read parameters ------------------%%
%%----------------------------------------------------------%%
fid = fopen('../bin/PARAMETER_MicroVic_flat.txt');
C   = textscan(fid, '%s','delimiter', '\n');
N         = str2num(C{1}{5});
Lx        = str2num(C{1}{12});
Ly        = str2num(C{1}{13});
dt        = str2num(C{1}{18});
Time      = str2num(C{1}{19});
nSeed     = str2num(C{1}{27});
dx        = str2num(C{1}{31});
dxy       = str2num(C{1}{32});
dtheta    = str2num(C{1}{33});
jumpPrint = str2num(C{1}{35});
fclose(fid);


%%------------------ 0.3) Initialisation ------------------%%
%%---------------------------------------------------------%%

nTime = floor(Time/dt+.5);
if (shouldPlotEnd==1)
    jumpPrint = nTime;
end

switch(choiceDensity)
  case 1
    %% Density in x
    n_x = floor(Lx/dx);
    intX = linspace(0,Lx,n_x+1);
    %% parameters
    zMin = -1;
    zMax = 3;
  case 2
    %% Density 2D
    nXgrid = floor(Lx/dxy);
    nYgrid = floor(Ly/dxy);
    %%intX = dxy*((1:nXgrid)-.5);
    %%intY = dxy*((1:nYgrid)-.5);
    intX = dxy*(0:nXgrid);
    intY = dxy*(0:nYgrid);
    %% colormap hot
    l_red   = [ones(1,38)  , (26:-1:1)/26];
    l_green = [ones(1,13)   , (26:-1:1)/26 , zeros(1,25)];
    l_blue  = [(12:-1:1)/12 , zeros(1,52)];
    A_hot   = [l_red',l_green',l_blue'];
    %% cmin,cmax
    zMin = 0;
    rhoMean = N/(nXgrid*nYgrid)/dxy/dxy;
    zMax = 4*rhoMean;
  case 3
    %% Density in θ
    nTheta   = floor(2*pi/dtheta);
    intTheta = linspace(-pi,pi,nTheta+1);
    %% The theoretical curve    
    intTheta2 = linspace(-pi,pi,200);
    d_Theoric = d/sqrt(nu);
    C = 1/(2*pi*besseli(0,1/d_Theoric));
    %% parameters
    zMin = -.1;
    zMax = 1.3;
end

clf

%% load Binary data
function l=loadBinary(nameFile)
%% read the binary data in the file 'nameFile' in 'l'
    fid = fopen(nameFile, 'r');
    fseek(fid, 4, 'bof');
    l = fread(fid, 'double')';
    fclose(fid);
endfunction
function matrixData=loadBinary2D(nameFile,numberRow)
%% Read the binary data in the file 'nameFile' in 'l'.
%% We also need to precise the number of rows (numberRow).
    fid = fopen(nameFile, 'r');
    fseek(fid, 4, 'bof');
    matrixData = fread(fid,[numberRow,Inf],'double');
    fclose(fid);
endfunction


%%---------------------------------------------------------------%%
%%------------------------  la boucle  --------------------------%%
%%---------------------------------------------------------------%%
tic

for iTime=0:jumpPrint:nTime	% Warning: rm images/*    
    %%---------------    A) load data     ---------------%%
    %%---------------------------------------------------%%
    %extension = [num2str(iTime,'%09d'), '.udat'];
    extension = ['seed',num2str(nSeed),'_',num2str(iTime,'%09d'), '.udat'];
    switch(choiceDensity)
      case 1
        %% Density in x
        dens1Dx   = loadBinary(['../data/dens1Dx_' ,extension]);
        theta_x   = loadBinary(['../data/theta1Dx_',extension]);
        u1Dx      = loadBinary(['../data/u1Dx_'    ,extension]);
        v1Dx      = loadBinary(['../data/v1Dx_'    ,extension]);
        e1Dx      = 1 - u1Dx.^2 - v1Dx.^2;
      case 2
        %% Density 2D
        rho2D     = loadBinary2D(['../data/rho2D_',extension],nXgrid+1);
        u2D       = loadBinary2D(['../data/u2D_'  ,extension],nXgrid+1);
        v2D       = loadBinary2D(['../data/v2D_'  ,extension],nXgrid+1);
      case 3
        %% Density θ
        densTheta = loadBinary(['../data/densTheta_',extension]);
    end
    %%---------------    B) plot data     ---------------%%
    %%---------------------------------------------------%%
    switch(choiceDensity)
      case 1
        %% Density in x
        plot(intX,Lx*dens1Dx, 'linewidth',2,...
             intX,u1Dx,'-.', ...
             intX,v1Dx)
        xlabel('x')
        axis([0 Lx zMin zMax])
        title(['t = ',num2str(iTime*dt,'%10.2f')])
        h = legend('rho','u','v');
        if shouldSave==1                % saving?
            l = ['images/density1Dx_' , num2str(iTime,'%09d') , '.jpg'];
            print(sprintf(l))           % ,'-S800,800','-F:12'
        end
      case 2
        %% Density 2D
        clf
        hold on
        pcolor(intX,intY,rho2D');                   %imagesc(intX,intY,rho2D'); set(gca,'YDir','normal')
        h = quiver(intX,intY, lengthArrow2D*u2D',lengthArrow2D*v2D', ...
                   'linewidth',1,'autoScale','off');% imagesc: u,v inverse
        set (h, 'maxheadsize', .3);
        hold off
        %% deco
        title(['Density and velocity at  t = ',num2str(iTime*dt,'%10.2f')],'fontsize',40)
        colormap(A_hot);
        caxis([zMin zMax])
        c = colorbar;
        xlabel('x','fontsize',36);
        ylabel('y','fontsize',36);
        set(gca,'FontSize',30)
        axis([0 Lx 0 Ly],'equal')
        if shouldSave==1                 % saving?
            l = ['images/density2D_',num2str(iTime,'%09d'),'.png'];
            print(sprintf(l),'-S1000,800');
        end
      case 3
        %% Density θ
        if (isTheoricCurveTheta==0)
            plot(intTheta,densTheta,'o-','linewidth',2.5,'markersize',10);
            legend('Numeric');
        else
            %% theta average
            thetaAverage = atan2(sum(sin(intTheta).*densTheta),sum(cos(intTheta).*densTheta));
            plot(intTheta,densTheta,'o-','linewidth',2.5,'markersize',10, ...
                 intTheta,densTheta,'o','markersize',0, ...
                 intTheta2,C*exp(cos(intTheta2-thetaAverage)/d_Theoric),'r','linewidth',3);
            legend('Numeric',' ','Theory');
        end
        xlabel('direction \theta','Interpreter','tex','fontsize',24);
        axis([-pi pi -.1 1.1])
        set(gca,'xtick',[-pi -pi/2 0 pi/2 pi]);
        set(gca,'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},'Interpreter','tex')
        set(gca,'FontSize',24)
        title(['t = ',num2str(iTime*dt,'%10.2f')],'fontsize',24)
        if shouldSave==1                % saving
            l = ['images/densityTheta_' , num2str(iTime,'%09d') , '.png'];
            print(sprintf(l));
        end
    end
    pause(.01)		% small break (to see something...)
end

toc
%%---------------------------------------------------------------%%

if (shouldSave==1)
    %% we make a movie
    timeNow = clock;
    extension = [num2str(timeNow(1),'%04d'),...
                 num2str(timeNow(2),'%02d'),...
                 num2str(timeNow(3),'%02d'),'_',...
                 num2str(timeNow(4),'%02d'),'h',...
                 num2str(timeNow(5),'%02d')];
    switch choiceDensity
      case(1)
        name = ['videos/density1DxMicro_',extension,'.avi'];
        system(['mencoder ''mf://images/density1Dx_*.png'' -mf fps=25 -o ',name,' -ovc lavc -lavcopts vcodec=mpeg4']);
      case(2)
        name = ['../videos/density2DMicro_',extension,'.avi'];
        system(['mencoder ''mf://images/density2D_*.png'' -mf fps=25 -o ',name,' -ovc lavc -lavcopts vcodec=mpeg4']);
      case(3)
        name = ['../videos/densityTheta_',extension,'.avi'];
        system(['mencoder ''mf://images/densityTheta_*.png'' -mf fps=25 -o ',name,' -ovc lavc -lavcopts vcodec=mpeg4']);
    end
end

%%-- video
%% system(['mencoder ''mf://images/*.png'' -mf fps=15 -o videos/bidon.avi -ovc lavc -lavcopts vcodec=mpeg4']);
%%-- crop
%% for name in density2D_*; do convert $name -crop 1100x900+40+0 $name; done
