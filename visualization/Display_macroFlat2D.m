#! /usr/bin/octave -qf

%%
%%--- Plot data computed with Vicsek (data must be in 'data')
%%

%% Parameters visualisation
choiceDensity = 2;		# 1: ρ(x),u(x)  2: ρ(x,y),u(x,y) 3: g(θ)
lengthArrow2D = .3;	        # for 2: ρ(x,y),u(x,y)
isTheoricCurveTheta = 0;	# for 3: g(θ)

shouldSave    = 0;
shouldPlotEnd = 0;


%%------------------ 0.1) Read parameters ------------------%%
%%----------------------------------------------------------%%
fid = fopen("../bin/PARAMETER_MicroVic_flat.txt");
for i=1:4
  temp = fgetl(fid);
endfor
N   = str2num(fgetl(fid));
temp = fgetl(fid);
c  = str2num(fgetl(fid));
nu = str2num(fgetl(fid));
d  = str2num(fgetl(fid));
rVision = str2num(fgetl(fid));
temp = fgetl(fid);
Lx = str2num(fgetl(fid));
Ly = str2num(fgetl(fid));
for i=1:4
  temp = fgetl(fid);
endfor
dt   = str2num(fgetl(fid));
Time = str2num(fgetl(fid));
for i=1:12
  temp = fgetl(fid);
endfor
dx = str2num(fgetl(fid));
dxy = str2num(fgetl(fid));
dtheta = str2num(fgetl(fid));
temp = fgetl(fid);
jumpPrint = str2num(fgetl(fid));
fclose(fid);


%%------------------ 0.3) Initialisation ------------------%%
%%---------------------------------------------------------%%

nTime = floor(Time/dt+.5);
if (shouldPlotEnd==1)
  jumpPrint = nTime;
endif

switch(choiceDensity)
  case 1
    %% Density in x
    n_x = floor(Lx/dx);
    intX = linspace(0,Lx,n_x+1);
    ## parameters
    zMin = -1;
    zMax = 3;
  case 2
    %% Density 2D
    nXgrid = floor(Lx/dxy);
    nYgrid = floor(Ly/dxy);
    intX = dxy*((1:nXgrid)-.5);
    intY = dxy*((1:nYgrid)-.5);
    ##intX = dxy*(0:nXgrid);
    ##intY = dxy*(0:nYgrid);
    %% colormap hot
    l_red   = [ones(1,38)  , (26:-1:1)/26];
    l_green = [ones(1,13)   , (26:-1:1)/26 , zeros(1,25)];
    l_blue  = [(12:-1:1)/12 , zeros(1,52)];
    A_hot   = [l_red',l_green',l_blue'];
    %% parameters
    zMin = 0;
    zMax = .05;
  case 3
    ## Density in θ
    nTheta   = floor(2*pi/dtheta);
    intTheta = linspace(-pi,pi,nTheta+1);
    ## The theoretical curve    
    intTheta2 = linspace(-pi,pi,200);
    d_Theoric = d/sqrt(nu);
    C = 1/(2*pi*besseli(0,1/d_Theoric));
    ## parameters
    zMin = -.1;
    zMax = 1.3;
endswitch

## init the figure
if (shouldSave==1)
  figure("visible","off")
  system("rm images/density*")
else
  figure("visible","on")
endif
clf

## load Binary data
function l=loadBinary(nameFile)
  ## read the binary data in the file 'nameFile' in 'l'
  fid = fopen(nameFile, 'r');
  fseek(fid, 4, 'bof');
  l = fread(fid, 'double')';
  fclose(fid);
endfunction
function matrixData=loadBinary2D(nameFile,numberRow)
  ## Read the binary data in the file 'nameFile' in 'l'.
  ## We also need to precise the number of rows (numberRow).
  fid = fopen(nameFile, 'r');
  fseek(fid, 4, 'bof');
  matrixData = fread(fid,[numberRow,Inf],'double');
  fclose(fid);
endfunction



%%---------------------------------------------------------------%%
%%---------------------------------------------------------------%%
%%------------------------  la boucle  --------------------------%%
%%---------------------------------------------------------------%%
%%---------------------------------------------------------------%%


tic

for iTime=0:jumpPrint:nTime	# Warning: rm images/*
  
  ##---------------------------------------------------##
  ##---------------    A) load data     ---------------##
  ##---------------------------------------------------##
  
  switch(choiceDensity)
    case 1
      %% Density in x
      %%-------------
      dens1Dx   = loadBinary(["../data/dens1Dx_",num2str(iTime,"%09d"), ".udat"]);
      theta_x   = loadBinary(["../data/theta1Dx_",num2str(iTime,"%09d"), ".udat"]);
      u1Dx      = loadBinary(["../data/u1Dx_",num2str(iTime,"%09d"), ".udat"]);
      v1Dx      = loadBinary(["../data/v1Dx_",num2str(iTime,"%09d"), ".udat"]);
      e1Dx      = 1 - u1Dx.^2 - v1Dx.^2;
    case 2
      %% Density 2D
      %%-----------
      rho2D     = loadBinary2D(["../data/rho2D_",num2str(iTime,"%09d"), ".udat"],nXgrid+1);
      u2D       = loadBinary2D(["../data/u2D_",num2str(iTime,"%09d"), ".udat"],nXgrid+1);
      v2D       = loadBinary2D(["../data/v2D_",num2str(iTime,"%09d"), ".udat"],nXgrid+1);
    case 3
      %% Density θ
      %%----------
      densTheta = loadBinary(["../data/densTheta_",num2str(iTime,"%09d"), ".udat"]);
  endswitch


  
  ##---------------------------------------------------##
  ##---------------    B) plot data     ---------------##
  ##---------------------------------------------------##

  switch(choiceDensity)
    case 1
      %% Density in x
      %%-------------
      plot(intX,Lx*dens1Dx,'linewidth',2,intX,u1Dx,'linewidth',2,intX,e1Dx,"-",'linewidth',2)
      xlabel("x",'fontsize',14)
      axis([0 Lx zMin zMax])
      title(["t = ",num2str(iTime*dt,"%10.2f")],'fontsize',14)
      legend("mass","theta","energy")
      %% save ?
      if shouldSave==1
	l = ["images/density1Dx_" , num2str(iTime,"%09d") , ".png"];
	print(sprintf(l),"-S800,800");
      endif
    case 2
      %% Density 2D
      %%-----------
      %% avoid saturation (ρ)
      rho2D = rho2D - (rho2D>zMax).*(rho2D-zMax);
      ##- %% ρ,u,v inside cells
      Kernel2 = .25*[1 1;1 1];
      f_u2D = rho2D .* u2D;
      f_v2D = rho2D .* v2D;
      rho2D_conv = conv2(rho2D,Kernel2);
      u2D_conv	 = conv2(f_u2D,Kernel2)./rho2D_conv;
      v2D_conv	 = conv2(f_v2D,Kernel2)./rho2D_conv;
      rho2D = rho2D_conv(2:(end-1),2:(end-1));
      u2D   = u2D_conv(2:(end-1),2:(end-1));
      v2D   = v2D_conv(2:(end-1),2:(end-1));
      ## plot
      clf
      hold on
      imagesc(intX,intY,rho2D');
      set(gca,'YDir','normal')
      h = quiver(intX,intY,lengthArrow2D*u2D',lengthArrow2D*v2D','linewidth',3,'autoScale','off'); # imagesc: u,v inverse
      set (h, "maxheadsize", .3);
      hold off
      ## deco
      title(["Density and velocity at  t = ",num2str(iTime*dt,"%10.2f")],'fontsize',14)
      axis([0 Lx 0 Ly],"equal")
      colormap(A_hot);
      caxis([zMin zMax])
      c = colorbar;
      xlabel("x",'fontsize',18);
      ylabel("y",'fontsize',18);
      %% save ?
      if shouldSave==1
	l = ["images/density2D_",num2str(iTime,"%09d"),".png"];
	##print(sprintf(l),"-S800,800");
	print(sprintf(l));
	close all
      endif
    case 3
      %% Density θ
      %%----------
      if (isTheoricCurveTheta==0)
	plot(intTheta,densTheta,"o-",'linewidth',2.5,'markersize',10);
	legend("Numeric");
      else
	## theta average
	thetaAverage = atan2(sum(sin(intTheta).*densTheta),sum(cos(intTheta).*densTheta));
	plot(intTheta,densTheta,"o-",'linewidth',2.5,'markersize',10,...
	     intTheta,densTheta,"o",'markersize',0,...
	     intTheta2,C*exp(cos(intTheta2-thetaAverage)/d_Theoric),"r",'linewidth',3);
	legend("Numeric"," ","Theory");
      endif
      xlabel('direction \theta','Interpreter','tex','fontsize',24);
      axis([-pi pi -.1 1.1])
      set(gca,'xtick',[-pi -pi/2 0 pi/2 pi]);
      set(gca,'xticklabel',{'-\pi','-\pi/2','0','\pi/2','\pi'},'Interpreter','tex')
      set(gca,'FontSize',24)
      title(["t = ",num2str(iTime*dt,"%10.2f")],'fontsize',24)
      %% save ?
      if shouldSave==1
	l = ["images/densityTheta_" , num2str(iTime,"%09d") , ".png"];
	print(sprintf(l));
      endif
      
  endswitch

  pause(.001)		% small break (to see something...)
endfor

%%---------------------------------------------------------------%%
%%---------------------------------------------------------------%%

toc

if (shouldSave==1)
  %% we make a movie
  timeNow = clock;
  extension = [num2str(timeNow(1),"%04d"),...
	       num2str(timeNow(2),"%02d"),...
	       num2str(timeNow(3),"%02d"),"_",...
	       num2str(timeNow(4),"%02d"),"h",...
	       num2str(timeNow(5),"%02d")];
  switch choiceDensity
    case(1)
      name = ["videos/density1DxMicro_",extension,".avi"];
      system(["mencoder ''mf://images/density1Dx_*.png'' -mf fps=25 -o ",name," -ovc lavc -lavcopts vcodec=mpeg4"]);
    case(2)
      name = ["../videos/density2DMicro_",extension,".avi"];
      system(["mencoder ''mf://images/density2D_*.png'' -mf fps=25 -o ",name," -ovc lavc -lavcopts vcodec=mpeg4"]);
    case(3)
      name = ["../videos/densityTheta_",extension,".avi"];
      system(["mencoder ''mf://images/densityTheta_*.png'' -mf fps=25 -o ",name," -ovc lavc -lavcopts vcodec=mpeg4"]);
  endswitch
else
  pause
endif

##-- video
## system(["mencoder ''mf://images/*.png'' -mf fps=15 -o videos/bidon.avi -ovc lavc -lavcopts vcodec=mpeg4"]);

##-- crop
## for name in density2D_*; do convert $name -crop 1100x900+40+0 $name; done
