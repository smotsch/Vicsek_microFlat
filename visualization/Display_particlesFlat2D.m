%%
%%--- Plot data computed with MicroVic_flat.
%%

%% Parameters visualisation
lengthArrow = .1;		% the arrows
headSize    = .1;

shouldSave  = 0;
shouldPlotEnd = 0;


%%------------------ 0.1) Read parameters ------------------%%
%%----------------------------------------------------------%%
fid = fopen('../bin/PARAMETER_MicroVic_flat.txt');
C   = textscan(fid, '%s','delimiter', '\n');
Lx        = str2num(C{1}{12});
Ly        = str2num(C{1}{13});
BC        = str2num(C{1}{18});
dt        = str2num(C{1}{18});
Time      = str2num(C{1}{19});
stringTrajectorySave = C{1}{28};
dx        = str2num(C{1}{32});
dxy       = str2num(C{1}{33});
dtheta    = str2num(C{1}{34});
jumpPrint = str2num(C{1}{36});
fclose(fid);

%%------------------ 0.2) Initialisation ------------------%%
%%---------------------------------------------------------%%
nTime = floor(Time/dt+.5);
%% mini fix
if (stringTrajectorySave=='T' || stringTrajectorySave=='1')
  isTrajectorySave = 1;
else
  isTrajectorySave = 0;
end
if (isTrajectorySave==0)
  jumpPrint = nTime;
end
if (shouldPlotEnd==1)
  jumpPrint = nTime;
end
if (BC==4)
  a = Lx/2;
  b = Ly/2;
  t = linspace(0,2*pi,100);
end


%% load Binary data
function l=loadBinary(nameFile)
  %% read the binary data in the file 'nameFile' in 'l'
  fid = fopen(nameFile, 'r');
  fseek(fid, 4, 'bof');
  l = fread(fid, 'double')';
  fclose(fid);
endfunction

%% init the figure
if (shouldSave==1)
  figure('visible','off')
  system('rm images/particles*')
else
  figure('visible','on')
end
clf


%%---------------------------------------------------------------%%
%%---------------------------------------------------------------%%
%%------------------------  la boucle  --------------------------%%
%%---------------------------------------------------------------%%
%%---------------------------------------------------------------%%

tic

for iStep=0:jumpPrint:nTime	% Warning: rm images/*

  %%---- A) load data ----%%
  particleX     = loadBinary(['../data/particleX_' , num2str(iStep,'%09d'),'.udat']);
  particleY     = loadBinary(['../data/particleY_' , num2str(iStep,'%09d'),'.udat']);
  particleTheta = loadBinary(['../data/particleTheta_' , num2str(iStep,'%09d'),'.udat']);
  %%---- B) plot data ----%%
  %% cross or arrow
  if (lengthArrow==0)
    plot(particleX,particleY,'+b')
  else
    h = quiver(particleX,particleY,lengthArrow*cos(particleTheta),lengthArrow*sin(particleTheta),'linewidth',2,'AutoScale','off');
    set(h,'maxheadsize', headSize);
  end
  if (BC==5)
    hold on
    plot(a+a*cos(t),b+b*sin(t),'linewidth',2)
    hold off
  end
  %% deco
  legend('off')
  axis([0 Lx 0 Ly],'equal');
  title(['Particles (N = ',num2str(length(particleX)),')  at  t = ',num2str(iStep*dt,'%10.2f')],'fontsize',24)
  xlabel('x','fontsize',26)
  ylabel('y','fontsize',26)
  set(gca,'FontSize',18)
  %% save ?
  if (shouldSave==1)
    l = ['images/particles_' , num2str(iStep,'%09d') , '.png'];
    print(sprintf(l))
  else
    pause(.01)
  end
end

toc

%%---------------------------------------------------------------%%
%%---------------------------------------------------------------%%

if (shouldSave==1)
  %% we make a movie
  timeNow = clock;
  extension = [num2str(timeNow(1),'%04d'),...
	       num2str(timeNow(2),'%02d'),...
	       num2str(timeNow(3),'%02d'),'_',...
	       num2str(timeNow(4),'%02d'),'h',...
	       num2str(timeNow(5),'%02d')];
  %% particles alone    
  name = ['videos/particles_',extension,'.avi'];
  system(['mencoder ''mf://images/particles_*.png'' -mf fps=20 -o ',name,' -ovc copy']);
end

%%-- video
%% mencoder 'mf://*.png' -mf fps=25 -o ../output.avi -ovc lavc -lavcopts vcodec=mpeg4

%%-- crop
%% for name in particles_*; do convert $name -crop 900x900+150+0 $name; done
%% for name in densityTheta_*; do convert $name -crop 1100x900+40+0 $name; done

%%-- melange
%% for i in `seq -f '%09g' 0 2 2000`; do
%%   nameA='particles_'$i'.png';
%%   nameB='densityTheta_'$i'.png';
%%   nameC='melange_'$i'.png';
%%   concatenatePng $nameA $nameB $nameC;
%% done

%% for i in `seq -f '%09g' 0 4 5000`; do
%%   nameA='melange_'$i'.png';
%%   cp $nameA '../image2/'$nameA;
%% done


%mencoder 'mf://melange*.png' -mf fps=15 -o ../output.avi -ovc lavc -lavcopts vcodec=mpeg4


break

%% estimation of the radial distribution
l = max([abs(particleX-5);abs(particleY-5)]);
intR = .25:.5:4.75;
temp = 0:.5:5;
dl = 4*(temp(2:(end)).^2 - temp(1:(end-1)).^2);
z = hist(l,intR)./dl;
plot(intR,z)
