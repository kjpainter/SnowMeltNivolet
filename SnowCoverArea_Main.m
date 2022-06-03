close all; clear all;
global mask mask2 incidence elevation n m cells nbors data prop frac means

greenwhite(1,1:3) = [0 0.2 0];
greenwhite(2,1:3) = [1 1 1];

%% 
% read data and set up two masks for the domain extent
mask = imread('Binary_Mask.tif');  
mask2= mask;
mask2(mask2<0) = NaN;   %substitute -3.4028235e+38 with NaN
mask(mask<0) = 0;       %substitute -3.4028235e+38 with 0

cells = sum(sum(mask)); %number of pixels in the region of interest (ROI)
boundarylength=(sum(sum(abs(diff(mask))))+sum(sum(abs(diff(mask'))))); %length of domain boundary
[m n] = size(mask); %rows and columns of mask

%calculate the number of neighbouring (domain) patches for each (domain) patch 
nbors = zeros(m,n);
nbors(2:m-1,2:n-1)=mask(1:m-2,2:n-1)+mask(3:m,2:n-1)+mask(2:m-1,1:n-2)+mask(2:m-1,3:n);
nbors(1,2:n-1)=mask(2,2:n-1)+mask(1,1:n-2)+mask(1,3:n);
nbors(m,2:n-1)=mask(m-1,2:n-1)+mask(m,1:n-2)+mask(m,3:n);
nbors(2:m-1,1)=mask(1:m-2,1)+mask(3:m,1)+mask(2:m-1,2);
nbors(2:m-1,n)=mask(1:m-2,n)+mask(3:m,n)+mask(2:m-1,n-1);
nbors(1,1)=mask(1,2)+mask(2,1);
nbors(1,n)=mask(1,n-1)+mask(2,n);
nbors(m,1)=mask(m-1,1)+mask(m,2);
nbors(m,n)=mask(m-1,n)+mask(m,n-1);

% read and convert incidence angle data
a = imread('INCIDENCEANGLE14m_Nivolet_DOY139_4326_filled.tif');
[r c val] = find(isnan(a)==false); 
incidence = -1*ones(m,n); 
for i = 1:length(r)
    incidence(r(i),c(i)) = abs(1-a(r(i),c(i))/90); 
end
[r c val] = find(isnan(a)==true & mask >0);
for i = 1:length(r)
    incidence(r(i),c(i)) = 0;
end  

% read and covert elevation data
b = imread('DEM14m_Nivolet_4326_filled.tif'); 
b(b<0) = NaN;
hmax = max(max(b)); %elevation raster maximum
hmin = min(min(b)); %elevation raster minimum
elevation = -1*ones(m,n);
for i = 1:m
    for j = 1:n
        if (isnan(b(i,j)) == false)     
            elevation(i,j) = (hmax-b(i,j))/(hmax-hmin);
        end
    end
end

% calculate domain mean values for elevation and incidence angle
esum = 0;
[ii jj] = find(incidence>=0);
for i = 1:length(ii)
    esum = esum+incidence(ii(i),jj(i));
end
means(1) = esum/length(ii);
hsum = 0;
[ii jj] = find(elevation>=0);
for i = 1:length(ii)
    hsum = hsum+elevation(ii(i),jj(i));
end
means(2) = hsum/length(ii);

% clear all temporary variables
clear a b c ii jj i j r read val hmax hmin hsum esum

%% 
% load in the snowcover masks for fitting
data(:,:,1)= imread('SnowCover20180525.tif'); 
data(:,:,2)= imread('SnowCover20180614.tif');
data(:,:,3)= imread('SnowCover20180619.tif'); 
data(:,:,4)= imread('SnowCover20180709.tif'); 
data(:,:,5)= imread('SnowCover20180729.tif');
dimdata = size(data);  
if (length(dimdata)<3) 
    dimdata = 1; 
else
    dimdata = dimdata(:,3); % n0layers
end

% measure fraction of snow cover for each snow mask.
for i = 1:dimdata
    prop(i)= sum(sum(data(:,:,i)));    
    l1 = find(abs(diff(data(:,:,i).*mask2))==1);  
    l2 = find(abs(diff((data(:,:,i).*mask2)'))==1);    
    frac(i)=(length(l1)+length(l2))/cells;
end

%% 
% set up parameter ranges and start the simulation loop
% Each model parameter ranges from lb to ub
% [rho alpha beta gamms pp qq rr]
lb = [2 0 0 0 0 0 0]; 
ub = [10 9 9 9 3 3 3];

% Number of sample parameter sets to test, number of runs for each sample
% set. Parameter sets randomly selected via Latin Hypercube Sampling
% default below is 25 sample sets, with one run of each. 
nsamples = 25;
nruns = 1;
samplepars = lhsdesign(nsamples,length(lb));
drawon = true;
for i = 1:nsamples
    xin = (lb+(ub-lb).*samplepars(i,:));

% If you specify the following parameter set, this represents the "best
% fit" parameter set from 5000 sample sets, averaged over 5 runs of each
% set
%     xin = [9.174 3.725 1.365 1.904 2.921 2.874 2.160];
    jj = 1;
    while jj < nruns+0.5                    
% call the stochastic solver, returning output corresponding to each snow cover mask and measures for computing error        
        [output runmeasure(i,jj,:)] = StochasticSolver(xin(1),xin(2),xin(3),xin(4),xin(5),xin(6),xin(7));
        lambda = 0.75;
        runerror(i,jj)=(lambda*sum(runmeasure(i,jj,1:2:end))+(1-lambda)*sum(runmeasure(i,jj,2:2:end)))/length(prop);
% if simulation exited with high error, give up on this set.
        if (runerror(i,jj) >0.9)
            runmeasure(i,jj:nruns,1:10) = 1;
            runerror(i,jj:nruns) = 1;
            jj = nruns;
        end
% drawdata        
        if (drawon)
            figure(1)
            set(1,'paperpositionmode','auto','position',[200 300 900 180])
            for ll = 1:5
                subplot('position',[0.005+(ll-1)*0.2 0.005 0.19 0.99])
                expand = ones(m+2,n+2);
                expand(2:end-1,2:end-1)=flipud(squeeze(output(:,:,ll))-mask+1);
                expandmask = zeros(m+2,n+2);
                expandmask(2:end-1,2:end-1) = flipud(mask);
                pcolor(expand),caxis([0 1]),colormap(greenwhite)
                set(gca,'xtick',[],'ytick',[]),shading flat, axis image
                axis off
                hold on
                contour(expandmask,[0.5 0.5],'k-','linewidth',2),axis image
                hold off
            end
        end
        optpars(i,:) = xin;
        fprintf('Sample set %d, Run %d, Error %5.3f, Pars %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n',i,jj,runerror(i,jj),xin) 
        jj = jj+1;
    end
end

%% 
% At the end of the simulation, we run and plot the best found set from the
% parameter space sampling above.

runplotbest = true;
if (runplotbest)
    [vals index] = sort(mean(runerror,2));
    xin = (lb+(ub-lb).*samplepars(index(1),:));
    [outputbest runbest] = StochasticSolver(xin(1),xin(2),xin(3),xin(4),xin(5),xin(6),xin(7));
    runerrorbest=(lambda*sum(runbest(1:2:end))+(1-lambda)*sum(runbest(2:2:end)))/length(prop);
    fprintf('Best sample set %d, Error %5.3f, Pars %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n',index(1),runerrorbest,xin) 
    
% drawdata        
    if (drawon)
        figure(1)
        set(1,'paperpositionmode','auto','position',[200 300 900 180])
        for ll = 1:5
            subplot('position',[0.005+(ll-1)*0.2 0.005 0.19 0.99])
            expand = ones(m+2,n+2);
            expand(2:end-1,2:end-1)=flipud(squeeze(outputbest(:,:,ll))-mask+1);
            expandmask = zeros(m+2,n+2);
            expandmask(2:end-1,2:end-1) = flipud(mask);
            pcolor(expand),caxis([0 1]),colormap(greenwhite)
            set(gca,'xtick',[],'ytick',[]),shading flat, axis image
            axis off
            hold on
            contour(expandmask,[0.5 0.5],'k-','linewidth',2),axis image
            hold off
        end
    end    
end







