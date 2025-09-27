close all

global usecamera_g folder_g video_g

% base units: length in mm, mass in g, time in seconds
% derived units:
% [rho]   = g / mm^3                     = 10^6 kg/m^3
% [acc.]  = mm / s^2                     = 10^-3 m / s^2
% [force] = g * mm / s^2   = 10^-6 N     = \mu N
% [press] = 10^-6 N / mm^2 = N / m^2     = Pa
% [sigma] = 10^-6 N / mm   = 10^-3 N / m = mN/m

numtimesteps = 107;  % number of time steps (for usecamera_g == 1)
deltat = 1;          % duration of each time step [s] 
nsample = 120;        % number of points that are detected on the interface
nsample_fit = 40;    % number points on sample line
nmodes = 16;         % highest Chebychev modes for fit (max 24)
ndof = 30;           % number of DOFs for fit
deltarho = 0.225e-3; % density difference [10^6 kg/m^3]
grav = 9.807e3;      % gravitational acceleration [mm/s^2]
rneedle = 1.9/2;     % diameter of the needle (for pixel to length ratio)
crop = 1;            % crop the image before processing
smpl = 1;            % = 1: simple algorithm to find interface
                     % = 0: fit a sigmoidal function 
                     %      f(x) = a/(1+exp(-b*(x-c)))+d;
flp = 0;             % = 0 if needle enters from the left (bubble)
                     % = 1 if needle enters from the right (drop)
                     
% use calcST.m to determine the values of zneedle, zneck and sp!
                     
% set the values of zneedle and zneck (after cropping)
zneedle = 8; 
zneck = 72;     

% set the values for cropping the figure
sp = [8 15 562 757]; % [xmin, ymin, xmax, ymax]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist('zneedle','var') || ~exist('zneedle','var')
  error(['The variables zneedle and zneck must be known',...
    '(use calcST to determine these values!)']);
end

if ~exist('sp','var') && crop == 1
  error(['The variables sp must be known when crop == 1',...
    '(use calcST to determine these values!)']);
end

% store the numerical parameters in an array
num_param = [nsample, nsample_fit, nmodes, ndof];

% if using files: get all filenames in folder_g (ending on .png)
if usecamera_g == 0
  [allfilenames,alltimes] = getFileNames(folder_g);
  numsteps = length(allfilenames);
else
  numsteps = numtimesteps;
end

% initialize some data
dat = zeros(numsteps,6);

% create figure for real-time plotting of data
figure(1); hold on

for ii = 1:numsteps

  tic
  
  disp(['At step ',num2str(ii)]) 
 
  % grab the image from the camera or from a file
  if usecamera_g
    im = step(video_g); % from camera
    time = (ii-1)*deltat;
  else
    disp(['Loading file ',allfilenames{ii}]);
    im = imread([folder_g,allfilenames{ii}]); % from file
    
    % get the actual time of the current file
    time = (ii-1)*deltat;
  end

  % check if color image; if yes: convert to greyscale
  if length(size(im)) == 3
    im = rgb2gray(im);
  else
    im = im;
  end
  
  % flip image with needle to the right
    im = imrotate(im,90);
    
  % flip the image from right to left
  if flp == 1; im = flip(im,2); end
  
  % crop the image
  if crop == 1; im = im(sp(2):sp(4), sp(1): sp(3),:); end
  
  % store the previous coordinates
  if ii > 1; zz_old = zz; rr_old = rr; end
  
  % get the coordinates of the interface
  [zz,rr,~,~,im,out] = ...
    getshape(im,num_param,0,smpl,rneedle*2,0,zneedle,zneck);
  
   % check for invalid values
  num_inv = sum(isnan(zz))+sum(isinf(zz))+sum(isnan(rr))+sum(isinf(rr));
  
  % copy prevous shape if getshape unsuccesful
  if out ~= 1 || num_inv > 0 
    
    warning('unable to find shape: using previous shape')
    zz = zz_old; rr = rr_old;

  else

    % scale the lenghts with the radius of the needle
    rr = rr/rneedle; zz = zz/rneedle;

    % substract average r-value to move shape to z-axis
    rr = rr - mean(rr);

    % same defition as Nagel
    zz = -(zz - zz(end));
  
  end

  % get continuous s around full shape
  ds = sqrt((zz(2:end)-zz(1:end-1)).^2+(rr(2:end)-rr(1:end-1)).^2);
  ss = zeros(length(rr),1);
  for i=1:length(ds)
      ss(1+i) = sum(ds(1:i));
  end

  % fit the shape as Cheby basis functions and diff/int matrices
  [kappas_nag,kappap_nag,zzfit_nag,rrfit_nag,ssfit_nag, ...
    psifit_nag,diffmat_nag,intmat_nag] = nagel(zz,rr,ss,num_param,0);

  % find volume and area by numerical integration
  volume = pi*intmat_nag*(rrfit_nag.^2.*sin(psifit_nag));
  area = 2*pi*intmat_nag*rrfit_nag;

  % convert the dimensionless values back to dimensional
  volume = volume*rneedle^3;
  area = area*rneedle^2;

  fmt = '%9.3f'; % format for printing output

  disp(['Volume = ',num2str(volume,fmt),' mm^3'])
  disp(['Area = ',num2str(area,fmt),' mm^2'])

  % calculate the best fitting Laplace shape
  [st,press,rrlaplace,zzlaplace] = ...
    makeIso(zzfit_nag,rrfit_nag,psifit_nag,diffmat_nag);

  % convert the dimensionless values back to dimensional
  tension = st*deltarho*grav*rneedle^2;
  pcap = press*deltarho*grav*rneedle;

  disp(['Interfacial tension = ',num2str(tension,fmt),' mN/m'])
  disp(['Capillary press = ',num2str(pcap,fmt),' Pa'])
  
  % store the data
  dat(ii,:) = [ii, time, tension, pcap, volume, area];
      
  % plot the data
  figure(1)
  subplot(2,2,1)       % add first plot in 2 x 2 grid
  plot(dat(1:ii,2),dat(1:ii,3),'-r');
  xlabel('time (sec)')
  ylabel('interfacial tension (mN/m)')
  %xlim([0, numsteps*deltat]);    
  %ylim([0, 100]);

  subplot(2,2,2)       % add second plot in 2 x 2 grid
  plot(dat(1:ii,2),dat(1:ii,4),'-r');
  xlabel('time (sec)')
  ylabel('calc pressure (Pa)')
%   xlim([0, numsteps*deltat]);    
%   ylim([0, 100]);

  subplot(2,2,3)       % add third plot in 2 x 2 grid
  plot(dat(1:ii,2),dat(1:ii,5),'-r');
  xlabel('time (sec)')    
  ylabel('volume (mm^3)')    
%   xlim([0, numsteps*deltat]);
%   ylim([0, 100]);

  subplot(2,2,4)       % add fourth plot in 2 x 2 grid
  plot(dat(1:ii,2),dat(1:ii,6),'-r');
  xlabel('time (sec)')
  ylabel('area (mm^2)')
%   xlim([0, numsteps*deltat]);  
%   ylim([0, 100]);
  
  % write the image to a file when using camera
  if usecamera_g
    timestamp = datestr(now,'HH.MM.SS.FFF');
    filename = ['test_figures/',timestamp,'.png'];
    imwrite(im,filename);
  end
  
  toc % print the elapsed time

  % pause the remaining time of this time step
  if usecamera_g == 1
    elap = toc;
    if elap > toc
      warning('timestep smaller than time needed to run script!')
    else
      pause(deltat-elap); % pause the remain time
    end
  end
  
  disp(' ')

end
