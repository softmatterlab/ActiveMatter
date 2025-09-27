close all

global usecamera_g filename_g video_g

% base units: length in mm, mass in g, time in seconds
% derived units:
% [rho]   = g / mm^3                     = 10^6 kg/m^3
% [acc.]  = mm / s^2                     = 10^-3 m / s^2
% [force] = g * mm / s^2   = 10^-6 N     = \mu N
% [press] = 10^-6 N / mm^2 = N / m^2     = Pa
% [sigma] = 10^-6 N / mm   = 10^-3 N / m = mN/

nsample = 100;        % number of points that are detected on the interface
nsample_fit = 80;    % number points on sample line
nmodes = 16;         % highest Chebychev modes for fit (max 24)
ndof = 30;           % number of DOFs for fit
deltarho = 0.225e-3; % 1e-3;     % density difference [10^6 kg/m^3]
grav = 9.807e3;      % gravitational acceleration [mm/s^2]
rneedle = 1.9/2;     % diameter of the needle (for pixel to length ratio)
crop = 1;            % crop the image before processing
smpl = 1;            % = 1: simple algorithm to find interface
                     % = 0: fit a sigmoidal function 
                     %      f(x) = a/(1+exp(-b*(x-c)))+d;
flp = 0;             % = 0 if needle enters from the left (bubble)
                     % = 1 if needle enters from the right (drop)
                     
% if these are known: set the values of zneedle and zneck (after cropping)
% zneedle = 68; zneck = 152;     

% if these are known: set the values for cropping the figure
% sp = [719 254 1870 1249]; % [xmin, ymin, xmax, ymax]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic

% store the numerical parameters in an array
num_param = [nsample, nsample_fit, nmodes, ndof];

% grab the image from the camera or from a file
if usecamera_g
  im = step(video_g); % from camera
else
  im = imread(filename_g); % from file
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

% crop the image to speed up the calculation
if crop
  
  if ~exist('sp','var')
    % display image with true aspect ratio
    figure; imshow(im); hold on; axis on; axis image
    p = ginput(2); 

    % get the x and y corner coordinates as integers
    sp(1) = min(floor(p(1)), floor(p(2))); %xmin
    sp(2) = min(floor(p(3)), floor(p(4))); %ymin
    sp(3) = max(ceil(p(1)), ceil(p(2)));   %xmax
    sp(4) = max(ceil(p(3)), ceil(p(4)));   %ymax

  end
  
  % crop the image
  im = im(sp(2):sp(4), sp(1): sp(3),:);    

end

% get the coordinates of the interface
if ~exist('zneedle','var')
  [zz,rr,zneedle,zneck,im,out] = getshape(im,num_param,1,smpl,rneedle*2,0);
else
  [zz,rr,~,~,im,out] = ...
    getshape(im,num_param,1,smpl,rneedle*2,0,zneedle,zneck);
end

% check for invalid values
num_inv = sum(isnan(zz))+sum(isinf(zz))+sum(isnan(rr))+sum(isinf(rr));

% check if the drop/bubble shape was found
if out ~= 1 || num_inv > 0 
  error('unable to find shape')
end

% scale the lenghts with the radius of the needle
rr = rr/rneedle; zz = zz/rneedle;

% substract average r-value to move shape to z-axis
rr = rr - mean(rr);

% same defition as Nagel
zz = -(zz - zz(end));

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

% display value of variable for calcST_time
disp(['zneedle= ', num2str(zneedle)])
disp(['zneck= ', num2str(zneck)])
disp(['sp= ', num2str(sp)])

% plot the experimental points and the best fitting Laplace shape
figure; scatter(rr,zz); hold on; 
plot(rrlaplace,zzlaplace,'r','LineWidth',3);
plot(-rrlaplace,zzlaplace,'r','LineWidth',3);
daspect([1 1 1]);
legend('exp.','Laplace fit');

% write the image to a file when using camera
if usecamera_g
  timestamp = datestr(now,'HH.MM.SS.FFF');
  filename = ['test_figures/',timestamp,'.png'];
  imwrite(im,filename);
end

toc
