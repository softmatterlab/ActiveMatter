function [centers,orientations]=ABP_detectEXP(image)

    % just ship with fixed parameters
    parameters.particleradius  =  6;
    parameters.radiustolerance =  3;
    parameters.picturemargin   = 10;
    
    centers = phase_hough_no_orientation(image); %Getting AP positions
    
    orientations = get_polarity(image,centers,parameters.particleradius); %Getting orientation
    
end

function [centers] = phase_hough_no_orientation(image, parameters)

    Image = double(image);
    radiusRange = parameters.particleradius-parameters.radiustolerance:0.5:parameters.particleradius+parameters.radiustolerance;
    %%%Phase Coding STARTs Here, adopted from chaccum.m, imfindcircles.m
    maxNumElemNHoodMat = 1e16; % Maximum number of elements in neighborhood matrix xc allowed, before memory chunking kicks in.
    % Calculate gradient
    hy = -fspecial('sobel');
    hx = hy';
    %Gx = imfilter(Image, hx, 'replicate','conv');
    Gx = conv2(Image, hx, 'same');
    %Gy = imfilter(Image, hy, 'replicate','conv');
    Gy = conv2(Image, hy, 'same');
    gradientImg = hypot(Gx, Gy);

    gradientImg(:,1)=0;
    gradientImg(:,end)=0;
    gradientImg(1,:)=0;
    gradientImg(end,:)=0;

    % Calculate gradient
    Gmax = max(gradientImg(:));
    edgeThresh = graythresh(gradientImg/Gmax);
    t = Gmax * cast(edgeThresh,'like',gradientImg);
    [Ey, Ex] = find(gradientImg > t); % Get edge pixels (positions)
    idxE = sub2ind(size(gradientImg), Ey, Ex); % Get edge pixels (Values)
    RR = -radiusRange;
    lnR = log(radiusRange);
    phi = ((lnR - lnR(1))/(lnR(end) - lnR(1))*2*pi) - pi;
    Opca = exp(sqrt(-1)*phi);
    w0 = Opca./(2*pi*radiusRange);
    xcStep = floor(maxNumElemNHoodMat/length(RR));
    lenE = length(Ex);
    [M, N] = size(Image);

    Ex_chunk = Ex(1:min(1+xcStep-1,lenE));
    Ey_chunk = Ey(1:min(1+xcStep-1,lenE));
    idxE_chunk = idxE(1:min(1+xcStep-1,lenE));
    xc = bsxfun(@plus, Ex_chunk, bsxfun(@times, -RR, Gx(idxE_chunk)./gradientImg(idxE_chunk))); % Eqns. 10.3 & 10.4 from Machine Vision by E. R. Davies
    yc = bsxfun(@plus, Ey_chunk, bsxfun(@times, -RR, Gy(idxE_chunk)./gradientImg(idxE_chunk)));
    xc = round(xc);
    yc = round(yc);
    w = repmat(w0, size(xc, 1), 1);
    %% Determine which edge pixel votes are within the image domain
    % Record which candidate center positions are inside the image rectangle.
    inside = (xc >= 1) & (xc <= N) & (yc >= 1) & (yc < M);
    % Keep rows that have at least one candidate position inside the domain.
    rows_to_keep = any(inside, 2);
    xc = xc(rows_to_keep,:);
    yc = yc(rows_to_keep,:);
    w = w(rows_to_keep,:);
    inside = inside(rows_to_keep,:);
    %% Accumulate the votes in the parameter plane
    xc = xc(inside); yc = yc(inside);
    Accumulator = abs(accumarray([yc(:), xc(:)], w(inside), [M, N]));
    %%%Phase Coding ENDs Here, ended with chaccum.m
    %Local Peak searching with own codes, replacing methods in imfindcircles
    %for speed
    Accumulator = FastlocalizePeak(Accumulator,parameters.particleradius,max(Accumulator(:))/10);
    Accumulator = FastSharpen(Accumulator,parameters.particleradius/2);

    Accumulator(:,1)=0; % prevents arrayoutofbounds errors in Peakfinder
    Accumulator(:,end)=0;
    Accumulator(1,:)=0;
    Accumulator(end,:)=0;

    [centers,~]=FastPeakFinder(Accumulator,max(Accumulator(:))*0.8,round(parameters.particleradius));
    centers(centers(:,1)<parameters.picturemargin | centers(:,1) > (N-parameters.picturemargin) ...
        | centers(:,2)<parameters.picturemargin | centers(:,2) > (M-parameters.picturemargin),:)=[];

end
 
function imsharp = FastSharpen(image,sigma)
    filter=exp(-(-sigma:sigma).^2/(2*sigma^2))/sqrt(2*pi*sigma^2);
    % To set pixels of intensity below centers.masklevel1old to zero
    %image=image.*(image>masklevel1); % faster than logical indexing (image(image<centers.masklevel1)=0)
    imsharp=conv2(image,filter'*filter,'same');
end

function imsharp = FastlocalizePeak(image,sigma,masklevel1)
    filter=exp(-(-sigma:sigma).^2/(2*sigma^2))/sqrt(2*pi*sigma^2);
    % To set pixels of intensity below centers.masklevel1old to zero
    image=image.*(image>masklevel1); % faster than logical indexing (image(image<centers.masklevel1)=0)
    conved=conv2(image,filter'*filter,'same')+0.0001;
    imsharp=image./conved;
end

function [peaks,intensity]=FastPeakFinder(image,masklevel1,sz)
    [nr,nc]=size(image);
    [ym, xm]=find(image>masklevel1);

    excluded=find(xm<1 | xm > (nc-1) | ym<1 | ym > (nr-1));
    xm(excluded)=[];
    ym(excluded)=[];
    disp(find(xm<2))
    Nm=length(xm);
    ind=zeros(Nm,1);
    x0=repmat(-1:1,3,1);
    x0=reshape(x0,1,9);
    y0=repmat(-1:1,1,3);

    for i=1:Nm    
        region=image(y0+ym(i)+nr*(x0+xm(i)-1)); % 1px thick region centered at pixel of interest
        ind(i)=sum((region(5)-region)>=0); % number of pixels of intensity smaller than or equal to the central one
    end

    xpk=xm(ind>8);
    ypk=ym(ind>8);

    %if size is specified, then get ride of pks within size of boundary
    if nargin==3 && numel(xpk)>0
       %throw out all pks within sz of boundary;
        excluded=find(xpk<sz | xpk > (nc-sz) | ypk<sz | ypk > (nr-sz));
        xpk(excluded)=[];
        ypk(excluded)=[];
    end

    %%%% Method 3
    % Fastest method so far (see below for other trials): Fairly optimised
    % Matrix method (no concatenation) +linear indexing + no averaging. Last operation costly but does not scale with number of peaks.
    impk=zeros(nr,nc);

    % Image only with detected peaks
    impk(ypk+nr*(xpk-1))=image(ypk+nr*(xpk-1));

    Npk=length(xpk);
    rsz=floor(sz/2);

    nn=2*rsz+1;
    xr=repmat(-rsz:rsz,nn,1);
    xr=reshape(xr,1,numel(xr));
    yr=repmat(-rsz:rsz,1,nn);
    for i=1:Npk   
        xi=xpk(i); yi=ypk(i); % create vraible to avoid repetitive access
        region=impk(yr+yi+nr*(xr+xi-1)); % cropped image centered at peak of interest
        m=find(region==max(max(region)));
        impk(yr+yi+nr*(xr+xi-1))=0;
        lm=m(1);
        xm=ceil(lm/nn)-(rsz+1); % check indexing here
        ym=mod(lm,nn)-(rsz+1); % check indexing here
        impk(ym+yi+nr*(xm+xi-1))=region(m(1)); % check indexing here
    end

    [rpk,cpk,intensity]=find(impk);
    peaks=[cpk,rpk];
end

function orientations=get_polarity(Image,centers,radius)
rcrop=ceil(radius);
polarities=zeros(size(centers));
% MaxMin=zeros(size(centers));
CedgeN=zeros(size(centers));

I=repmat(-rcrop:rcrop,2*rcrop+1,1); % indices of pixels on cropped image
J=repmat((-rcrop:rcrop)',1,2*rcrop+1); % indices of pixels on cropped image
Disk=sqrt(I.^2+J.^2)<radius;
X=I(Disk);
Y=J(Disk);

for index=1:size(centers,1)
    %%%%%%%%%%%% PART 1.    %%%%%%%%%%%%
    % Cropping the image around the particle's centre of mass
    %tic
    CM=centers(index,:); %center of mass in pixels
    cropped=Image((-rcrop:rcrop)+CM(2),(-rcrop:rcrop)+CM(1)); % cropped image
    
    % Forming a cloud of points representing intensity within radius
    Z=double(cropped(Disk));
    
    %%%%%%%%%%%% PART 2.    %%%%%%%%%%%%
    
    % Dark edge thresholding + centre of mass
    edge=find(Z>sum(Z)/numel(Z));
    Xedge=X(edge); Yedge=Y(edge); Zedge=Z(edge);
    Cedge=[Xedge'*Zedge,Yedge'*Zedge]/sum(Zedge);
    Rg=sqrt(((Xedge-Cedge(1)).^2+(Yedge-Cedge(2)).^2)'*Zedge/sum(Zedge));

    % Selection of light gradient inside the particle between the two caps.
    % Region centred at centre of the particle and of radius Rg.
    grad=find((X.^2+Y.^2)<Rg^2); % particles in light gradient
    Xgrad=X(grad); Ygrad=Y(grad); Zgrad=Z(grad);
    Cgrad=sum(Zgrad)/numel(Zgrad); % centre of mass of light gradient
    
    %%%%%%%%%%%% PART 3.    %%%%%%%%%%%%
    
    % Normal to the gradient plane using covariance matrix method
    %C=cov([X(grad)-Cgrad(1),Y(grad)-Cgrad(2),Z(grad)-Cgrad(3)]); %
    %covariance matrix ("cov" slow)
    %Covariance matrix built term by term is faster than using "cov". Average is replaced by
    %sum()/numel()
    Zgrad=Zgrad-Cgrad;
    N=numel(Xgrad)-1;
    C(1,2)=Xgrad'*Ygrad/N-sum(Xgrad)*sum(Ygrad)/N^2;
    C(1,3)=Xgrad'*Zgrad/N-sum(Xgrad)*sum(Zgrad)/N^2;
    C(2,3)=Ygrad'*Zgrad/N-sum(Ygrad)*sum(Zgrad)/N^2;
    C(1,1)=Xgrad'*Xgrad/N-sum(Xgrad)*sum(Xgrad)/N^2;
    C(2,2)=Ygrad'*Ygrad/N-sum(Ygrad)*sum(Ygrad)/N^2;
    C(3,3)=Zgrad'*Zgrad/N-sum(Zgrad)*sum(Zgrad)/N^2;
    C(2,1)=C(1,2); C(3,1)=C(1,3); C(3,2)=C(2,3);
    
    [V,~]=eig(C);
    
    %%%%%%%%%%%% PART 4.    %%%%%%%%%%%%
    
    P=V(:,1)-V(3,1)*[0;0;1]; % projection of first eigenvector onto x-y plane = polarity
    polarities(index,1:2)=[P(1),P(2)];
    
    CedgeN(index,:)=Cedge;
end

% Correct orientation is aligned with the max-min vector
normp=sqrt(polarities(:,1).^2+polarities(:,2).^2);
polarities(:,1)=polarities(:,1)./normp;
polarities(:,2)=polarities(:,2)./normp;
DotProd=polarities(:,1).*CedgeN(:,1)+polarities(:,2).*CedgeN(:,2);
% DotProd=polarities(:,1).*MaxMin(:,1)+polarities(:,2).*MaxMin(:,2);
polarities(DotProd<0,1:2)=-polarities(DotProd<0,1:2);


orientations=acos(polarities(:,1));
orientations(polarities(:,2)<0)=2*pi-orientations(polarities(:,2)<0);

end





