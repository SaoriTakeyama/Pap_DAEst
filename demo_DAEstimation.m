%% This code 

clear
close all

%% image loading
n_bands = 14;   % band number
prec_byte = 4;
hdsz = 4;       % Requesting size for header data
offs_sz = 183;     % Offset for image size data
offs_thumb = 192;
offs_eof = 1568;

fname = ''; %Please input a MS image name of a Papanicolaou stain specimen
fd = fopen( fname, 'rb' );
fseek( fd, offs_sz, 'bof' );
hdata = fread( fd, hdsz, 'uint32' );

h = hdata(1);
w = hdata(2);
fprintf('%d bands, %d x %d, %d bits\n', n_bands, w, h, prec_byte*8 );

imsize = h * w * n_bands;
imbytes = imsize * prec_byte;

info = dir(fname);
fsize = info.bytes;
hsize = fsize - imbytes;
fprintf('File size = %d, Image size = %d, Header size = %d\n', fsize, imbytes, hsize );

frewind(fd);
fseek( fd, offs_thumb, 'bof' );
[~, cnt] = fread( fd, [w h], 'uint32' );
fprintf('%d data read (Request: %d).\n', cnt, w*h );

frewind(fd);
fseek( fd, hsize-offs_eof, 'bof' );
msd = fread( fd, imsize, 'float32' );
msim = reshape( msd, [n_bands w h]);

fclose( fd );

white = csvread("app_func_data/white_absorb.csv"); % load white data
msim = msim - repmat(reshape(white,14,1),[1,size(msim,[2,3])]);

bandimg2 = zeros(h,w,14);
for i = 1:14
   bandimg2(:,:,i) = squeeze(msim(i,:,:))';
end

xyz401 = csvread("app_func_data/xyz.csv"); % load the tranformation matrix from MS image (14 bands) to XYZ
xyz14 = xyz401(61:20:321,2:4)';
bandimg = reshape(bandimg2,w*h,14); 

msim2 = 10.^(-msim);
bandimg4 = zeros(h,w,14);
for i = 1:14
   bandimg4(:,:,i) = squeeze(msim2(i,:,:))';
end
bandimg3 = reshape(bandimg4,w*h,14);

gxyz = horzcat(sum(bandimg3.*xyz14(1,:),2),sum(bandimg3.*xyz14(2,:),2),sum(bandimg3.*xyz14(3,:),2));

ref = ones(1,14);
refilly = ref.*xyz14(2,:);
refillx = ref.*xyz14(1,:);
refillz = ref.*xyz14(3,:);
refillysum = sum(refilly);
refillxsum = sum(refillx);
refillzsum = sum(refillz);
k = 1.0/refillysum;
kx = 1.0/refillxsum;
kz = 1.0/refillzsum;

Xw = 0.3127;
Yw = 0.3290;
Zw = 1-Xw-Yw;
Xw = Xw*2.4;
Yw = Yw*2.4;
Zw = Zw*2.4;
gxyz_norm(:,2) = Yw.*k.*gxyz(:,2);
gxyz_norm(:,1) = Xw.*kx.*gxyz(:,1);
gxyz_norm(:,3) = Zw.*kz.*gxyz(:,3);

gxyz_norm(gxyz_norm>1) = 1.0;

xyz2srgbmat = [3.2404542 -1.5371385 -0.4985314;...
         -0.9692660  1.8520108  0.0415560;...
          0.0556434 -0.2040259  1.0572276];              

sRGB = (xyz2srgbmat*gxyz_norm')';

bgX = gammahosei(sRGB,h,w); % Gamma correction
bgX2 = reshape(bgX,h,w,3);

figure
imshow(bgX2)
save_name = ['result/',fname(1:end-4),'_RGB.png'];
imwrite(bgX2, save_name);

Ha = csvread("app_func_data/Ha_0823.csv"); %load the matrix H (premeasurement)

H14 = Ha(1:14,3:6)';
H = H14;

Hinv = pinv(H);
C2 = bandimg*Hinv;
C = reshape(C2,[h,w,4]);
C(C<0) = 0;

eosin = C(:,:,1);
hema = C(:,:,2);
LG = C(:,:,3);
OG = C(:,:,4);

%background selection
disp('Please select background area.')
rect = getrect;
X = rect(1,1);
Y = rect(1,2);
width = rect(1,3);
height = rect(1,4);
X = round(X);
Y = round(Y);
if X < 0
    X = 0;
end
if Y < 0
    Y = 0;
end
width = round(width);
height = round(height);
X2 = X + width;
Y2 = Y + height;
if X2 > w
    width = w - X;
    X2 = w;
end
if Y2 > h
    height = h - Y;
    Y2 = h;
end
bgX2(Y,X:X2,:) = 0;
bgX2(Y2,X:X2,:) = 0;
bgX2(Y:Y2,X,:) = 0;
bgX2(Y:Y2,X2,:) = 0;

M1 = eosin(Y:Y2, X:X2);
r = reshape(M1,[(width+1)*(height+1),1]);
triE = prctile(r,99.99);
eosin(eosin<triE) = 0.0;

M2 = hema(Y:Y2, X:X2);
r = reshape(M2,[(width+1)*(height+1),1]);
triH = prctile(r,99.99);
hema(hema<triH) = 0.0;

M3 = LG(Y:Y2, X:X2);
r = reshape(M3,[(width+1)*(height+1),1]);
triL = prctile(r,99.99);
LG(LG<triL) = 0.0;

M4 = OG(Y:Y2, X:X2);
r = reshape(M4,[(width+1)*(height+1),1]);
triO = prctile(r,99.99);
OG(OG<triO) = 0.0;

eosin = eosin*0.499350650005020;
hema = hema*0.529574866150920;
LG = LG*0.392146052341812;
OG = OG*0.823913588791684;

r = reshape(eosin,[w*h,1]);
eosinmax = prctile(r,99.99);
if eosinmax < 1
    eosinmax = 1;
end

r = reshape(hema,[w*h,1]);
hemax = prctile(r,99.99);
if hemax < 1
    hemax = 1;
end

r = reshape(LG,[w*h,1]);
LGmax = prctile(r,99.99);
if LGmax < 1
    LGmax = 1;
end

r = reshape(OG,[w*h,1]);
OGmax = prctile(r,99.9);
if OGmax < 1
    OGmax = 1;
end

low=0;
figure
subplot(2,2,1)
imshow(eosin,[low eosinmax]);
title('eosin');
colorbar;
subplot(2,2,2)
imshow(hema,[low hemax]);
title('hema');
colorbar;
subplot(2,2,3)
imshow(LG,[low LGmax]);
title('LG');
colorbar;
subplot(2,2,4)
imshow(OG,[low OGmax]);
title('OG');
colorbar;

DA_img = cat(3,hema,eosin,LG,OG);
color_img = bgX2;
save('result/DA_data.mat','DA_img','color_img')