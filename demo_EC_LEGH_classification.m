clear
close all

load 'result/DA_data.mat'
mode = 'Fischer_2'; % select the method solving the decision boundaly: 'Fischer_2', or 'Fischer_3'

% Data generation
[rows, cols, bands] = size(DA_img);
block_img = zeros(floor(rows/10),floor(cols/10),bands);
for n = 1:floor(rows/10)
    for m = 1:floor(cols/10)
        patch = DA_img((n-1)*10+1:(n-1)*10+10,(m-1)*10+1:(m-1)*10+10,:);
        H_area = patch(:,:,1) ~= 0;
        H_area = sum(H_area(:));
        S_exist = sum(patch,3) == 0;
        S_exist = sum(S_exist(:));
        patch = sum(sum(patch,1),2)/100;
        block_img(n,m,:) = patch;
        if (H_area >= 50) || (S_exist > 80)
            block_img(n,m,1) = 1;
        else
            block_img(n,m,1) = 0;
        end
    end
end
Label = block_img(:,:,1);
test_data = block_img(:,:,2:4);
test_data = reshape(test_data,[n*m,3]);

% classify the data into EC or LEGH
switch mode
    case 'manual'
        load 'app_func_data/Manual.mat'
        Est = @(x)(x(3,:)>0.005) .* (x(3,:)> ab(1)*x(1,:)+ab(2));
        Est_test = Est(test_data');
    case 'Fischer_2'
        load 'app_func_data/Fischer_2.mat'
        Est_test = predict(MdlLinear,test_data(:,[1 3]));
    case 'Fischer_3'
        load 'app_func_data/Fischer_3.mat'
        Est_test = predict(MdlLinear,test_data);
end

if strcmp(mode,'manual') == 0
    Est_test = strcmp(Est_test,'LEGH');
end

Est_test = double(reshape(Est_test,[n,m]));
Est_test(Label == 1) = 2;
map = [138/255,43/255,226/255;1,1,0;0,128/255,0];
C = labeloverlay(color_img(1:n*10,1:m*10,:),categorical(imresize(Est_test+1,10,'nearest'),[1,2,3]),...
    'Transparency',0.7,'Colormap',map);
imshow(C)
LEGH_ratio_e = Est_test;
LEGH_ratio_e(Est_test == 2) = 0;
LEGH_ratio_e = sum(LEGH_ratio_e(:))/(n*m-sum(Label(:)));
disp(['LEGH ratio: ', num2str(LEGH_ratio_e,4)])
