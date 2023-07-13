%% M1-drived bbNAP
clc;clear
load('D:\Desktop\m\data\ESAoccci_data\V5.0\monthly\bin1\chl.mat');
load('D:\Desktop\m\data\ESAoccci_data\V5.0\monthly\bin1\bbp.mat');
chl(chl<=0) = NaN;bbp(bbp<=0) = NaN;
nap_B18 = [];k = [];r2_B18 = [];
for i = 1:360
    for j = 1:181
        X = [];Y = [];aa1 = [];aa2 = [];num = 0;
        x = chl(i,j,:);y = bbp(i,j,:);
        x = x(:);y = y(:);x1 = x;y1 = y;
        n = 3;
        x(x1 >= (nanmean(x1)+n*nanstd(x1))) = NaN;x(x1 <= (nanmean(x1)-n*nanstd(x1))) = NaN;
        y(y1 >= (nanmean(y1)+n*nanstd(y1))) = NaN;y(y1 <= (nanmean(y1)-n*nanstd(y1))) = NaN;%%用三个标准偏差去除异常值
        x(find(isnan(y))) = NaN;y(find(isnan(x))) = NaN;
        for ii = 1:276
            if isnan(x(ii)) == 0
                num = num+1;
                aa1(num) = x(ii);aa2(num) = y(ii);
            end
        end
        if num >= 20
            aa = aa1';bb = aa2';
            AA = [ones(size(bb)),aa];
            [b,bint,r,rint,stats] = regress(bb,AA);
            nap_B18(i,j) = b(1);
            k(i,j)   = b(2);
            r2_B18(i,j)   = stats(1);
        else if num < 20
                nap_B18(i,j)     = NaN;
                k(i,j)           = NaN;
                r2_B18(i,j)      = NaN;
            end
        end
    end
end
nap_B18(k<=0) = NaN;
r2_B18(k<=0) = NaN;

save('B18.mat','nap_B18','r2_B18');

