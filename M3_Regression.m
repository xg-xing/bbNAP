% %  M3-derived bbNAP
clc;clear
load('D:\Desktop\Cod\Data\Global data\Chl.mat');
load('D:\Desktop\Cod\Data\Global data\bbp.mat');
load('D:\Desktop\Cod\Data\Global data\SST.mat');
load('D:\Desktop\Cod\Data\Global data\PARg.mat');
nap_M3 = [];chlcmin = []; r2 = [];
for i = 1:360
    for j = 1:173
        num = 0;
        X = [];Y = [];
        aa1 = [];aa2 = []; aa3 = []; aa4 = [];
        x = bbp(i,j+8,:);
        y = 0.0155+0.00005*exp(0.215*sst(i,j+8,:));%%Behrenfeld 2005
        z = exp(-3*ig(i,j,:));
        w = chl(i,j+8,:);
        x = x(:);z = z(:);w = w(:);y = y(:);
        x1 = x;y1 = y;z1 = z;w1 = w;
        n = 3;
        x(x1 >= (nanmean(x1)+n*nanstd(x1))) = NaN;x(x1 <= (nanmean(x1)-n*nanstd(x1))) = NaN;
        y(y1 >= (nanmean(y1)+n*nanstd(y1))) = NaN;y(y1 <= (nanmean(y1)-n*nanstd(y1))) = NaN;%%用三个标准偏差去除异常值
        z(z1 >= (nanmean(z1)+n*nanstd(z1))) = NaN;z(z1 <= (nanmean(z1)-n*nanstd(z1))) = NaN;
        w(w1 >= (nanmean(w1)+n*nanstd(w1))) = NaN;w(w1 <= (nanmean(w1)-n*nanstd(w1))) = NaN;
        x(find(isnan(z))) = NaN;w(find(isnan(z))) = NaN;y(find(isnan(z))) = NaN;
        z(find(isnan(x))) = NaN;w(find(isnan(x))) = NaN;y(find(isnan(x))) = NaN;
        x(find(isnan(w))) = NaN;z(find(isnan(w))) = NaN;y(find(isnan(w))) = NaN;
        x(find(isnan(y))) = NaN;z(find(isnan(y))) = NaN;w(find(isnan(y))) = NaN;
        a0=[0.01,0.01];
        for ii = 1:276
            if isnan(x(ii)) == 0
                num      = num+1;
                aa1(num) = x(ii);aa2(num) = y(ii);aa3(num) = z(ii);aa4(num) = w(ii);
            end
        end
        if num >= 20
            X                = [aa1;aa2;aa3]';
            Y                = aa4';
            func             = @(a,X)13000*(X(:,1)-a(1)).*(a(2)+(X(:,2)-a(2)).*X(:,3));
            [a,resnorm]      = lsqcurvefit(func,a0,X,Y,[-0.1 0],[0.1 0.1]);
            r2(i,j)      = 1 - resnorm/sum((aa4 - mean(aa4)).^2);%%求R2指数
            nap_M3(i,j)     = a(1);
            chlcmin(i,j) = a(2);
        else if num < 20
                nap_M3(i,j)     = NaN;
                chlcmin(i,j) = NaN;
                r2(i,j)      = NaN;
            end
        end
    end
end

chlcmin(nap_M3<0) = NaN;
r2(nap_M3<0)      = NaN;
nap_M3(nap_M3<0)     = NaN;

% save('D:\Desktop\bbNAP_2022\cod-nap\m3_G87.mat','nap_m3G87','chlcmin_m3G87');

