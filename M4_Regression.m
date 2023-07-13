%% M4-derived bbNAP
clc;clear
load('D:\Desktop\Cod\Data\Global data\Chl.mat');
load('D:\Desktop\Cod\Data\Global data\bbp.mat');
load('D:\Desktop\Cod\Data\Global data\SST.mat');
load('D:\Desktop\Cod\Data\Global data\PARg.mat');
nap_M4 = [];chlcmin = [];chlcmax = []; r2 = [];
for i = 1:360
    for j = 1:173
        num = 0;X = [];Y = [];
        aa1 = []; aa3 = []; aa4 = [];
        x = bbp(i,j+8,:);z = exp(-3*ig(i,j,:)); w = chl(i,j+8,:);
        x(find(isnan(z))) = NaN;w(find(isnan(z))) = NaN;
        z(find(isnan(x))) = NaN;w(find(isnan(x))) = NaN;
        x(find(isnan(w))) = NaN;z(find(isnan(w))) = NaN;
        x = x(:);z = z(:);w = w(:);
        a0=[0.01,0.01,0.01];
        for ii = 1:276
            if isnan(x(ii)) == 0
                num      = num+1;
                aa1(num) = x(ii);aa3(num) = z(ii);aa4(num) = w(ii);
            end
        end
        if num >= 20
            X                = [aa1;aa3]';
            Y                = aa4';
            func             = @(a,X)13000*(X(:,1)-a(1)).*(a(2)+a(3)*X(:,2));
            [a,resnorm]      = lsqcurvefit(func,a0,X,Y,[-0.1 0 0],[0.1 0.1 0.1]);
            sstot            = sum((aa4 - nanmean(aa4)).^2);
            r2(i,j)      = 1 - resnorm/sstot;%%求R2指数
            nap_M4(i,j)      = a(1);
            chlcmin_m2(i,j)  = a(2);
            chlcmax(i,j)     = a(2)+a(3);
        else if num < 20
                nap_M4(i,j)     = NaN;
                chlcmin_m2(i,j) = NaN;
                chlcmax(i,j)    = NaN;
                r2(i,j)         = NaN;
            end
        end
    end
end

chlcmin_m2(nap_M4<0) = NaN;
r2(nap_M4<0)         = NaN;
nap_M4(nap_M4<0)     = NaN;
% save('m2.mat','nap_m2','r2_m2','chlcmin_m2');

