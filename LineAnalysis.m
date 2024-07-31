%
% This software is released under the GPL v3. It is provided AS-IS and no
% warranty is given.
%
% Author: Ye Zhiwei, 2024
%

function [Results] =  LineAnalysis(xs,ys, xw, yw,Magsize, gauss_num)
%LINEANALYSIS ��SRͼ��WFͼ��line�����ݽ��и�˹���
% gauss_num, number of gauss fcn to fit.
%   Magsize �Ŵ���
% test
% addr_SR = 'E:\Super_resolution_data\YE Zhiwei\2017.04.08\5.Lyso-R\500nM\9\2\SR.csv';
% addr_F = 'E:\Super_resolution_data\YE Zhiwei\2017.04.08\5.Lyso-R\500nM\9\2\F.csv';
% Magsize = 10;
% % fid1 = fopen(addr_SR);
% % fid2 = fopen(addr_F);
% % ori_d1 = textscan(fid1,'%f,%f','HeaderLines',1);
% % ori_d2 = textscan(fid2,'%f,%f','HeaderLines',1);
% % fclose(fid1);
% % fclose(fid2);
% % xw=ori_d2{1};
% % yw=ori_d2{2}; 
% % % filter out noise
% % % yw = yw-min(yw);
% % ys=ori_d1{2};
% % xs=ori_d1{1};

%Normalization
nys = ys./max(ys);
nyw = yw./max(yw);

[ fit_type,fit_option ] = GaussianFitInitializer( max(xw), gauss_num );
[ fit_typeS,fit_optionS ] = GaussianFitInitializer( max(xs), gauss_num );
[ fit_optionS ] = StartPointInitializer( xs, nys, gauss_num, fit_optionS );
[fw,goffw] = fit(xw, nyw,fit_type, fit_option);
[fs,goffs] = fit(xs, nys, fit_typeS, fit_optionS);
disp(goffw.adjrsquare); disp(goffs.adjrsquare);
width_WF = zeros(gauss_num, 1); center_WF = zeros(gauss_num, 1);
width_SR = zeros(gauss_num, 1); center_SR = zeros(gauss_num, 1);
for i = 1 : gauss_num
   eval(['width_WF(i) = fw.c' num2str(i) ' * 2 * sqrt(log(2));']);
   eval(['center_WF(i) = fw.b' num2str(i) ';']);
   eval(['width_SR(i) = fs.c' num2str(i) ' * 2 * sqrt(log(2));']);
   eval(['center_SR(i) = fs.b' num2str(i) ';']);
end
Results = struct();
Results.width_WF = width_WF; Results.center_WF = center_WF;
Results.width_SR = width_SR; Results.center_SR = center_SR;
%%
xf=[0:10000];
xf=xf';
xf=xf.*(max(xw)/10000);
ywf=fw(xf);
ysf=fs(xf);

[ xw_bin, nyw_bin ] = Binning_WideFieldData(xw, nyw, Magsize);

% %% ����������Ϻ����ֱ�Ϊfw��fs���õ�����igor��ͼ������ �²�һ����˹���
% gau_n_w = 1;
% gau_n_s = 1;
% %s = fitoptions('gauss1');
% %f = fittype('gauss1');
% mygaussfit = @(a,b,c,d,x) a+(b-a).*exp(-(x-c).*(x-c)./(2.*d.*d));
% myfittype = fittype(mygaussfit); myfitoption = fitoptions(myfittype);
% myfitoption.Lower = [0 -Inf xs(1) xs(1)]; myfitoption.StartPoint = [0.05 1 median(xs) max(xs)/5];
% myfitoption.Upper = [1 Inf max(xs) median(xs)];
% [fw,goffw] = fit(xw,nyw,myfittype,myfitoption);
% [fs,goffs] = fit(xs,nys,myfittype,myfitoption);
% 
% % if goffw.adjrsquare > 0.8
% %    disp('Wide-field data fitting is ok!');clear gof3;
% % else
% %    s = fitoptions('gauss2');
% %    f = fittype('gauss2');
% %    s.Lower = [0 -Inf 0 0 -Inf 0];
% %    s.Upper = [1 Inf Inf 1 Inf Inf];
% %    [fw,goffw] = fit(xw,nyw,f,s);
% %    gau_n_w = 2;
% %    if (goffw.adjrsquare > 0.9)
% %        disp('Wide-field data fitting is ok!');
% %    else
% %        disp('Wide-field data fitting is wrong!');
% %    end
% % end
% 
% if goffs.adjrsquare > 0.9
%     disp('Super-resolution data fitting is ok!');
% else
%     s = fitoptions('gauss2');
%     f = fittype('gauss2');
%     s.Lower = [0 -Inf 0 0 -Inf 0];
%     s.Upper = [1 Inf Inf 1 Inf Inf];
%     [fs,goffs] = fit(xs,nys,f,s);
%     gau_n_s = 2;
%     s = fitoptions('gauss2');
%     f = fittype('gauss2');
%     s.Lower = [0 -Inf 0 0 -Inf 0];
%     s.Upper = [1 Inf Inf 1 Inf Inf];
%     [fw,goffw] = fit(xw,nyw,f,s);
%     gau_n_w = 2;
%     if goffs.adjrsquare > 0.9
%         disp('Super-resolution data fitting is ok!');
%     else
%         disp('Super-resolution data fitting is wrong!');
%     end
%             disp(['The R^2 is ', num2str(goffs.adjrsquare, 4)]);
% end
%% Save plot Width;
% Results = struct();
% if gau_n_w == 1
%     width_WF1 = fw.d * 2.35482;
%     center_WF1 = fw.c;
%     Results.width_WF = width_WF1;
%     Results.center_WF = center_WF1;
% else if gau_n_w ==2
%     width_WF1 = fw.c1 * 2 * sqrt(log(2));
%     width_WF2 = fw.c2 * 2 * sqrt(log(2));
%     center_WF1 = fw.b1;
%     center_WF2 = fw.b2;
%     Results.width_WF = [width_WF1 width_WF2];
%     Results.center_WF = [center_WF1, center_WF2];
%     end
% end
% 
% if gau_n_s == 1
%     width_SR1 = fs.d  * 2.35482;
%     Results.width_SR = width_SR1;
%     Results.center_SR = fs.c;
%     else if gau_n_s == 2
%             width_SR1 = fs.c1 * 2 * sqrt(log(2));
%             width_SR2 = fs.c2 * 2 * sqrt(log(2));
%             center_SR1 = fs.b1; center_SR2 = fs.b2;
%             Results.width_SR = [width_SR1 width_SR2];
%             Results.center_SR = [center_SR1, center_SR2];
%         end
% end

% %%
% %��������binһ��
% samecount = find(diff([nyw;1.05])) - find(diff([-0.05;nyw]))+1;
% xnumels = numel(xw);
% nyw_bin = zeros(ceil(xnumels/ Magsize),1);
% xw_bin = zeros(ceil(xnumels/ Magsize),1);
% diffxw = xw(2) - xw(1);
% xw_bin(1) = ( sum(xw(1:samecount(1))) + sum((1:floor( (Magsize - samecount(1)) ))*(-1)*diffxw) ) / Magsize;
% nyw_bin(1) = sum(nyw(1:samecount(1))) / samecount(1);
% k = samecount(1);    i = 2;
% while 1
%     j = k + Magsize; x1 = round(k+1); x2 = round(j);
%     if j <= xnumels
%         nyw_bin(i) = sum(nyw(x1:x2))/Magsize;
%         xw_bin(i) = sum(xw(x1:x2))/Magsize;
%         k = j; i = i +1;
%     else
%         nyw_bin(i) = sum(nyw(x1:end))/(xnumels - k);
%         xw_bin(i) = ( sum(xw(x1:end)) + sum((1:floor( (Magsize - xnumels+k) ))*diffxw + xw(end)) )/Magsize;
%         break;
%     end
% end
% xs = xs - min(xw_bin);  xf = xf - min(xw_bin); xw_bin = xw_bin - min(xw_bin);
% 
% %ת��umΪnm fix me !
% xs=xs.*1000;xw_bin=xw_bin.*1000;xf=xf.*1000;

Results.xs = xs; Results.ys = nys; Results.xw = xw_bin; Results.yw = nyw_bin;
Results.xf = xf; Results.ysf = ysf; Results.ywf = ywf;

h=figure;
ax = axes;
hold on;
linewidth_plot=0.5;
%mgraphics(h, [2.1 1.8]);% Fix me!
%mgraphics(h, [3.5 3]);% Fix me!

p1 = plot (ax, xs,nys,'.','color',Colorslib('red')/256,'MarkerSize',5);
p2 = plot (ax, xf,ysf, '-', 'color',Colorslib('red')/256, ...
     'LineWidth',linewidth_plot);
p3 = plot (ax, xw_bin,nyw_bin,'x','color',Colorslib('grey')/256/1.8,'MarkerSize',3);
p4 = plot (ax, xf,ywf, '-', 'color',Colorslib('grey')/256/1.8, ...
     'LineWidth',linewidth_plot);
 
 
[L1 L2 L3 L4]=legend([''],['SR'],[''],['WF']);
set(L1,'Fontsize',5);
set(L1,'box','off');

    for i1 = 1:4
    L2(i1).FontSize = 5;
    end  

%����һ��,����������legend��Ϊ���ɼ�
set (L2(5),'Visible','off');
set (L2(6),'Visible','off');
set (L2(1),'Visible','off');
set (L2(9),'Visible','off');
set (L2(10),'Visible','off');
set (L2(3),'Visible','off');

Pos1 = get (L2(1),'Position');
Pos2 = get (L2(2),'Position');
Pos3 = get (L2(3),'Position');
Pos4 = get (L2(4),'Position');

Pos4(2) = Pos2 (2);
Pos2(2) = Pos1 (2);

Posx1=get(L2(5),'XData');
Posy1=get(L2(5),'YData');
Posx2=get(L2(7),'XData');
Posy2=get(L2(7),'YData');
Posx3=get(L2(9),'XData');
Posy3=get(L2(9),'YData');
Posx4=get(L2(11),'XData');
Posy4=get(L2(11),'YData');

set(L2(7),'YData',Posy1);
set(L2(11),'YData',Posy2);
set(L2(2),'Position',Pos2);
set(L2(4),'Position',Pos4);

%��С�ڶ���ͼ��line�Ĵ�С
Posx2=get(L2(7),'XData');
Posy2=get(L2(7),'YData');
Posx2(2)=(Posx2(2)-Posx2(1))/4+Posx2(1);
set(L2(7),'XData',Posx2);
%�ڶ���ͼ��textλ�ø��Ÿı�
Pos2=get(L2(2),'Position');
Pos2(1)=Pos2(1)-((Posx2(2)-Posx2(1))*4+Posx2(1))+Posx2(2);
set(L2(2),'Position',Pos2);


%��С���ĸ�ͼ��line�Ĵ�С
Posx4=get(L2(11),'XData');
Posy4=get(L2(11),'YData');
Posx4(2)=(Posx4(2)-Posx4(1))/4+Posx4(1);
set(L2(11),'XData',Posx2);
%���ĸ�ͼ��textλ�ø��Ÿı�
Pos4=get(L2(4),'Position');
Pos4(1)=Pos4(1)-((Posx4(2)-Posx4(1))*4+Posx4(1))+Posx4(2);
set(L2(4),'Position',Pos4);

%�ı�Legendλ��
PosL = get(L1, 'Position');
PosL(1)= 0.70;
PosL(2) = 0.65;
set(L1,'Position',PosL);
hold off;
Results.figure = h;
end

