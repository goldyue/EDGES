X1 = csvread('stdata2_allcell_order.csv',1,1);
X2 = csvread('Z:\jyzhao\R-4.2.2\stplus中data2\10fold_pre_order.csv',1,1);
X2_pre = csvread('Z:\jyzhao\R-4.2.2\stplus中data2\stdata2_allcell_order.csv',1,1);
X3 = csvread('RNAdata1X31.csv',0,0);
%分折
X1 = readmatrix('D:\jyzhao\R-4.2.2\data1\分折\data1_X12.csv','OutputType','string');
X2 = readmatrix('D:\jyzhao\R-4.2.2\data1\分折\RNAdata1_X22.csv','OutputType','string');
X3 = readmatrix('D:\jyzhao\R-4.2.2\data1\分折\RNAdata1_X32_top2000.csv','OutputType','string');

%MATLAB中数据专用
%X1 = readmatrix('stdata用转置.csv', 'OutputType', 'string');
X1 = readmatrix('stdata_osmallen_X11.csv', 'OutputType', 'string');
X2 = readmatrix('RNAdata_osmallen_X21.csv', 'OutputType', 'string');
X3 = readmatrix('RNAdata14_X33.csv','OutputType', 'string');

%文件夹数据
X1 = readmatrix('D:\jyzhao\R-4.2.2\merfish_moffit\MERFISH\MERFISH_X1.csv','OutputType','string');
X2 = readmatrix('D:\jyzhao\R-4.2.2\merfish_moffit\Moffit_RNA\RNAdata_X2.csv','OutputType','string');
X3 = readmatrix('D:\jyzhao\R-4.2.2\merfish_moffit\Moffit_RNA\datamm_X3_top2000.csv','OutputType','string');
%X1 = readmatrix('D:\jyzhao\R-4.2.2\stplus中data2osmzei\stdata_data2X11.csv','OutputType','string');
%X2 = readmatrix('D:\jyzhao\R-4.2.2\stplus中data2osmzei\RNAdata_data2X21.csv','OutputType','string');
%X3 = readmatrix('D:\jyzhao\R-4.2.2\data17\RNAdata17_X33_top2000.csv','OutputType','string');
%分折
X1 = readmatrix('D:\jyzhao\R-4.2.2\data4\stdata\stdata4_X1.csv','OutputType','string');
X2 = readmatrix('D:\jyzhao\R-4.2.2\data4\RNAdata\RNAdata4_X2.csv','OutputType','string');
X3 = readmatrix('D:\jyzhao\R-4.2.2\data4\RNAdata\RNAdata4_X3.csv','OutputType','string');
%预测
X1_pre = csvread('D:\jyzhao\R-4.2.2\stplus中data2osmzei\matlab2020\缺项验证\6427细胞\edges_pre.csv',1,1);
X1_row = csvread('D:\jyzhao\R-4.2.2\stplus中data2osmzei\stdata2_6427_order_X1.csv',1,1);

%X1_row = csvread('D:\jyzhao\R-4.2.2\data9\edges_row2000_order.csv',1,1);
%X1_row = csvread('D:\jyzhao\R-4.2.2\stplus中data2osmzei\stdata2_allcell_order_X1.csv',1,1);
%X1_row = csvread('D:\jyzhao\R-4.2.2\data14\X1_3折_order.csv',1,1);
%X1_row = csvread('D:\jyzhao\R-4.2.2\data17\X1_3折_order.csv',1,1);
%X1_row = csvread('D:\jyzhao\R-4.2.2\data1\分折\X1_order_row.csv',1,1);
%X1_row = csvread('D:\jyzhao\R-4.2.2\merfish_moffit\tangram_roworder.csv',1,1);
%X1_row = csvread('D:\jyzhao\R-4.2.2\data4\stdata4_X1order.csv',1,1);
%X1_row = csvread('D:\jyzhao\R-4.2.2\data_starmap_allen\tangram_row_order2000.csv',1,1);
%X1_row = csvread('D:\jyzhao\R-4.2.2\osm_allen\X1_roworder.csv',1,1);


%L1 =csvread('D:\jyzhao\R-4.2.2\osm_allen\osmFISH\knn_osmfish.csv',1,1);%D-A
%L1_duliang = diag(diag(L1));
%L1 = L1_duliang - L1;

%L1 = csvread('knn_data14第二版.csv',1,1);
%L1_duliang = diag(diag(L1));
%L1 = L1_duliang - L1;

%L1 = csvread('D:\jyzhao\R-4.2.2\data1\L1.csv',1,1);
%L1_duliang = csvread('D:\jyzhao\R-4.2.2\data1\L1_duliang.csv',1,1);
%L1(:,31579) = [];
%L1 = L1_duliang - L1;邻接矩阵

%L1 = csvread('D:\jyzhao\R-4.2.2\data_starmap_allen\Starmap\knn_starmap.csv',1,1);%knn邻接矩阵
%L1_duliang = diag(diag(L1));
%L1 = L1_duliang - L1;

%L1 = csvread('D:\jyzhao\R-4.2.2\data9\knn_data9_less.csv',1,1);%应该是D-A
%L1_duliang = diag(diag(L1));
%L1 = L1_duliang - L1;

%L1 = csvread('D:\jyzhao\R-4.2.2\merfish_moffit\MERFISH\knn_merfish.csv',1,1);%应该是knn邻接矩阵
%L1_duliang = diag(diag(L1));
%L1 = L1_duliang - L1;

%L1 = csvread('D:\jyzhao\R-4.2.2\data4\stdata\knn_data4.csv',1,1);%应该是knn邻接矩阵
%L1_duliang = diag(diag(L1));
%L1 = L1_duliang - L1;

%L1 = csvread('D:\jyzhao\R-4.2.2\stplus中data2osmzei\osm_knn.csv',1,1);%应该是knn邻接矩阵
%rowSums = sum(L1, 2);
%L1_duliang = diag(rowSums);

%L1 = csvread('D:\jyzhao\R-4.2.2\stplus中data2osmzei\matlab2020\新标准化\缺项验证\knn_data2_6427.csv',1,1);%应该是knn邻接矩阵
%L1_duliang = diag(diag(L1));
%L1 = L1_duliang - L1;

%L1 =csvread('D:\jyzhao\R-4.2.2\data17\knn_data17.csv',1,1);%D-A
%L1_duliang = diag(diag(L1));
%L1 = L1_duliang - L1;

L1 =csvread('C:\Users\13089\OneDrive\Desktop\data16\knn_data16.csv',1,1);%D-A
L1_duliang = diag(diag(L1));
L1 = L1_duliang - L1;%邻接

L1_new = (L1_duliang^(-0.5))*L1*(L1_duliang^(-0.5));
% 获取 L1 矩阵的维度信息
[m, n] = size(L1_duliang);

% 创建和 L1 矩阵同维度的单位矩阵
I = eye(m, n);
L1_final = I - L1_new;


L2 = csvread('D:\jyzhao\R-4.2.2\data9\Allen_VISp\L2network01.csv',1,1);
%L2 = csvread('D:\jyzhao\R-4.2.2\data17\L2network.csv',1,1);

%L2 = csvread('L2_network_data14.csv',1,1);
%L2_linjie(L2~=0) = 1;
L2_diag = diag(diag(L2));
L2_linjie = L2-L2_diag;

L2_duliang = sum(L2_linjie,2);
L2_duliang = diag(L2_duliang);
L2_new = (L2_duliang^(-0.5))*L2_linjie*(L2_duliang^(-0.5));
[s, q] = size(L2_duliang);

% 创建和 L1 矩阵同维度的单位矩阵
E = eye(s, q);
L2_final = E - L2_new;

lambda1 = 10^(-5);
%lambda1 = 10^(-5);
lambda2 = 10^(1);
%lambda3 = 10^(1);
lambda3 = 10^(1);

gamma1 = 10^(-1);
%gamma1 = 10^(-1);
%gamma2 = 10^(-4);
gamma2 = 10^(-4);

seita1 = 10^(0);
seita2 = 10^(0);

tol = 10^(-7);
k = 10;
K = 20;
t = 10;
iterMax = 1000;

function[X1_div,X2_div,X3_div] = data1_divide()
% 10 pieces of datarandomly generate
s = rng;
rng(s);
X1_sort = sortrows(X1);
X2_sort = sortrows(X2);
rowrank_X1 = randperm(size(X1_sort, 1));
X1_luan = X1_sort(rowrank_X1,:);
X2_luan = X2_sort(rowrank_X1,:);

m = size(X1_luan,1);
fold_num = floor(m/k);
yushu = mod(m,k);

fold = ones(k,1)*nan;

for i=1:k
    if i==1
        fold(i) = fold_num+yushu;
    else
        fold(i) = fold_num;
    end
end

X1_divide = mat2cell(X1_luan,fold);
X2_divide = mat2cell(X2_luan,fold);

%store 10 copies of data
X1_div = cell(k,1);
X2_div = cell(k,1);
X3_div = cell(k,1);
stX1_predict_row = cell(k,1);
X1_predict_rownames = cell(k,1);
p_all = cell(k,1);
pcc = cell(k,1);


    for i = 1:k
           sum_X1 = [];
           sum_X2 = [];
           sum_X3 = [];
           sum_stX1_pre = [];
 
           for j = 1:k 
               if j~=i
                   sum_X1 = [sum_X1;X1_divide{j}];
                   sum_X2 = [sum_X2;X2_divide{j}];
               else
                   sum_X3 = [X3;X2_divide{j}];
                   sum_stX1_pre = [sum_stX1_pre;X1_divide{j}];
               end
           end
  
            X1_div{i} = sum_X1;
            X2_div{i} = sum_X2;
            X3_div{i} = sum_X3;
            stX1_predict_row{i} = sum_stX1_pre;
            X1_predict_rownames{i} = stX1_predict_row{i}(:,1);

             filetitle_X = ['stdata_data4X1' num2str(i) '.xlsx'];
             filetitle_Y = ['RNAdata_data4X2' num2str(i) '.xlsx'];
             filetitle_Z = ['RNAdata_data4X3' num2str(i) '.xlsx'];
             filetitle_X_predict = ['stdata_data4pre' num2str(i) '.xlsx'];
             filetitle_X_predict_name = ['stdata_data4prename' num2str(i) '.xlsx'];
% %         
%         filetitle_X = ['stdata14X1' num2str(i) '.csv'];
%         filetitle_Y = ['RNAdata14X2' num2str(i) '.csv'];
%         filetitle_Z = ['RNAdata14X3' num2str(i) '.csv'];
%           filetitle_X_predict = ['stdata1pre' num2str(i) '.csv'];
%         filetitle_X_predict_name = ['stdata14prename' num2str(i) '.csv'];
        
        
                  dataTmp_X = X1_div{i,1};
                  dataTmp_Y = X2_div{i,1};
                  dataTmp_Z = X3_div{i,1};
                  dataTmp_X_predict = stX1_predict_row{i,1};
                  dataTmp_X_predict_name = X1_predict_rownames{i,1};
%         
%           csvwrite(filetitle_X,dataTmp_X);
%           csvwrite(filetitle_Y,dataTmp_Y);
%           csvwrite(filetitle_Z,dataTmp_Z);
%           csvwrite(filetitle_X_predict,dataTmp_X_predict);
%           csvwrite(filetitle_X_predict_name,dataTmp_X_predict_name);
%         
                  xlswrite(filetitle_X,dataTmp_X);
                  xlswrite(filetitle_Y,dataTmp_Y);
                  xlswrite(filetitle_Z,dataTmp_Z);
                  xlswrite(filetitle_X_predict,dataTmp_X_predict);
                  xlswrite(filetitle_X_predict_name,dataTmp_X_predict_name);

        %[W1,W2,H1,H2,x1_predict,p] = EDGES_final(sum_X1{i},sum_X2{i},sum_X3{i},L1,L2,lambda1,lambda2,gamma1,gamma2,seita1,seita2,tol,K,iterMax);

        %p_all{i} = p;
        %pcc{i} = x1_predict;

        %stX1_predict_row{i}(:,1) = [];
        %pcc{i} = corr((x1_predict(size(X3,1)+1,size(x1_predict)))',(str2double(stX1_predict_row{i}))');
        %filetitle_pcc = ['Z:\','pcc' '.xlsx'];
        %dataTmp_pcc = pcc{i,1};
        %xlswrite(filetitle_pcc,dataTmp_pcc);

    end
end


%X1 = readmatrix('stdata_norm_havename.csv', 'OutputType', 'string');
%X2 = readmatrix('data1_X2.csv', 'OutputType', 'string');
%X3 = readmatrix('data1_X3.csv','OutputType', 'string');
%L1 = csvread('L1.csv',1,1);
%L1_linjie = csvread('knn_fin.csv',1,1);
%L1_duliang = csvread('L1_duliang.csv',1,1);
%L1_new = (L1_duliang^(-0.5))*L1_linjie*(L1_duliang^(-0.5));
%L1_final = L1_duliang - L1_new;

%L2 = csvread('L2network_norm.csv',1,1);
%L2_linjie(L2~=0) = 1;
%L2_diag = diag(diag(L2));
%L2_linjie = L2-L2_diag;

%L2_duliang = sum(L2_linjie,2);
%L2_duliang = diag(L2_duliang);
%L2_new = (L2_duliang^(-0.5))*L2_linjie*(L2_duliang^(-0.5));
%L2_final = L2_duliang - L2_new;
