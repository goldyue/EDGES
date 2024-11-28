function [W1,W2,H1,H2,x1_predict,p] = EDGES_final(X,Y,Z,L1,L2,lambda1,lambda2,gamma1,gamma2,seita1,seita2,tol,K,iterMax)     
    %[X,Y,Z,W1,W2,H1,H2] = rowdata(X1,X2,X3);
    profile('on');
    t1 = clock;	
    %X1 = importdata('D:\jyzhao\R-4.2.2\merfish_moffit\stdata_X18.csv');
    %X = X1.data;
    %X1 = readmatrix('D:\jyzhao\R-4.2.2\merfish_moffit\stdata_X18.csv','OutputType', 'string');
     %X = str2double(sortrows(X1));
      %X(:,1) = []; 
    X1 = readmatrix('C:\Users\13089\OneDrive\Desktop\data16\stdatax11.csv','OutputType', 'string');
    X2 = readmatrix('C:\Users\13089\OneDrive\Desktop\data16\RNA_X21.csv','OutputType','string');
    X3 = readmatrix('C:\Users\13089\OneDrive\Desktop\data16\RNA_X31top2000.csv','OutputType','string');

    X = str2double(sortrows(X1));
    Y = str2double(sortrows(X2));
    Z = str2double(X3);
    X(:,1) = []; 
    Y(:,1) = [];
    Z(:,1) = [];
    
    %X(:,all(X==0,1)) = [];
    
    row_W1 = size(X,1);
    row_W2 = size(Z,1);
    col_H1 = size(X,2);
    col_H2 = size(Y,2);
    
    s= rng;
    
    %t1 = clock;
    
    rng(s);
    W1 = rand(row_W1,K); 
    W2 = rand(row_W2,K);
    H1 = rand(K,col_H1); 
    H2 = rand(K,col_H2); 
    
    e1k = eye(1,K); 
    e1kt = e1k';
   
    [delta_init1,~,~,~] = intermedia_v(W1,W2,H1,H2,X,Y,Z,L1_final,lambda1,lambda3,gamma1,gamma2,K);

    delta2 = delta_init1;

    stop = [];
    %delt1all = [];
    %delt2all = [];
    X1_lossall = [];
    X2_lossall = [];
    X3_lossall = [];
    H1_lossall = [];
    H2_lossall = [];
    W1_loss = [];
    W2_loss = [];
    H_lossall = [];
   total_lossall = [];
     
    
        for iter = 1:500
            disp(iter);

            %[W1tX1,W1tW1H1] = H1_update(W1,H1,X);
            %h1 = H1.*((W1tX1)./(W1tW1H1));
            %H1 = h1;
 
            %[W1tW1H2,W1tX2,W2tX3,W2tW2H2] = H2_update(W1,W2,H2,Z,Y);
            %h2 = H2.*((W1tX2 + gamma2*W2tX3)./(W1tW1H2 + gamma2*W2tW2H2));
            %H2 = h2;
            
            [X1H1t,X2H2t,W1H1H1t,W1H2H2t] = W1_update(W1,H1,H2,X,Y);
            w1 = W1.*((gamma1*X1H1t + X2H2t)./(gamma1*W1H1H1t + W1H2H2t));
            W1 = w1;
            
            [W1tX1,W1tW1H1,H1L1,ekkH1] = H1_update(W1,H1,X,L1_final,K);
            h1 = H1.*((gamma1*W1tX1- lambda1*H1L1)./(gamma1*W1tW1H1+ lambda3*ekkH1));
            %分子- lambda1*H1L1
            %分母+ lambda3*ekkH1
            H1 = h1;
            
            [X3H2t,W2H2H2t] = W2_update(Z,W2,H2);
            w2 = W2.*((X3H2t)./(W2H2H2t));
            W2 = w2;            
             
            [W1tW1H2,W1tX2,W2tX3,W2tW2H2,ekkH2] = H2_update(W1,W2,H2,Z,Y,K);
            h2 = H2.*((W1tX2 + gamma2*W2tX3)./(W1tW1H2+ lambda3*ekkH2));
            %- lambda2*H2L2
            %分母+ lambda3*ekkH2
            
            H2 = h2;  
            
            %compute stop value
            
            %H1t = H1';
            
            H1H1t = H1*H1';
            H2H2t = H2*H2';
            
            
            X1_loss = (norm((X - W1*H1),'fro'))^2;
            X2_loss = (norm((Y - W1*H2),'fro'))^2;
            X3_loss = (norm((Z - W2*H2),'fro'))^2;
            H1_loss = trace(H1*L1_final*H1');
            %H2_loss = trace(H2*L2_final*H2');
            %W1_loss = trace(W1W1t);
            %W2_loss = trace(W2W2t);
            H_total_loss = e1k*H1H1t*e1kt + e1k*H2H2t*e1kt;
            
            total_loss = gamma1*X1_loss + X2_loss + gamma2*X3_loss+ lambda1*H1_loss + lambda3*H_total_loss; 
            % + lambda1*H1_loss 
            %+ lambda2*H2_loss
            %+ lambda3*H_total_loss
            
            stop_control1 = abs((delta2 - total_loss)/(delta_init1 - total_loss)); 
            stop_value = abs(stop_control1);

            stop = [stop;stop_value];
            %delt1all = [delt1all;total_loss];
            %delt2all = [delt2all;delta2];
            X1_lossall = [X1_lossall;X1_loss];
            X2_lossall = [X2_lossall;X2_loss];
            X3_lossall = [X3_lossall;X3_loss];
            H1_lossall = [H1_lossall;H1_loss];
            %H2_lossall = [H2_lossall;H2_loss];
            H_lossall = [H_lossall;H_total_loss];
            total_lossall = [total_lossall;total_loss];
            
                    if stop_value < tol
                        break;
                    end
                        delta2 = total_loss;
        end
    %p = W1*H1;
    p = W2*H1;
    p = p(2001:2051,:);
    csvwrite('C:\Users\13089\OneDrive\Desktop\data16\edges_pre1fold.csv',p)
    %csvwrite('C:\Users\13089\OneDrive\Desktop\data16\edges_pre未检测不分折.csv',pp)
    %csvwrite('D:\jyzhao\R-4.2.2\data9\matlab2020\缺项验证\参-110-1-4\edges_未检测preL3.csv',pp)
    %csvwrite('C:\Users\13089\OneDrive\Desktop\data16\edges_W1_1fold.csv',W1)
    %csvwrite('C:\Users\13089\OneDrive\Desktop\data16\edges_W2_1fold.csv',W2)
    %csvwrite('C:\Users\13089\OneDrive\Desktop\data16\edges_H1_1fold.csv',H1)
    %csvwrite('C:\Users\13089\OneDrive\Desktop\data16\edges_H2_1fold.csv',H2)
    %csvwrite('D:\jyzhao\R-4.2.2\stplus中data2osmzei\matlab2020\缺项验证\edges_pre_未检测L3.csv',pp)
    t2 = clock;
    profile('off');
    memory_report = profile('info');
    disp(memory_report);
    t = t + etime(t2,t1);
    
    plot(stop','r*-');
    
    subplot(3,3,1)
    plot(X1_lossall','r*-')
    title('X1')
    
    subplot(3,3,2)
    plot(X2_lossall','r*-')
    title('X2')
    
    subplot(3,3,3)
    plot(X3_lossall','r*-')
    title('X3')
    
    subplot(3,3,4)
    plot(H1_lossall','r*-')
    title('L1')
    
    %subplot(4,4,5)
    %plot(H2_lossall','r*-')
    %title('L2')
       
    subplot(3,3,5)
    plot(H_lossall','r*-')
    title('H_total')
    
    subplot(3,3,6)
    plot(total_lossall','r*-')
    title('total')
    
    p = W2*H1;
    x1_predict = W1*H1;

        %filetitle_x1_predict = ['Z:\','stdata_yuce_x1''.xlsx'];
        %xlswrite(filetitle_x1_predict,x1_predict);
        %filetitle_p = ['Z:\','stdata_kong''.xlsx'];
        %xlswrite(filetitle_p,p);
        
        %data1
        pt1 = p(348:388,:);
        pt2 = p(348:381,:);
        pt3 = p(348:381,:);
        pt4 = p(348:381,:);
        pt5 = p(348:381,:);
        pt6 = p(348:381,:);
        pt7 = p(348:381,:);
        pt8 = p(348:381,:);
        pt9 = p(348:381,:);
        pt10 = p(348:381,:);
        %data1 2000
        pt1 = p(2001:2041,:);
        pt2 = p(2001:2034,:);
        pt3 = p(2001:2034,:);
        pt4 = p(2001:2034,:);
        pt5 = p(2001:2034,:);
        pt6 = p(2001:2034,:);
        pt7 = p(2001:2034,:);
        pt8 = p(2001:2034,:);
        pt9 = p(2001:2034,:);
        pt10 = p(2001:2034,:);
        %data2
        p1 = p(15043:15053,:);
        p2 = p(15043:15053,:);
        p3 = p(15043:15053,:);
        %data2 2000genes
        p1 = p(2001:2011,:);
        p2 = p(2001:2011,:);
        p3 = p(2001:2011,:);
         %data9 2000
        pt1 = p(2001:2031,:);
        pt2 = p(2001:2024,:);
        pt3 = p(2001:2024,:);
        pt4 = p(2001:2024,:);
        pt5 = p(2001:2024,:);
        pt6 = p(2001:2024,:);
        pt7 = p(2001:2024,:);
        pt8 = p(2001:2024,:);
        pt9 = p(2001:2024,:);
        pt10 = p(2001:2024,:);
         %data14
        ppp1 = p(10070:10097,:);
        ppp2 = p(10070:10096,:);
        ppp3 = p(10070:10096,:);
        %data14_2000
        ppp1 = p(2001:2028,:);
        ppp2 = p(2001:2027,:);
        ppp3 = p(2001:2027,:);
        %data17
        pp1 = p(30956:30969,:);
        pp2 = p(30956:30969,:);
        pp3 = p(30956:30969,:);
         %data172000
        pp1 = p(2001:2014,:);
        pp2 = p(2001:2014,:);
        pp3 = p(2001:2014,:);
        %oam_allen
        pppp1 = p(30965:30975,:);
        pppp2 = p(30965:30975,:);
        pppp3 = p(30965:30975,:);
        %oam_allen2000
        pppp1 = p(2001:2011,:);
        pppp2 = p(2001:2011,:);
        pppp3 = p(2001:2011,:);
        %starmap_allen 2000
        pt1 = p(2001:2099,:);
        pt2 = p(2001:2099,:);
        pt3 = p(2001:2099,:);
        pt4 = p(2001:2099,:);
        pt5 = p(2001:2099,:);
        pt6 = p(2001:2101,:);
        pt7 = p(2001:2101,:);
        pt8 = p(2001:2101,:);
        pt9 = p(2001:2101,:);
        pt10 = p(2001:2101,:);
        
        pt_final = [pt1;pt2];
        pt_final2 = [pt_final;pt3];
        pt_final3 = [pt_final2;pt4];
        pt_final4 = [pt_final3;pt5];
        pt_final5 = [pt_final4;pt6];
        pt_final6 = [pt_final5;pt7];
        pt_final7 = [pt_final6;pt8];
        pt_final8 = [pt_final7;pt9];
        pt_final9 = [pt_final8;pt10];
        
        csvwrite('D:\jyzhao\R-4.2.2\merfish_moffit\edges_pre_1fold2000.csv',pt1)
        csvwrite('D:\jyzhao\R-4.2.2\merfish_moffit\edges_pre_2fold2000.csv',pt2)
        csvwrite('D:\jyzhao\R-4.2.2\merfish_moffit\edges_pre_3fold2000.csv',pt3)
        csvwrite('D:\jyzhao\R-4.2.2\merfish_moffit\edges_pre_4fold2000.csv',pt4)
        csvwrite('D:\jyzhao\R-4.2.2\merfish_moffit\edges_pre_5fold2000.csv',pt5)
        csvwrite('D:\jyzhao\R-4.2.2\merfish_moffit\edges_pre_6fold2000.csv',pt6)
        csvwrite('D:\jyzhao\R-4.2.2\merfish_moffit\edges_pre_7fold2000.csv',pt7)
        csvwrite('D:\jyzhao\R-4.2.2\merfish_moffit\edges_pre_8fold2000.csv',pt8)
        csvwrite('D:\jyzhao\R-4.2.2\merfish_moffit\edges_pre_9fold2000.csv',pt9)
        csvwrite('D:\jyzhao\R-4.2.2\merfish_moffit\edges_pre_10fold2000.csv',pt10)
        
        boxplot(diag(corr(X1_row',X1_pre')))
        boxplot(diag(corr(X',x1_predict')))
        
        %所有结果图合并
        boxplot(diag(corr(pt_final1,X)),diag(corr(pt_final1,X)),'Labels',{'mu = china','mu = usa'},'Widths',0.5,'DataLim',[-1,1])
        title('Benchmarking')
        
        %所有结果图合并
            % Create plots
            t = tiledlayout(7,1);
            ax1 = nexttile;
            plot(ax1,x1,y1)
            ax2 = nexttile;
            stem(ax2,x2,y2)
            % Link the axes
            linkaxes([ax1,ax2],'x');
            % Add shared title and axis labels
            title(t,'Benchmarking')
            xlabel(t,'method')
            ylabel(t,'PCC')

            % Move plots closer together
            xticklabels(ax1,{})
            t.TileSpacing = 'compact';
end