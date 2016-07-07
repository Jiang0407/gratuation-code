function debug()
    close all
    clear all
    clc  

    load Copy_of_BBXX-0707-4.mat
    pathL='F:/桌面杂/0707/imageL/imgL';
    pathR='F:/桌面杂/0707/imageR/imgR';
    [m,n,~]=size(warp_img);
%% 确定重叠区域在原始图像中的位置
    sti=uint8(sti_mask(:,:,3));
    sti(find(sti==12))=255;
    sti=(medfilt2(sti));
    img_edge0=edge(sti);
    B=[0 1 0; 1 1 1; 0 1 0];
    pano_edge=imdilate(img_edge0,B)*255;
    imshow(pano_edge); 

    %% 确定原始图像中的重叠区域
    %  XXX_edge    存放的分别是重叠区域的轮廓
    %  XXX_region  存放的分别是重叠区域,其中重叠区域为0，非重叠区域为1；
    img_name{1}=strcat(strcat(pathL,num2str(100)),'.bmp');
    img_name{2}=strcat(strcat(pathR,num2str(100)),'.bmp');
    img=cell(num_img,1);
    for i=1:num_img
        img{i}=imread( img_name{i} );
    end
    [ml,nl,~]=size(img{1});
    [mr,nr,~]=size(img{2});
    imgL_edge=zeros(ml,nl);
    imgR_edge=zeros(mr,nr);

    for i=1:m
        for j=1:n
            if img_edge0(i,j)==1
                if sti_mask(i,j,3)==1
                    imgL_edge(sti_mask(i,j,1),sti_mask(i,j,2))=1;
                elseif sti_mask(i,j,3)==12
                    imgL_edge(ble_mask(i,j,1),ble_mask(i,j,2))=1;
                end
                if sti_mask(i,j,3)==2
                    imgR_edge(sti_mask(i,j,1),sti_mask(i,j,2))=1;
                elseif sti_mask(i,j,3)==12
                    imgR_edge(ble_mask(i,j,3),ble_mask(i,j,4))=1;
                end            
            end
        end
    end
    imgL_edge=imdilate(imgL_edge,B);
    imgR_edge=imdilate(imgR_edge,B); 
    imgL_region=find_region(imgL_edge,1);
    imgR_region=find_region(imgR_edge,2);
    if 0
        figure,imshow(imgL_region,[]),title('imgL_region');
        figure,imshow(imgR_region,[]),title('imgR_region');
        figure,imshow(imgcat(img{1},imgL_edge*255),[]),title('imgL_edge');
        figure,imshow(imgcat(img{2},imgR_edge*255),[]),title('imgR_edge'); 
    end
%% 测试运动检测程序
    img_name{1}=strcat(strcat(pathL,num2str(100)),'.bmp');
    img_name{2}=strcat(strcat(pathR,num2str(100)),'.bmp');
    imgL00=imread( img_name{1} );
    imgR00=imread( img_name{2} );
    for imgN=150:150
        fprintf( 'imgN =%d \n',imgN);
        img_name{1}=strcat(strcat(pathL,num2str(imgN)),'.bmp');
        img_name{2}=strcat(strcat(pathR,num2str(imgN)),'.bmp');
        img=cell(num_img,1);
        for i=1:num_img
            img{i}=imread( img_name{i} );
        end 
        
        Check_MovingTarget(imgL00,img{1},~imgL_region);
        imgL00=img{1};
        
%         Check_MovingTarget(imgR00,img{2},~imgR_region);
%         imgR00=img{2};
    end
%% 开始主循环
    kk=1;
    last_IMG=[];
    StartNum=140;
    img_name1=strcat(strcat(pathL,num2str(StartNum-1)),'.bmp');
    img_name2=strcat(strcat(pathR,num2str(StartNum-1)),'.bmp');
    imgL00=imread(img_name1);
    imgR00=imread(img_name2);
    for imgN=StartNum:3:300
        fprintf( 'imgN =%d \n ',imgN);
        
        img_name{1}=strcat(strcat(pathL,num2str(imgN)),'.bmp');
        img_name{2}=strcat(strcat(pathR,num2str(imgN)),'.bmp');
        
        img=cell(num_img,1);
        for i=1:num_img
            img{i}=imread( img_name{i} );
        end 
        IMG=uint8(zeros(m,n,3)); 
        tic;
        for i=1:m
            for j=1:n
                switch sti_mask(i,j,3)
                    case 1
                        IMG(i,j,:)=img{1}(sti_mask(i,j,1),sti_mask(i,j,2),:);
                    case 12
                        %IMG(i,j,:)=sti_mask(i,j,1)*img{1}(ble_mask(i,j,1),ble_mask(i,j,2),:)+sti_mask(i,j,2)*img{2}(ble_mask(i,j,3),ble_mask(i,j,4),:);
                        IMG(i,j,:)=0.5*img{1}(ble_mask(i,j,1),ble_mask(i,j,2),:)+0.5*img{2}(ble_mask(i,j,3),ble_mask(i,j,4),:);
                    case 2
                        IMG(i,j,:)=img{2}(sti_mask(i,j,1),sti_mask(i,j,2),:);           
                end
            end
        end
        fprintf( 'imgN =%fs \n ',toc);
        [retL,imgL] = Check_MovingTarget(imgL00,img{1},~imgL_region);
        [retR,imgR] = Check_MovingTarget(imgR00,img{2},~imgR_region);
        imgL00=img{1};
        imgR00=img{2};
        if retL==1 && retR==1 && kk>1 
            ret=Check_MovingTarget(last_IMG,IMG,~uint8((pano_edge/255)));
            figure(20),imshow(ret);            
            pause(0.01); 
        end
        last_IMG=IMG; 
        
        IMG=imgcat(IMG,pano_edge);
        figure;imshow(IMG);         
        kk=kk+1;
%         imwrite(IMG,strcat(strcat('C:/Users/Administrator/Desktop/0702/reslut/ret',num2str(imgN)),'.bmp'));
    end
end
%% 检查重叠区域是否存在运动物体 
function [ret,img]=Check_MovingTarget(new_img0,old_img0,region)
    B=ones(4);     
    new_img=rgb2gray(new_img0).*uint8(region);
    old_img=rgb2gray(old_img0).*uint8(region);    
    img=edge(new_img-old_img);
    img=imdilate(img,B);img=imdilate(img,B);
    img=imerode(img,B);img=imerode(img,B);img=imerode(img,B);
    img=medfilt2(img,[3 3]); img=medfilt2(img,[5 5]); 
    img=imdilate(img,B);    
    if 0
        subplot('Position',[0.02 0.6 0.3 0.3]),imshow(new_img0,[]),title('new img'); 
        subplot('Position',[0.35 0.6 0.3 0.3]),imshow(old_img0,[]),title('old img');
        subplot('Position',[0.68 0.6 0.3 0.3]),imshow(new_img,[]),title('new img * uint8(region)');
        subplot('Position',[0.02 0.1 0.3 0.3]),imshow(old_img,[]),title('old img * uint8(region)');
        subplot('Position',[0.35 0.1 0.3 0.3]),imshow(new_img-old_img,[]),title('new_img - old_img'); 
        subplot('Position',[0.68 0.1 0.3 0.3]),imshow(img,[]),title('result');
    end
    position_w=find(sum(img,1)>0);
    p_w_min=min(position_w);
    p_w_max=max(position_w);
    position_h=find(sum(img,2)>0);
    p_h_min=min(position_h);
    p_h_max=max(position_h);    
    if (p_h_max-p_h_min)*(p_w_max-p_w_min)<200*200
        [m,n]=size(img);
        img_edge=zeros(m,n);
        
        SIZE=100;
        if p_w_min-SIZE<1,p_w_min=1;else p_w_min=p_w_min-SIZE; end
        if p_h_min-SIZE<1,p_h_min=1;else p_h_min=p_h_min-SIZE; end
        if p_w_max+SIZE>n,p_w_max=n;else p_w_max=p_w_max+SIZE; end
        if p_h_max+SIZE>m,p_h_max=m;else p_h_max=p_h_max+SIZE; end
        
        img_edge(p_h_min:p_h_max,p_w_min)=ones(p_h_max-p_h_min+1,1)'*255;
        img_edge(p_h_min:p_h_max,p_w_max)=ones(p_h_max-p_h_min+1,1)'*255;
        img_edge(p_h_min,p_w_min:p_w_max)=ones(p_w_max-p_w_min+1,1)*255;
        img_edge(p_h_max,p_w_min:p_w_max)=ones(p_w_max-p_w_min+1,1)*255;
        img=new_img-old_img+uint8(img_edge);
        ret=1;
    else
        ret=0;
    end 
%      figure(50),imshow(img);
end

%% 计算图像中的重叠区域
function region=find_region(edge,flag)
    [m,n]=size(edge);
    region=ones(m,n);
    Max=max(max(edge));
    if flag==1
       
        for i=1:m  
            k=1;
            tmp=[];
            for j=1:n
                if(edge(i,j)==Max)
                    tmp(k)=j;
                    k=k+1;
                end
            end
            mi=min(tmp);
            ma=max(tmp);
            if mi==0 & ma~=0
                mi=1;
            end
            region(i,mi:n)=zeros(1,n-mi+1);
        end
    end
    if flag==2
        k=1;
        for i=1:m 
            tmp=[];
            for j=1:n
                if(edge(i,j)==Max)
                    tmp(k)=j;
                    k=k+1;
                end
            end
            mi=min(tmp);
            ma=max(tmp);
            if mi==0 & ma~=0
                mi=1;
            end
            region(i,mi:ma)=zeros(1,ma-mi+1);
        end
    end
end

%% 运动区域标记
function ret=imgcat(img,edge)
    edge=uint8(edge);
    img(:,:,1)=img(:,:,1)+edge;
    img(:,:,2)=img(:,:,2)-edge;
    img(:,:,2)=img(:,:,2)-edge;
    ret=img;
end