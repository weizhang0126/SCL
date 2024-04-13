%N=128 CA-SCL
clc
clear
close all
addpath('GA/');
addpath('code/');
n=7;
N=2^n;
L=32;
crc_length=8;
K1=N/2;
K=K1+crc_length;
SNR=2;
R=1/2;
ListSNR=1:0.4:3;
ts=10000;
Len=length(ListSNR);
sigma=10^(-SNR/20);
%% GA构造
[channels, ~] = GA(sigma, N);
[~, channel_ordered] = sort(channels, 'descend');%降序
info_bits = sort(channel_ordered(1 : K), 'ascend');
frozen_bits = ones(N , 1);
frozen_bits(info_bits) = 0;%把1当成冻结比特
info_bits_logical = logical(mod(frozen_bits + 1, 2))';
CBR=double(info_bits_logical);
xx=CBR;
BLER=zeros(1,Len);
for i=1:Len
    tic
    count=0;
    snr=ListSNR(i);
    sigma=1/(sqrt(2*R))*10^(-snr/20);
    for tt=1:ts
    u=randi([0,1],1,K1);
    g=[1,1,1,1,1,1,1,0,0];%x8
    [G,H]=crc_generator_matrix(g,K1);
    crc_check=G(:,K1+1:end);
    u_crc=mod(u*crc_check,2);
    %拼接
    uu=[u,u_crc];
    %信息位
    for kk=1:K
        xx(info_bits(kk))=uu(kk);
    end
    %编码，调制，加噪，转换成对数域
    u1=encode(xx);
    u2=1-2*u1;
    n1=randn(1,N);
    nn=n1*sigma;
    y=u2+nn;
    llr=2/(sigma^2)*y;
    ud=SCL(llr,CBR,N,n,L);
    u_d=zeros(K,L);
    flag=-1;
    for l=1:L
        for kk=1:K
            u_d(kk,l)=ud(info_bits(kk),l);
        end
         %第一分段通过CRC
        if sum(mod(H*u_d(:,l),2))==0
            flag=l;
            break;
        end
    end
    if flag==-1
        count=count+1;
    end
    end
    %计算误帧率
    BLER(i)=count/ts;
    toc
end
semilogy(ListSNR,BLER,'b-x',LineWidth=1);
xlabel("EbNo");
ylabel("BLER");
legend("CA-SCL");
title("FER仿真结果");