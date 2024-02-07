clear;
clc;

zmax=3;
zxmin=-zmax;
zxmax=zmax;
zymin=-zmax;
zymax=zmax;
Nzx=1e5;
Nzy=1e5;
x0=linspace(zxmin,zxmax,Nzx);
y0=linspace(zymin,zymax,Nzy);
z0x=zeros(1,Nzx*Nzy);
z0y=zeros(1,Nzx*Nzy);
%%% start points
for i=1:Nzy
z0x(Nzx*(i-1)+1:Nzx*i)=linspace(zxmin,zxmax,Nzx);
z0y(Nzx*(i-1)+1:Nzx*i)=y0(1,i)*ones(1,Nzx);
end
z0=[z0x;z0y];


h=0.001;
stepnum=1;
alpha=1;
sigma=2;
epsilong=0.8;
Nz=Nzx*Nzy;

Deltaonestep=zeros(2,Nz);

%%% Generate Levy noise
M=stblrnd(alpha/2,1,2*(h*cos(pi*alpha/4))^(2/alpha),0,1,Nz);
Normal=randn(2,Nz);
Bh=sqrt(h)*randn(2,Nz);
Levyh=[sigma*sqrt(M).*Normal(1,:);sigma*sqrt(M).*Normal(2,:)];

PK=0.2;
Pm=1.2;
zfx=(z0x-(z0x.^3)/3-z0y)*h+(1+z0y).*Bh(1,:)+Bh(2,:)+Levyh(1,:);
zfy=PK*(Pm+z0x.^2)*h+z0x.*Bh(2,:)+Levyh(2,:);

zf=[zfx;zfy];
M1=size(zf,2);
Deltaonestep=zf;


%%% Identify alpha and sigma
R=sqrt(Deltaonestep(1,:).^2+Deltaonestep(2,:).^2);
q=epsilong;
m=5;
N=2;% 3rounds,n0,n1,n2
Nroundsnum=zeros(1,N+1);% n0,n1,n2

totalnum=size(Deltaonestep,2);

for k=1:N+1
    Nroundpoints=find(R>=q*(m^(k-1)) & R<q*(m^k));
    Nroundsnum(1,k)=size(Nroundpoints,2);
end


nratio=Nroundsnum(1,1)./Nroundsnum(1,2:end);
pos=1:N;
alpha1=log(nratio)./(pos*log(m));
alpha0=sum(alpha1)/N;

cnalpha=alpha0*gamma((2+alpha0)/2)/(2^(1-alpha0)*pi*gamma(1-alpha0/2));
pos=0:N;
sigmak=(q^alpha0*m.^(alpha0*pos).*Nroundsnum*alpha0/(h*Nz*2*pi*cnalpha*(1-m^(-alpha0)))).^(1/alpha0);
sigma0=sum(sigmak)/(N+1);

%%% Identify the drift, R<=epsilong
inepsilong=(R<epsilong);
z02=z0(:,inepsilong);
Deltaonestep2=Deltaonestep(:,inepsilong);

%%% Function dictinary,until cubic polynomials.
%%% 1,x,y,x^2,xy,y^2,x^3,x^2y^1,x^1y^2,y^3
Functionnum=10;
%%% Cross validation?
k=randn(1,size(z02,2));
[m,n]=sort(k);%生成0-Size(data,2)的随机排序记为n
testsetnum=floor(size(z02,2)/10);
z02train=z02(:,n>testsetnum);
z02test=z02(:,n<=testsetnum);

z02xtrain=z02train(1,:);
z02ytrain=z02train(2,:);
z02xtest=z02test(1,:);
z02ytest=z02test(2,:);

Deltaonestep2train=Deltaonestep2(:,n>testsetnum);
Deltaonestep2test=Deltaonestep2(:,n<=testsetnum);

Deltaonestep2xtrain=Deltaonestep2train(1,:);
Deltaonestep2ytrain=Deltaonestep2train(2,:);
Deltaonestep2xtest=Deltaonestep2test(1,:);
Deltaonestep2ytest=Deltaonestep2test(2,:);

M2=size(z02,2);
bix=(M2/(M1*h))*Deltaonestep2train(1,:)';
biy=(M2/(M1*h))*Deltaonestep2train(2,:)';
bixtestreal=(M2/(M1*h))*Deltaonestep2test(1,:)';
biytestreal=(M2/(M1*h))*Deltaonestep2test(2,:)';
%%% for 循环导致计算速度过慢，改写为矩阵形式
%for i=1:M2
%    i
%Atrain(i,:)=[1,z02train(1,i),z02train(2,i),z02train(1,i)^2,z02train(1,i)*z02train(2,i),z02train(2,i)^2,z02train(1,i)^3,z02train(1,i)^2*z02train(2,i),...
%    z02train(1,i)*z02train(2,i)^2,z02train(2,i)^3];
%end

Atrain=[(z02xtrain.^0)',z02xtrain',z02ytrain',(z02xtrain.^2)',(z02xtrain.*z02ytrain)',(z02ytrain.^2)',...
    (z02xtrain.^3)',(z02xtrain.^2.*z02ytrain)',(z02ytrain.^2.*z02xtrain)',(z02ytrain.^3)'];

%for j=1:size(z02test,2)
    
%    Atest(j,:)=[1,z02test(1,j),z02test(2,j),z02test(1,j)^2,z02test(1,j)*z02test(2,j),z02test(2,j)^2,z02test(1,j)^3,z02test(1,j)^2*z02test(2,j),...
%    z02test(1,j)*z02test(2,j)^2,z02test(2,j)^3];

%end
Atest=[(z02xtest.^0)',z02xtest',z02ytest',(z02xtest.^2)',(z02xtest.*z02ytest)',(z02ytest.^2)',...
    (z02xtest.^3)',(z02xtest.^2.*z02ytest)',(z02ytest.^2.*z02xtest)',(z02ytest.^3)'];

%%% identify bx
steperrorbx=999999999*ones(1,Functionnum);
indexnumbx=1:1:10;
stepcoefficientbx=zeros(Functionnum,Functionnum);%%%第i列表示第i次稀疏迭代中原函数基的系数

for i=1:Functionnum

AAA1=Atrain(:,indexnumbx);
AAA2=Atest(:,indexnumbx);
Cx=inv(AAA1'*AAA1)*AAA1'*bix;
bixtest=AAA2*Cx;
steperrorbx(1,i)=sqrt(sum((bixtestreal-bixtest).^2))/size(z02test,2);
stepcoefficientbx(indexnumbx,i)=Cx;
deletenum=find(abs(Cx)==min(abs(Cx)));
indexnumbx(:,deletenum)=[];

end

bestiteratenum=find(steperrorbx==min(steperrorbx));%寻找对应误差最小的迭代次数
cofficientbx=stepcoefficientbx(:,bestiteratenum);

%%% identify by
steperrorby=999999999*ones(1,Functionnum);
indexnumby=1:1:10;
stepcoefficientby=zeros(Functionnum,Functionnum);%%%第i列表示第i次稀疏迭代中原函数基的系数

for i=1:Functionnum

AAA1=Atrain(:,indexnumby);
AAA2=Atest(:,indexnumby);
Cy=inv(AAA1'*AAA1)*AAA1'*biy;
biytest=AAA2*Cy;
steperrorby(1,i)=sqrt(sum((biytestreal-biytest).^2))/size(z02test,2);
stepcoefficientby(indexnumby,i)=Cy;
deletenum=find(abs(Cy)==min(abs(Cy)));
indexnumby(:,deletenum)=[];

end

bestiteratenum=find(steperrorby==min(steperrorby));%寻找对应误差最小的迭代次数
cofficientby=stepcoefficientby(:,bestiteratenum);

%%% Identify the diffusion term/the covariance matrix
Bxx=(M2/(M1*h))*Deltaonestep2xtrain.^2-pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
Byy=(M2/(M1*h))*Deltaonestep2ytrain.^2-pi*sigma0^alpha0*epsilong^(2-alpha0)*cnalpha/(2-alpha0);
Bxy=(M2/(M1*h))*Deltaonestep2ytrain.*Deltaonestep2xtrain;
Bxx=Bxx';
Byy=Byy';
Bxy=Bxy';

dxx=inv(Atrain'*Atrain)*Atrain'*Bxx;
dyy=inv(Atrain'*Atrain)*Atrain'*Byy;
dxy=inv(Atrain'*Atrain)*Atrain'*Bxy;