clear;
clc;

xmax=pi;
xmin=-pi;
Nzx=1e6;
x0=linspace(xmin,xmax,Nzx);


%%% start points
z0=x0;


h=0.001;
stepnum=1;
alpha=1;
sigma=2;
sigmaB=0.5;
epsilong=0.2;
r=0.2;
Nz=Nzx;

Deltaonestep=zeros(1,Nz);

%%% Generate Levy noise
M=stblrnd(alpha/2,1,2*(h*cos(pi*alpha/4))^(2/alpha),0,1,Nz);
Normal=randn(1,Nz);
Bh=sqrt(h)*randn(1,Nz);
Levyh=[sigma*sqrt(M).*Normal(1,:)];


%zfx=(z0x-z0x.^3-5*z0x.*z0y.^2)*h+(1+z0y).*Bh(1,:)+Bh(2,:)+Levyh(1,:);
zfx=((1-cos(z0))+(1+cos(z0))*r)*h+sigmaB*(1+cos(1*z0)).*Bh(1,:)+Levyh(1,:);


zf=[zfx];
M1=size(zf,2);
Deltaonestep=zf;


%%% Identify alpha and sigma
R=sqrt(Deltaonestep(1,:).^2);
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

cnalpha=alpha0*gamma((1+alpha0)/2)/(2^(1-alpha0)*(pi^(1/2))*gamma(1-alpha0/2));
pos=0:N;
sigmak=(q^alpha0*m.^(alpha0*pos).*Nroundsnum*alpha0/(h*Nz*2*cnalpha*(1-m^(-alpha0)))).^(1/alpha0);
sigma0=sum(sigmak)/(N+1);

%%% Identify the drift, R<=epsilong
inepsilong=(R<epsilong);
z02=z0(:,inepsilong);
Deltaonestep2=Deltaonestep(:,inepsilong);

%%% Function dictinary,until cubic polynomials and trigonometric functions.
%%% 1,x,x^2,x^3,sinx,cosx,sin2x,cos2x,sin3x,cos3x
Functionnum=10;
%%% Cross validation
k=randn(1,size(z02,2));
[m,n]=sort(k);%生成0-Size(data,2)的随机排序记为n
testsetnum=floor(size(z02,2)/10);
z02train=z02(:,n>testsetnum);
z02test=z02(:,n<=testsetnum);

z02xtrain=z02train(1,:);
%z02ytrain=z02train(2,:);
z02xtest=z02test(1,:);
%z02ytest=z02test(2,:);

Deltaonestep2train=Deltaonestep2(:,n>testsetnum);
Deltaonestep2test=Deltaonestep2(:,n<=testsetnum);

Deltaonestep2xtrain=Deltaonestep2train(1,:);
%Deltaonestep2ytrain=Deltaonestep2train(2,:);
Deltaonestep2xtest=Deltaonestep2test(1,:);
%Deltaonestep2ytest=Deltaonestep2test(2,:);

M2=size(z02,2);
bix=(M2/(M1*h))*Deltaonestep2train(1,:)';
%biy=(M2/(M1*h))*Deltaonestep2train(2,:)';
bixtestreal=(M2/(M1*h))*Deltaonestep2test(1,:)';
%biytestreal=(M2/(M1*h))*Deltaonestep2test(2,:)';


Atrain=[(z02xtrain.^0)',z02xtrain',(z02xtrain.^2)',(z02xtrain.^3)',(sin(z02xtrain))',(cos(z02xtrain))',...
    (sin(2*z02xtrain))',(cos(2*z02xtrain))',(sin(3*z02xtrain))',(cos(3*z02xtrain))'];

Atest=[(z02xtest.^0)',z02xtest',(z02xtest.^2)',(z02xtest.^3)',(sin(z02xtest))',(cos(z02xtest))',...
    (sin(2*z02xtest))',(cos(2*z02xtest))',(sin(3*z02xtest))',(cos(3*z02xtest))'];

steperrorbx=999999999*ones(1,Functionnum);
indexnum=1:1:Functionnum;
stepcoefficient=zeros(Functionnum,Functionnum);%%%第i列表示第i次稀疏迭代中原函数基的系数

for i=1:Functionnum

AAA1=Atrain(:,indexnum);
AAA2=Atest(:,indexnum);
Cx=inv(AAA1'*AAA1)*AAA1'*bix;
bixtest=AAA2*Cx;
steperrorbx(1,i)=sqrt(sum((bixtestreal-bixtest).^2))/size(z02test,2);
stepcoefficient(indexnum,i)=Cx;
deletenum=find(abs(Cx)==min(abs(Cx)));
indexnum(:,deletenum)=[];

end

bestiteratenum=find(steperrorbx==min(steperrorbx))%寻找对应误差最小的迭代次数
cofficient=stepcoefficient(:,bestiteratenum)

%%% identify the covariance matrix
omiga11x=(M2/(M1*h))*(Deltaonestep2train(1,:).^2)';
omiga11xtestreal=(M2/(M1*h))*(Deltaonestep2test(1,:).^2)';
s11x=2*(sigma0^alpha0)*cnalpha*(1/(2-alpha0))*epsilong^(2-alpha0);
B11x=omiga11x-s11x;
B11xtestreal=omiga11xtestreal-s11x;
%%% basis:1,x,x^2,sinx,cosx,cosx^2
Functionnum2=6;

A2train=[(z02xtrain.^0)',(z02xtrain.^1)',(z02xtrain.^2)',...
    (sin(z02xtrain))',(cos(z02xtrain))',(cos(z02xtrain).^2)'];

A2test=[(z02xtest.^0)',(z02xtest.^1)',(z02xtest.^2)',...
    (sin(z02xtest))',(cos(z02xtest))',(cos(z02xtest).^2)'];

steperrorbx2=999999999*ones(1,Functionnum2);
indexnum2=1:1:Functionnum2;
stepcoefficient2=zeros(Functionnum2,Functionnum2);%%%第i列表示第i次稀疏迭代中原函数基的系数

for i=1:Functionnum2

AAA1=A2train(:,indexnum2);
AAA2=A2test(:,indexnum2);
Cx=inv(AAA1'*AAA1)*AAA1'*B11x;
B11xtest=AAA2*Cx;
steperrorbx2(1,i)=sqrt(sum((B11xtestreal-B11xtest).^2))/size(z02test,2);
stepcoefficient2(indexnum2,i)=Cx;
deletenum2=find(abs(Cx)==min(abs(Cx)));
indexnum2(:,deletenum2)=[];

end

bestiteratenum2=find(steperrorbx2==min(steperrorbx2))%寻找对应误差最小的迭代次数
cofficient2=stepcoefficient2(:,bestiteratenum2)

%%% save('r=-0.2.mat','stepcoefficient','bestiteratenum','cofficient','stepcoefficient2','bestiteratenum2','cofficient2');