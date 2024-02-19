clear;
clc;

r=-0.05;
theta=acos((1+r)/(1-r));
pitheta=theta/pi;
J=2000;
N=2*J+1;% index number of x=0 is (N+1)/2
meshx=linspace(-pi,pi,N);
% meshxNum=(1:1:N)-(J+1)*ones(1,N);
h=meshx(1,2)-meshx(1,1);
numthetapositive=(N+1)/2+floor(theta/h)+1;
numthetanegative=(N-1)/2-floor(theta/h)+1;
exitregion=[1001,numthetapositive];
A=meshx(1,exitregion(1,1));
B=meshx(1,exitregion(1,2));
nA=exitregion(1,1);
nB=exitregion(1,2);

meanexittimeB0=[meshx];
alphaB0=0.95:0.05:1.3;

for alpha=0.95:0.05:1.3

% solve KU=F
K=zeros(N,N);
F=zeros(N,1); % mean exit time F
F2=zeros(N,1); % exit distribution
sigmaB=0; % intensity of Bt
%alpha=1;
sigmaL=1;
k=sigmaL^alpha; % k:the jump frequency coefficient
cnalpha=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*(pi^(1/2))*gamma(1-alpha/2));
Drift=@(x) (1-cos(x))+r*(1+cos(x))
Diffusion=@(x) sigmaB*(1+cos(x))
Driftx=Drift(meshx);
Diffusionx=Diffusion(meshx);

for i=1:nA

  K(i,i)=1;

end

for i=nB:N

  K(i,i)=1;
  F2(i,1)=1;

end

for i=nA+1:nB-1

  %%显示进度
  [i,numthetapositive-1]

  F(i,1)=-1;
  K(i,i+1)=K(i,i+1)+Driftx(1,i)/(2*h);
  K(i,i-1)=K(i,i+1)-Driftx(1,i)/(2*h);

  K(i,i+1)=K(i,i+1)+(1/2)*(1/h^2)*Diffusionx(1,i)^2;
  K(i,i-1)=K(i,i-1)+(1/2)*(1/h^2)*Diffusionx(1,i)^2;
  K(i,i)=K(i,i)-(1/h^2)*Diffusionx(1,i)^2;

  %%K(i,i)=K(i,i)-(k*cnalpha/alpha)*(1/(-pi+meshx(1,i))^alpha+1/(pi-meshx(1,i))^alpha);
  K(i,i)=K(i,i)-(k*cnalpha/alpha)*(1/(meshx(1,i)-A)^alpha+1/(B-meshx(1,i))^alpha);

  for j=nA-i:1:nB-i

    if j~=0

      K(i,i+j)=K(i,i+j)+k*cnalpha/(abs(meshx(1,i+j)-meshx(1,i)))^(1+alpha);
      K(i,i)=K(i,i)-k*cnalpha/(abs(meshx(1,i+j)-meshx(1,i)))^(1+alpha);

    end

  end

  %F2(i,1)=-1*k*cnalpha*(1/alpha)*(1/(theta-meshx(1,i))^alpha);

end

U=inv(K)*F;
%P=inv(K)*F2;
meanexittimeB0=[meanexittimeB0;U'];

end

%%save('B0.mat','sigmaB','sigmaL','alphaB0','meanexittimeB0');
subplot(1,2,1);
for i=2:size(meanexittimeB0,1)

set(0,'defaultfigurecolor','w');
plot(meanexittimeB0(1,:),meanexittimeB0(i,:),'linewidth',2);
grid on
hold on

end
legend('α = 0.95','α = 0.1','α = 1.05','α = 1.1','α = 1.15','α = 1.2','α = 1.25','α = 1.3','AutoUpdate','off','FontSize',12,'FontName','Times New Roman');
xlim([-pi,pi]);
xlabel('Action potential   θ','FontSize',14,'FontName','Times New Roman');
ylabel('Mean first exit time   u','FontSize',14,'FontName','Times New Roman');
text(-2.7,0.09,'(a)','FontSize',20,'FontName','Times New Roman');
hold on
scatter([-pi/2],[0],'g','filled');
hold on
scatter([-theta],[0],'r','filled');
hold on
scatter([theta],[0],'g','filled');
text(-0.5,-0.004,'θ-','FontSize',12,'FontName','Times New Roman');
text(0.4,-0.004,'θ+','FontSize',12,'FontName','Times New Roman');
text(-pi/2-0.2,-0.004,'-π/2','FontSize',12,'FontName','Times New Roman');


subplot(1,2,2);
exittimelocal=meanexittimeB0(2:size(meanexittimeB0,1),nA+1:nB-1);
exittimelocal=1./exittimelocal;
frequencyB0=[meshx(1,nA+1:nB-1);exittimelocal];

for i=2:size(frequencyB0,1)

set(0,'defaultfigurecolor','w');
plot(frequencyB0(1,:),frequencyB0(i,:),'linewidth',2);
grid on
hold on

end
legend('α = 0.95','α = 0.1','α = 1.05','α = 1.1','α = 1.15','α = 1.2','α = 1.25','α = 1.3','AutoUpdate','off','FontSize',12,'FontName','Times New Roman');
xlim([-pi,pi]);
ylim([0,1000]);
xlabel('Action potential   θ','FontSize',14,'FontName','Times New Roman');
ylabel('Mean spiking frequency   ω','FontSize',14,'FontName','Times New Roman');
text(-2.7,900,'(b)','FontSize',20,'FontName','Times New Roman');
hold on
scatter([-pi/2],[0],'g','filled');
hold on
scatter([-theta],[0],'r','filled');
hold on
scatter([theta],[0],'g','filled');

text(-0.5,-40,'θ-','FontSize',12,'FontName','Times New Roman');
text(0.4,-40,'θ+','FontSize',12,'FontName','Times New Roman');
text(-pi/2-0.2,-40,'-π/2','FontSize',12,'FontName','Times New Roman');
