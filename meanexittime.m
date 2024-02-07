clear;
clc;

r=-0.2;
theta=acos((1+r)/(1-r));
pitheta=theta/pi;
J=2000;
N=2*J+1;% index number of x=0 is (N+1)/2
meshx=linspace(-pi,pi,N);

meshxNum=(1:1:N)-(J+1)*ones(1,N);

h=meshx(1,2)-meshx(1,1);
numthetapositive=(N+1)/2+floor(theta/h)+1;
num0=N+1;
%Numthetapositive=numthetapositive-(J+1);
% numthetanegative=(N+1)/2-floor(theta/h)-1;
% meshx(numthetapositive)
% meshx(numthetanegative)

% solve KU=F
K=zeros(N,N);
F=zeros(N,1);
sigmaB=0.5; % intensity of Bt
alpha=1;
cnalpha=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*(pi^(1/2))*gamma(1-alpha/2));
Drift=@(x) (1-cos(x))+r*(1+cos(x))
Diffusion=@(x) sigmaB*(1+cos(x))
Driftx=Drift(meshx);
Diffusionx=Diffusion(meshx);


%% K(1,1)=1;
for i=1:10
  K(i,i)=1;
end


for i=numthetapositive:N
  K(i,i)=1;
end

for i=11:1:numthetapositive-1

  F(i,1)=-1;

  K(i,i+1)=0.5*(1/h)*Driftx(1,i)+0.5*(1/h^2)*Diffusionx(1,i)^2;
  K(i,i-1)=0.5*(1/h)*Driftx(1,i)+0.5*(1/h^2)*Diffusionx(1,i)^2;
  K(i,i)=-(1/h^2)*Diffusionx(1,i)^2;

  for k=1-i:1:numthetapositive-i

    if k~=0
      
      K(i,i+k)=K(i,i+k)+cnalpha*h*(1/abs(meshx(1,i+k)-meshx(1,i))^(1+alpha));
      K(i,i)=K(i,i)-cnalpha*h*(1/abs(meshx(1,i+k)-meshx(1,i))^(1+alpha));

    end

  end

end

U=inv(K)*F;

plot(meshx,U);