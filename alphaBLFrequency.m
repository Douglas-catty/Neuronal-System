%sigmaB=[0,2]
%alpha=[0.95,1,1.1,1.5]
%sigmaL=[1,2]
r=-0.05;
theta=acos((1+r)/(1-r));
pitheta=theta/pi;
J=2000;
N=2*J+1;% index number of x=0 is (N+1)/2
meshx=linspace(-pi,pi,N);
h=meshx(1,2)-meshx(1,1);
numthetapositive=(N+1)/2+floor(theta/h)+1;
numthetanegative=(N-1)/2-floor(theta/h)+1;
exitregion=[1001,numthetapositive];
A=meshx(1,exitregion(1,1));
B=meshx(1,exitregion(1,2));
nA=exitregion(1,1);
nB=exitregion(1,2);

meanexittimealpha10=[];
meanexittimealpha13=[];
meanexittimealpha16=[];
meanexittimealpha19=[];

alpha=0.95;

for sigmaB=0:0.1:2
    [sigmaB,2]
    for sigmaL=1:0.1:2
% solve KU=F
K=zeros(N,N);
F=zeros(N,1); % mean exit time F

%sigmaB=0; % intensity of Bt
%sigmaL=1;

k=sigmaL^alpha; % k:the jump frequency coefficient
cnalpha=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*(pi^(1/2))*gamma(1-alpha/2));
Drift=@(x) (1-cos(x))+r*(1+cos(x));
Diffusion=@(x) sigmaB*(1+cos(x));
Driftx=Drift(meshx);
Diffusionx=Diffusion(meshx);

for i=1:nA

  K(i,i)=1;

end

for i=nB:N

  K(i,i)=1;
 
end

for i=nA+1:nB-1


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

end

U=inv(K)*F;
list=[sigmaB,sigmaL,U(numthetanegative,1)];
meanexittimealpha10=[meanexittimealpha10;list];

   end
end


alpha=1;

for sigmaB=0:0.1:2
    [sigmaB,2]
    for sigmaL=1:0.1:2
% solve KU=F
K=zeros(N,N);
F=zeros(N,1); % mean exit time F

%sigmaB=0; % intensity of Bt
%sigmaL=1;

k=sigmaL^alpha; % k:the jump frequency coefficient
cnalpha=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*(pi^(1/2))*gamma(1-alpha/2));
Drift=@(x) (1-cos(x))+r*(1+cos(x));
Diffusion=@(x) sigmaB*(1+cos(x));
Driftx=Drift(meshx);
Diffusionx=Diffusion(meshx);

for i=1:nA

  K(i,i)=1;

end

for i=nB:N

  K(i,i)=1;
 
end

for i=nA+1:nB-1


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

end

U=inv(K)*F;
list=[sigmaB,sigmaL,U(numthetanegative,1)];
meanexittimealpha13=[meanexittimealpha13;list];

   end
end


alpha=1.05;

for sigmaB=0:0.1:2
    [sigmaB,2]
    for sigmaL=1:0.1:2
% solve KU=F
K=zeros(N,N);
F=zeros(N,1); % mean exit time F

%sigmaB=0; % intensity of Bt
%sigmaL=1;

k=sigmaL^alpha; % k:the jump frequency coefficient
cnalpha=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*(pi^(1/2))*gamma(1-alpha/2));
Drift=@(x) (1-cos(x))+r*(1+cos(x));
Diffusion=@(x) sigmaB*(1+cos(x));
Driftx=Drift(meshx);
Diffusionx=Diffusion(meshx);

for i=1:nA

  K(i,i)=1;

end

for i=nB:N

  K(i,i)=1;
 
end

for i=nA+1:nB-1


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

end

U=inv(K)*F;
list=[sigmaB,sigmaL,U(numthetanegative,1)];
meanexittimealpha16=[meanexittimealpha16;list];

   end
end



alpha=1.1;

for sigmaB=0:0.1:2
    [sigmaB,2]
    for sigmaL=1:0.1:2
% solve KU=F
K=zeros(N,N);
F=zeros(N,1); % mean exit time F

%sigmaB=0; % intensity of Bt
%sigmaL=1;

k=sigmaL^alpha; % k:the jump frequency coefficient
cnalpha=alpha*gamma((1+alpha)/2)/(2^(1-alpha)*(pi^(1/2))*gamma(1-alpha/2));
Drift=@(x) (1-cos(x))+r*(1+cos(x));
Diffusion=@(x) sigmaB*(1+cos(x));
Driftx=Drift(meshx);
Diffusionx=Diffusion(meshx);

for i=1:nA

  K(i,i)=1;

end

for i=nB:N

  K(i,i)=1;
 
end

for i=nA+1:nB-1


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

end

U=inv(K)*F;
list=[sigmaB,sigmaL,U(numthetanegative,1)];
meanexittimealpha19=[meanexittimealpha19;list];

   end
end

frequencyalpha10=meanexittimealpha10;
frequencyalpha10(:,3)=1./frequencyalpha10(:,3);
frequencyalpha13=meanexittimealpha13;
frequencyalpha13(:,3)=1./frequencyalpha13(:,3);
frequencyalpha16=meanexittimealpha16;
frequencyalpha16(:,3)=1./frequencyalpha16(:,3);
frequencyalpha19=meanexittimealpha19;
frequencyalpha19(:,3)=1./frequencyalpha19(:,3);


%%% plot
subplot(4,2,1)
x1=0:0.1:2;
y1=1:0.1:2;
[x2,y2]=meshgrid(x1,y1);
z2=griddata(meanexittimealpha10(:,1),meanexittimealpha10(:,2),meanexittimealpha10(:,3)...
    ,x2,y2,'v4');
set(0,'defaultfigurecolor','w');
surf(x2,y2,z2);
shading flat
shading interp
hold on
C=contour(x2,y2,z2);
colormap('jet')
colorbar
view(0,90);
xlabel('Intensity of Gaussian white noise σB','FontSize',14,'FontName','Times New Roman');
ylabel('Intensity of Lévy noise σL','FontSize',14,'FontName','Times New Roman');

subplot(4,2,2)
z2=griddata(frequencyalpha10(:,1),frequencyalpha10(:,2),frequencyalpha10(:,3)...
    ,x2,y2,'v4');
set(0,'defaultfigurecolor','w');
surf(x2,y2,z2);
shading flat
shading interp
hold on
C=contour(x2,y2,z2);
colormap('jet')
colorbar
view(0,90);

subplot(4,2,3)
z2=griddata(meanexittimealpha13(:,1),meanexittimealpha13(:,2),meanexittimealpha13(:,3)...
    ,x2,y2,'v4');
set(0,'defaultfigurecolor','w');
surf(x2,y2,z2);
shading flat
shading interp
hold on
C=contour(x2,y2,z2);
colormap('jet')
colorbar
view(0,90);

subplot(4,2,4)
z2=griddata(frequencyalpha13(:,1),frequencyalpha13(:,2),frequencyalpha13(:,3)...
    ,x2,y2,'v4');
set(0,'defaultfigurecolor','w');
surf(x2,y2,z2);
shading flat
shading interp
hold on
C=contour(x2,y2,z2);
colormap('jet')
colorbar
view(0,90);

subplot(4,2,5)
z2=griddata(meanexittimealpha16(:,1),meanexittimealpha16(:,2),meanexittimealpha16(:,3)...
    ,x2,y2,'v4');
set(0,'defaultfigurecolor','w');
surf(x2,y2,z2);
shading flat
shading interp
hold on
C=contour(x2,y2,z2);
colormap('jet')
colorbar
view(0,90);

subplot(4,2,6)
z2=griddata(frequencyalpha16(:,1),frequencyalpha16(:,2),frequencyalpha16(:,3)...
    ,x2,y2,'v4');
set(0,'defaultfigurecolor','w');
surf(x2,y2,z2);
shading flat
shading interp
hold on
C=contour(x2,y2,z2);
colormap('jet')
colorbar
view(0,90);

subplot(4,2,7)
z2=griddata(meanexittimealpha19(:,1),meanexittimealpha19(:,2),meanexittimealpha19(:,3)...
    ,x2,y2,'v4');
set(0,'defaultfigurecolor','w');
surf(x2,y2,z2);
shading flat
shading interp
hold on
C=contour(x2,y2,z2);
colormap('jet')
colorbar
view(0,90);

subplot(4,2,8)
z2=griddata(frequencyalpha19(:,1),frequencyalpha19(:,2),frequencyalpha19(:,3)...
    ,x2,y2,'v4');
set(0,'defaultfigurecolor','w');
surf(x2,y2,z2);
shading flat
shading interp
hold on
C=contour(x2,y2,z2);
colormap('jet')
colorbar
view(0,90);