%%混合润滑模型M-Wway
%%润滑油-SAE30
%中心油膜厚度and4];
clear
close
clc
%中心油膜厚度hc
W1=[2.2e-5 1.1e-4];
for kk=2
W=W1(kk);%载荷系数（线接触)
U=1e-11;deta=2e-5;G=5700;R=0.02;V=0.01;
Hc=2.691*W^(-0.135)*U^(0.705)*G^(0.556)*(1+0.2*deta^1.222*V^0.223*W^(-0.229)*U^(-0.748)*G^(-0.842));
hc=R*Hc
%载荷分布系数
La=0.005*W^(-0.408)*U^(-0.088)*G^0.103*(log(1+4470*deta^6.015*V^1.168*W^0.485*U^(-3.741)*G^(-2.898)))
%温度升高dT
f=[];jj=1;
dV=[];
for S=0:0.01:2;
B=0.05;E=228e9;
%边界润滑系数
fc(jj)=(-0.1+22.28*S)*exp(-181.46*S)+0.1;
ur=0.13;us=S*ur;
k1=60.5;cp1=434;
% k2=47;cp2=460;
k2=0.145;cp2=1880;
den1=7850;den2=888;
lamata=0.091;Z=0.57;mu_o=0.35;Kt=0.045;
F=W*(B*R*E);
L=sqrt(8*R*F/pi/B/E);
b=L;
Pe1=us*b*den1*cp1/2/k1;
Pe2=us*b*den2*cp2/2/k2;
P=F/2/b/B;
Ph=P*(1-La/100);
q=fc(jj)*us*P*La/100+us*lamata*P*(1-La/100);
% A=3*L*q;
% B=sqrt(pi)*(k1*sqrt(1+Pe1)+k2*sqrt(1+Pe2));
% dA=2*L*(fc*P*La/100+lamata*P*(1-La/100));
% dB=sqrt(pi)*(b*den1*cp1/4/sqrt(1+Pe1)+b*den2*cp2/4/sqrt(1+Pe2));
% dT(jj)=ur*(dA*B-A*dB)/(B^2);
dT=2*L*q/(sqrt(pi)*(k1*sqrt(1+Pe1)-k2*sqrt(1+Pe2)));
tau=lamata*P*(1-La/100);
mu_avg=mu_o*exp((log(mu_o)+9.67)*(-1+(1+(5.1e-9)*Ph)^Z)-Kt*dT);
dV(jj)=dT+293;
f(jj)=(La/100)*fc(jj)+(tau/P)*(1-exp(-mu_avg*us/tau/hc));
jj=jj+1;
end
plot(f);
hold on
end

