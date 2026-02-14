clc; clear; close all;

% ---- 固定参数 ----
% tR = sqrt(3);
tL =20/10;
beta = pi/4;
alpha=pi/3;
delta1=1;
delta2_min = -5; delta2_max = 5;
tR_min = -3; tR_max   =6;
N1 = 2000; N2 = 2000;   

[delta2, tR] = meshgrid(linspace(delta2_min, delta2_max, N1), ...
                        linspace(tR_min, tR_max, N2));


C =-2 * tL .* tR .* cos(alpha - beta) + 2 * delta2.^2+delta1.^2 ; 
C0 =-2 * tL .* tR .* cos(alpha - beta) + 2 * delta2.^2;
D=2*sqrt((tR.^2 -  delta2.^2).*(tL.^2 - delta2.^2));
mu= tL.^2.*cos(alpha).^2./(tL.^2 - delta2.^2) + tR.^2*cos(beta).^2./(tR.^2 - delta2.^2) - 1;
kappa= (cos(alpha).*tL./(tL.^2 - delta2.^2) + cos(beta).*tR./(tR.^2 - delta2.^2)).*delta2.*delta1;
eta = delta2.^2.*delta1.^2./(tR.^2 - delta2.^2) + delta2.^2.*delta1.^2./(tL.^2 - delta2.^2);
Xe = (-2*delta2.^2 + 2*cos(alpha - beta).*tL.*tR)./D;
Ee = -delta1.*tL.*tR.*sin(alpha - beta)./(delta2.*(sin(alpha).*tL - sin(beta).*tR));
U = 4*(Ee.*cos(beta).*tR + delta2.*delta1).*(Ee.*tL.*cos(alpha) + delta2.*delta1)./D.^2 - 1;
GAMMA = Ee.^2.*mu./D + 2*Ee.*kappa./D + eta./D + C./D;

Delta = (U.*Xe - GAMMA).^2 - 4*Xe.^3.*GAMMA;
Y = 1/2.*(U.*Xe - GAMMA + sqrt(Delta))./Xe.^2;
Z = 1/2.*(U.*Xe - GAMMA - sqrt(Delta))./Xe.^2;

W1 = @(x) x + sqrt(x.^2 - 1);
W2 = @(x) x - sqrt(x.^2 - 1);

W  = @(x) max(abs(W1(x)), abs(W2(x)));
W_Xe = W(Xe);

Wc_Xe = W_Xe;
W_Y = W(Y);
W_Z = W(Z);

Wm_YZ = max(W_Y, W_Z);
 figure
 

[zxy, h] = contourf(delta2, tR,  Wm_YZ - Wc_Xe, [0 0]);

hold on;box on; grid on;
sitez0=find(zxy(1,:)==0);
N1=zxy(2,sitez0(1));N2=zxy(2,sitez0(2));N3=zxy(2,sitez0(3));N4=zxy(2,sitez0(4));
N5=zxy(2,sitez0(5));N6=zxy(2,sitez0(6));N7=zxy(2,sitez0(7));
N8=zxy(2,sitez0(8));N9=zxy(2,sitez0(9));N10=zxy(2,sitez0(10));
N11=zxy(2,sitez0(11));N12=zxy(2,sitez0(12));N13=zxy(2,sitez0(13));
N14=zxy(2,sitez0(14));N15=zxy(2,sitez0(15));N16=zxy(2,sitez0(16));
d21=zxy(1,sitez0(1)+1:sitez0(1)+N1);d22=zxy(1,sitez0(2)+1:sitez0(2)+N2);
d23=zxy(1,sitez0(3)+1:sitez0(3)+N3);d24=zxy(1,sitez0(4)+1:sitez0(4)+N4);
d25=zxy(1,sitez0(5)+1:sitez0(5)+N5);d26=zxy(1,sitez0(6)+1:sitez0(6)+N6);
d27=zxy(1,sitez0(7)+1:sitez0(7)+N7);d28=zxy(1,sitez0(8)+1:sitez0(8)+N8);
d29=zxy(1,sitez0(9)+1:sitez0(9)+N9);d210=zxy(1,sitez0(10)+1:sitez0(10)+N10);
d211=zxy(1,sitez0(11)+1:sitez0(11)+N11);d212=zxy(1,sitez0(12)+1:sitez0(12)+N12);
d213=zxy(1,sitez0(13)+1:sitez0(13)+N13);d214=zxy(1,sitez0(14)+1:sitez0(14)+N14);
d215=zxy(1,sitez0(15)+1:sitez0(15)+N15);d216=zxy(1,sitez0(16)+1:sitez0(16)+N16);
tr1=zxy(2,sitez0(1)+1:sitez0(1)+N1);tr2=zxy(2,sitez0(2)+1:sitez0(2)+N2);
tr3=zxy(2,sitez0(3)+1:sitez0(3)+N3);tr4=zxy(2,sitez0(4)+1:sitez0(4)+N4);
tr5=zxy(2,sitez0(5)+1:sitez0(5)+N5);tr6=zxy(2,sitez0(6)+1:sitez0(6)+N6);
tr7=zxy(2,sitez0(7)+1:sitez0(7)+N7);tr8=zxy(2,sitez0(8)+1:sitez0(8)+N8);
tr9=zxy(2,sitez0(9)+1:sitez0(9)+N9);tr10=zxy(2,sitez0(10)+1:sitez0(10)+N10);
tr11=zxy(2,sitez0(11)+1:sitez0(11)+N11);tr12=zxy(2,sitez0(12)+1:sitez0(12)+N12);
tr13=zxy(2,sitez0(13)+1:sitez0(13)+N13);tr14=zxy(2,sitez0(14)+1:sitez0(14)+N14);
tr15=zxy(2,sitez0(15)+1:sitez0(15)+N15);tr16=zxy(2,sitez0(16)+1:sitez0(16)+N16);


fillcolor='r';
fill([d21  d21(end) d21(1)],[tr1 tr1(1) tr1(1)],fillcolor,FaceAlpha='0.4');
fill([d22  d22(end) d22(1)],[tr2 tr2(end)  tr2(end)],fillcolor,FaceAlpha='0.4');
fill([d23  d23(1) ],[tr3 tr3(end)],fillcolor,FaceAlpha='0.4');
fill([d24  d24(end) ],[tr4 tr4(1)],fillcolor,FaceAlpha='0.4');
fill([d25 ],[tr5 ],fillcolor,FaceAlpha='0.4');
fill([d26 ],[tr6 ],fillcolor,FaceAlpha='0.4');
fill([d29 ],[tr9 ],fillcolor,FaceAlpha='0.4');
fill([d212 ],[tr12 ],fillcolor,FaceAlpha='0.4');
xlabel('$\delta_2$','Interpreter','latex', 'FontSize', 14);
ylabel('$t_R$','Interpreter','latex' ,'FontSize', 14);
xticks([ -4  -2  0  2  4 ]);

