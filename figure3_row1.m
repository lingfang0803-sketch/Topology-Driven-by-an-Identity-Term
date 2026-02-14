clc; clear; close all;


tL =20/10;
alpha=pi/3;
beta = pi/4;
delta1=0;

delta2_min = -4; delta2_max = 4;
tR_min = -2; tR_max   =4;
N1 = 2000; N2 = 2000;   

[delta2, tR] = meshgrid(linspace(delta2_min, delta2_max, N1), ...
                        linspace(tR_min, tR_max, N2));



C0 =-2 * tL .* tR .* cos(alpha - beta) + 2 * delta2.^2;
D=2*sqrt((tR.^2 -  delta2.^2).*(tL.^2 - delta2.^2));

subplot(2,3,1);hold on; box on;grid on;


 
[zxy,h2] = contour(delta2, tR, D.^2-C0.^2, [0 0], '-','LineWidth', 2,'LineColor', [0.4940 0.1840 0.5560]);
h2.DisplayName = 'C^2-D^2,';

sitez0=find(zxy(1,:)==0);
N1=zxy(2,sitez0(1));N2=zxy(2,sitez0(2));
d21=zxy(1,sitez0(1)+1:sitez0(1)+N1);d210=find(d21==min(abs(d21)));
tr1=zxy(2,sitez0(1)+1:sitez0(1)+N1);
d22=zxy(1,sitez0(2)+1:sitez0(2)+N2);d220=find(d22==min(abs(d22)));
tr2=zxy(2,sitez0(2)+1:sitez0(2)+N2);


fill(d21,tr1,'c',FaceAlpha='0.7');
fill(d22,tr2,'c',FaceAlpha='0.7');
fill([delta2_max d21(1:d210) d22(d220:end) delta2_max],[tR_min tr1(1:d210) tr2(d220:end) tR_max],'g',FaceAlpha='0.3');%右下到右上
fill([delta2_min d22(1:d220) d21(d210:end) delta2_min],[tR_max tr2(1:d220) tr1(d210:end) tR_min],'g',FaceAlpha='0.3');%左上到左下
xlabel('$\delta_2$','Interpreter','latex', 'FontSize', 14);
ylabel('$t_R$','Interpreter','latex' ,'FontSize', 14);
plot(d21,tr1,d22,tr2,'LineWidth',2,'Color',[0.4940 0.1840 0.5560]);



tR_fixed =0.5;            
yline(tR_fixed,'b--',LineWidth=2);
delta2_fixed=1;
plot(delta2_fixed,tR_fixed,'kx','MarkerSize',10,LineWidth=2);


subplot(2,3,4);
hold on; box on;
N_sites = 100;                   

delta2_scan = linspace(delta2_min, delta2_max, 800); 
nD = numel(delta2_scan);
 


H0 = ham(tL,tR_fixed,alpha,beta,N_sites,delta1,delta2_scan(1));
nStates = size(H0,1);
spec_abs = NaN(nStates, nD);
spec_real = NaN(nStates, nD); 
spec_imag = NaN(nStates, nD);
Eedge = NaN(1,nD);
cheb_N=100;
for ii = 1:nD
    d2 = delta2_scan(ii);
    H = ham(tL,tR_fixed,alpha,beta,N_sites,delta1,d2);
    E = eig(H);                    
   
    [~, ord] = sort(real(E));
    E = E(ord);
    
    if numel(E) ~= nStates
        Etmp = NaN(nStates,1);
        Etmp(1:min(numel(E),nStates)) = E(1:min(numel(E),nStates));
        E = Etmp;
    end
    spec_abs(:,ii)  = abs(E);
    spec_real(:,ii) = real(E);
    spec_imag(:,ii) = imag(E);
 
  C0 =-2 * tL * tR_fixed * cos(alpha - beta) + 2 * d2^2+delta1^2 ;       
  D0=2*sqrt((tR_fixed^2 -  d2^2)*(tL^2 - d2^2));
  mu0= tL^2*cos(alpha)^2/(tL^2 - d2^2) + tR_fixed^2*cos(beta)^2/(tR_fixed^2 - d2^2) - 1;
  kappa0= (cos(alpha)*tL/(tL^2 - d2^2) + cos(beta)*tR_fixed/(tR_fixed^2 - d2^2))*d2*delta1;
  E0 = -delta1*tL*tR_fixed *sin(alpha - beta)/(d2*(sin(alpha)*tL - sin(beta)*tR_fixed ));
    U0 = 4*(E0*cos(beta)*tR_fixed + d2*delta1)*(E0*tL*cos(alpha) + d2*delta1)/D0^2 - 1;
    eta0 = d2^2*delta1^2/(tR_fixed^2 - d2^2) + d2^2*delta1^2/(tL^2 - d2^2);
    GAMMA0 = E0^2*mu0/D0 + 2*E0*kappa0/D0 + eta0/D0 + C0/D0;
    Xe0 = (-2*d2^2 + 2*cos(alpha - beta)*tL*tR_fixed)/D0;
    Delta0 = (U0*Xe0 - GAMMA0)^2 - 4*Xe0^3*GAMMA0;
    Y0 = 1/2*(U0*Xe0 - GAMMA0 + sqrt(Delta0))/Xe0^2;
    Z0 = 1/2*(U0*Xe0 - GAMMA0 - sqrt(Delta0))/Xe0^2;


    if max(abs(chebyshevU(cheb_N,Y0)),abs(chebyshevU(cheb_N,Z0))) < abs(chebyshevU(cheb_N,Xe0))
        Eedge(ii) = E0;
    else
        Eedge(ii)=NaN;
    end

   
end

for k = 1:nStates
    plot(delta2_scan, spec_abs(k,:), 'k.', MarkerSize=3);
end
  % plot(delta2_scan, abs(Eedge), 'r.', MarkerSize=5);
xlabel('$\delta_2$','Interpreter','latex', 'FontSize', 18); ylabel('$|E|$ ','Interpreter','latex', 'FontSize', 18);
xlim([delta2_scan(1), delta2_scan(end)]);
 ylim([-1,7]);
grid on;
title(sprintf('|E|  vs \\delta_2 (t_R = %.2f)', tR_fixed));




subplot(2,3,2);hold on; box on;grid on;
% clearvars;
syms E A B Gamma  % 定义符号变量

N = 21; 

tR =0;    delta2=0.9;
G=25;
kai=4*cos(alpha)*cos(beta)*tL*tR;
mu=(tL^2*cos(alpha)^2)/(tL^2-delta2^2)+(tR^2*cos(beta)^2)/(tR^2-delta2^2)-1
 r=((tR^2-delta2^2)/(tL^2-delta2^2))^(1/4);

C=(-2*tL*tR*cos(alpha-beta)+2*delta2^2+delta1^2);
D=2*sqrt((tR^2-delta2^2)*(tL^2-delta2^2));

B2=(4*cos(alpha)*cos(beta)*tL*tR)/(D^2);
B1=(4*delta2*delta1*(tL*cos(alpha)+tR*cos(beta)))/(D^2);
B0=-1+(4*delta1^2*delta2^2)/D^2;


Gamma2=((tL^2*(cos(alpha))^2)/(tL^2-delta2^2)+(tR^2*(cos(beta))^2)/(tR^2-delta2^2)-1)/D;
Gamma1=(2*delta2*delta1*((tL*cos(alpha))/(tL^2-delta2^2)+(tR*cos(beta))/(tR^2-delta2^2)))/D;
Gamma0=(1/(tL^2-delta2^2)+1/(tR^2-delta2^2))*delta2^2*delta1^2/D+C/D;



% ---- （E）----
A =vpa(-(-E^2 + C)/D) ;
B =vpa(B2*E^2+B1*E+B0) ;
Gamma = vpa(Gamma2*E^2+Gamma1*E+Gamma0);

tic
% ---- 初始化符号向量 ----
lambda = sym(zeros(1, N+2));
mu = sym(zeros(1, N+2));
nu = sym(zeros(1, N+2));

% 初始条件
lambda(1) = 0; mu(1) = 0; nu(1) = 1;
lambda(2) = 0; mu(2) = 1; nu(2) = 0;

% ---- 递推符号 ----
for n = 2:N+1
    lambda(n+1) = vpa(simplify( 2*A*lambda(n)+ 2*mu(n) - lambda(n-1) ));
    mu(n+1)     = vpa(simplify( -2*B*lambda(n) + 2*nu(n) - mu(n-1) ));
    nu(n+1)     = vpa(simplify( 2*Gamma*lambda(n) - nu(n-1)) );
end
% 
MuE = vpa(simplify(-lambda(N+2)));



% 求Nu
lambda = sym(zeros(1, N+2));
mu = sym(zeros(1, N+2));
nu = sym(zeros(1, N+2));

lambda(1) = 0; mu(1) = 1; nu(1) = 0;
lambda(2) = 1; mu(2) = 0; nu(2) = 0;

% ---- 递推符号----
for n = 2:N+1
    lambda(n+1) = vpa(simplify(2*A*lambda(n)+ 2*mu(n) - lambda(n-1)));
    mu(n+1)     = vpa(simplify(-2*B*lambda(n) + 2*nu(n) - mu(n-1)) );
    nu(n+1)     = vpa(simplify(2*Gamma*lambda(n)- nu(n-1)));
end

NuE = vpa(simplify(-lambda(N+2)));

%
K=tL*exp(1i*alpha)- tR*exp(1i*beta);
P1=-D*(exp(1i*(alpha+beta))*tL*tR*delta1^2+K^2*delta2^2);

P2=((tR^2-delta2^2)*tL^2*exp(1i*alpha)^2+(tL^2-delta2^2)*tR^2*exp(1i*beta)^2)*delta1^2+(2*delta1^2-C)*delta2^2*K^2;

OBCE=vpa(MuE*P2+NuE*P1);

Etime=toc

 E_exac = vpa(root(expand(OBCE),E));   

 [E_exact, idx] = sort_complex_group(E_exac, 1e-10);

H = ham(tL,tR,alpha, beta,N,delta1,delta2);
[V_num,DD]=eig(H);
Enum=diag(DD);

[E_num, idx] = sort_complex_group(Enum, 1e-10);
xlabel('Re(E)',FontSize=18);
ylabel('Im(E)',FontSize=18);
set (gcf, 'color', 'white');



Ee_edge=-delta1*tL*tR*sin(alpha - beta)/(delta2*(sin(alpha)*tL - sin(beta)*tR));
% plot(real(Ee_edge), imag(Ee_edge), 'r*', 'MarkerSize', 8,'LineWidth',1);

Hk=zeros(2,2);
Ek=[];
K = linspace(0,2*pi,500);
for k=K 
      Hk(1,1)=tL*exp(1i*k)*exp(1i*alpha)+tR*exp(-1i*k)*exp(1i*beta);
      Hk(2,2)=tL*exp(1i*k)*exp(-1i*alpha)+tR*exp(-1i*k)*exp(-1i*beta);
      Hk(1,2)=delta2*(exp(1i*k)+exp(-1i*k))+delta1;
      Hk(2,1)=delta2*(exp(1i*k)+exp(-1i*k))+delta1;
      Ek0=eig(Hk);
      Ek=[Ek,Ek0];
end

exchange=Ek(1,1);
Ek(1,1)=Ek(2,1);
Ek(2,1)=exchange;
 exchange=Ek(1,251:500);
Ek(1,251:500)=Ek(2,251:500);
Ek(2,251:500)=exchange;
m = linspace(0.8,1,size(Ek(1,:),2));
patch(real(Ek(1,:)),imag(Ek(1,:)), m, 'edgecolor','flat','facecolor','none');hold on
patch(real(Ek(2,:)),imag(Ek(2,:)), m, 'edgecolor','flat','facecolor','none');



 z1 =sym(zeros(2*N, 1));
 z2=sym(zeros(2*N, 1));
 z3 =sym(zeros(2*N, 1));
 z4 =sym(zeros(2*N, 1));
for j=1:length(E_exact)
    E=E_exact(j);
b=((-2*E*tL*cos(alpha) - 2*delta2*delta1)/(tL^2 - delta2^2));
c=((2*tL*tR*cos(alpha - beta) + E^2 - delta1^2 - 2*delta2^2)/(tL^2 - delta2^2));
d=((-2*E*cos(beta)*tR - 2*delta2*delta1)/(tL^2 - delta2^2));
e=((tR^2 - delta2^2)/(tL^2 - delta2^2));

TT=[1, b, c, d, e];
           z=roots(TT);
           EQ1=vpa(z(1) + z(2) + z(3) + z(4)-((2*E*tL*cos(alpha) + 2*delta2*delta1)/(tL^2 -delta2^2)))
           EQ2=vpa(z(1)*z(2)*z(3) + z(1)*z(2)*z(4) + z(1)*z(3)*z(4) + z(2)*z(3)*z(4) -((2*E*cos(beta)*tR + 2*delta2*delta1)/(tL^2 - delta2^2)))
           z1(j)=z(1);
           z2(j)=z(2);
           z3(j)=z(3);
           z4(j)=z(4);

end

for ii=1:N
    CC=[z1(ii),z2(ii),z3(ii),z4(ii)];
    [res,ind]=sort(vpa(abs(CC)));
    z1(ii)=CC(ind(1));
    z2(ii)=CC(ind(2));
    z3(ii)=CC(ind(3));
    z4(ii)=CC(ind(4));
end

Gp= @(z)tL*exp(1i*alpha).*z + tR*exp(1i*beta)./z;
Gs= @(z)tL*exp(-1i*alpha).*z + tR*exp(-1i*beta)./z;
delta= @(z)delta2.*z +delta2./z + +delta1;
phi1 = delta(z1)./(E_exact - Gp(z1));
phi2 = delta(z2)./(E_exact - Gp(z2));
phi3 = delta(z3)./(E_exact - Gp(z3));
phi4 = delta(z4)./(E_exact - Gp(z4));
phi=@(z)delta(z)./(E_exact - Gp(z));
c4=(-phi(z1).*(phi(z2) - phi(z3)).*z1.^(N + 1) + phi(z2).*(phi(z1) - phi(z3)).*z2.^(N + 1)- phi(z3).*(phi(z1) - phi(z2)).*z3.^(N + 1))./(phi(z3).*phi(z4).*(phi(z1) - phi(z2)));
c3=(phi(z1).*(phi(z2) - phi(z4)).*z1.^(N + 1) - phi(z2).*(phi(z1) - phi(z4)).*z2.^(N + 1)+ phi(z4).*(phi(z1) - phi(z2)).*z4.^(N + 1))./(phi(z3).*phi(z4).*(phi(z1) - phi(z2)));
c2=(-phi(z1).*(phi(z3) - phi(z4)).*z1.^(N + 1) + phi(z3).*(phi(z1) - phi(z4)).*z3.^(N + 1)- phi(z4).*(phi(z1) - phi(z3)).*z4.^(N + 1))./(phi(z3).*phi(z4).*(phi(z1) - phi(z2)));
c1=(phi(z2).*(phi(z3) - phi(z4)).*z2.^(N + 1) - phi(z3).*(phi(z2) - phi(z4)).*z3.^(N + 1)+ phi(z4).*(phi(z2) - phi(z3)).*z4.^(N + 1))./(phi(z3).*phi(z4).*(phi(z1) - phi(z2)));

phiA=@(n) phi1.* z1.^n.*c1+phi2.* z2.^n.*c2+phi3.* z3.^n.*c3+phi4.* z4.^n.*c4;
phiB=@(n) z1.^n.*c1+ z2.^n.*c2+ z3.^n.*c3 + z4.^n.*c4;

V_exac=zeros(2*N,2*N);
for j = 2: 2 : 2*N
    V_exac(j-1,:)=phiA(j/2);
    V_exac(j,:)  =phiB(j/2);
end
%归一化
TTT=(sqrt(sum(abs(V_exac).^2)));
for ii=1:2*N
    V_exac(:,ii)=V_exac(:,ii)/TTT(ii);
end



for j = 1:2*N
    psi = V_exac(:, j); 
    

    left_weight = norm(psi(1:N))^2; 
    right_weight = norm(psi(N+1:end))^2; 
   

    if left_weight > 0.7 * (left_weight + right_weight) % 左局域
        color = [0.93,0.69,0.13];
    elseif right_weight > 0.7 * (left_weight + right_weight) % 右局域
        color = [0.47,0.67,0.19];
    else % 扩展态或双边局域
        color = 'k';
    end
    
   
    plot(real(E_exact(j)), imag(E_exact(j)), '.', 'Color', color, 'MarkerSize', 10);
 
    colors{j} = color; 

 
end
% ---- IPR ----
IPR = sum(abs(V_exac).^4, 1) ./ (sum(abs(V_exac).^2, 1).^2);


IPR_min = min(IPR);
IPR_max = max(IPR);

marker_size = 8 + 100 * (IPR - IPR_min) ./ (IPR_max - IPR_min + eps);

scatter(real(E_exact), imag(E_exact),  marker_size, [0.93,0.69,0.13], 'filled');

subplot(2,3,3);hold on; box on;grid on;



for j = 1:2*N
    
   plot(1:2*N, abs(V_exac(:, j)),'Color', colors{j},LineWidth=1);
end
xlabel('Site number');  ylabel('eigenvector number'); zlabel('Intensity');
xlim([1 2*N]);
xlabel("\sl n"); ylabel("|\psi|^2");
set(gca,'FontName','Times New Roman','FontSize', 12 ,'LineWidth', 1);


subplot(2,3,5);hold on; box on;grid on;


%
M_full    =2000;
theta_all = (2*pi/M_full) * (1:(M_full-1));  
theta_all(M_full/2) = [];            
M         = numel(theta_all);

for j = 1:M
    theta = theta_all(j);
   
 a0 = (exp(-1i*theta) - 1)^2 * (4*cos(beta)^2*exp(-1i*theta) *tR^2 - (tR^2 - delta2^2)*(exp(-1i*theta) + 1)^2) * (tR^2 - delta2^2);

a1 = 4 * delta2 * delta1 * (tR^2 - delta2^2) * (exp(-1i*theta) - 1) * (exp(-2*1i*theta) - 1);  

a2 = 4*cos(alpha)*cos(beta) * tL * tR * ( exp(1i*theta) * (exp(-2*1i*theta) + 1) * (exp(-1i*theta) - 1)^2 ) * (tR^2 - delta2^2) ...
    + 4*cos(beta)^2 * (exp(-1i*theta) - 1)^2 * (delta1^2 + 2*delta2^2 - 2* tL * tR* cos(alpha - beta))* tR^2 ...
    - 4*(exp(-1i*theta) - 1)^2 * delta2^2 * delta1^2;

a3 = 4*delta2 * delta1 * (2*1i*sin(theta) + exp(-2*1i*theta) - 1) * (-2*cos(alpha )*cos(beta) * tL * tR + cos(2*beta) * tR^2 + delta2^2);

a4 = 16*(cos(beta)*cos(alpha) * tL * tR * (2*tL* tR*cos(alpha - beta) - delta1^2 - 2*delta2^2) + delta1^2*delta2^2)*(cos(theta) - 1) ...
    - 8*(tL^2*(tR^2 - delta2^2)*cos(alpha)^2 + cos(beta)^2 * tR^2*(tL - delta2) * (tL + delta2))*(2*cos(theta) +1) *(cos(theta) - 1) ...
    -8*sin(theta)^2*(tR^2-delta2^2)*(tL^2-delta2^2);

a5 = 4*delta2 * delta1 *  (2*1i*sin(theta) - exp(2*1i*theta) + 1) * (2*cos(alpha )*cos(beta) * tL * tR - cos(2*alpha) * tL^2 - delta2^2);

a6 = 4*cos(alpha)*cos(beta) * tL * tR* ( exp(-1i*theta) * (exp(2*1i*theta) + 1) * (exp(1i*theta) - 1)^2 )* (tL^2 - delta2^2)  ...
    + 4*cos(alpha)^2*(exp(1i*theta) - 1)^2 * (delta1^2 + 2*delta2^2 - 2* tL * tR* cos(alpha - beta))* tL^2 ...
    - 4*(exp(1i*theta) - 1)^2 * delta2^2 * delta1^2;

a7 = 4 * delta2 * delta1 * (tL^2 - delta2^2)*(exp(1i*theta) - 1) * (exp(2*1i*theta) - 1);

a8 = (exp(1i*theta) - 1)^2 * (4*cos(alpha)^2*exp(1i*theta) *tL^2 - (tL^2 - delta2^2)*(exp(1i*theta) + 1)^2) * (tL^2 - delta2^2);

    p = [a8, a7, a6, a5, a4, a3, a2, a1, a0];

    z_roots = roots(p);    
   [z_sort,idx]=sort(z_roots,'ComparisonMethod','abs');

     z_all(:,j)=z_sort;
   
end
%% 画图：8 条分支 + 单位圆
theta1 = linspace(0, 2*pi, 400);
unit_x = cos(theta1);
unit_y = sin(theta1);

% plot(real(z_branch(:)), imag(z_branch(:)), 'bo', 'MarkerSize', 2);
zk1 = z_all(1, :);
zk2 = z_all(2, :);
zk3 = z_all(3, :);
zk4 = z_all(4, :);
zk5 = z_all(5, :);
zk6 = z_all(6, :);
zk7 = z_all(7, :);
zk8 = z_all(8, :);

 plot(real(zk1), imag(zk1), 'r.', MarkerSize=5);
  plot(real(zk2), imag(zk2), 'b.', MarkerSize=5); 
 plot(real(zk3), imag(zk3), 'g.', MarkerSize=5);
  plot(real(zk4), imag(zk4), 'm.', MarkerSize=8);
plot(real(zk5), imag(zk5), 'c.', MarkerSize=5);
plot(real(zk6), imag(zk6), '.', 'Color',[0.96,0.47,0.16],MarkerSize=5);
plot(real(zk7), imag(zk7), 'k.', MarkerSize=5);
plot(real(zk8), imag(zk8), '.', 'Color',[0.52,0.09,0.82],MarkerSize=5);
plot(unit_x, unit_y, 'k--', 'LineWidth', 1.2);
 % legend('z1','z2','z3','z4','z5','z6','z7','z8','BZ',fontsize=20);
% axis equal;
xlabel('Re(\beta)'); ylabel('Im(\beta)');
xlim([-1.5,1.5]);
ylim([-1.5,1.5]);
grid on;
%%numerical
N = 60; 
H = ham(tL,tR,alpha, beta,N,delta1,delta2);
[V_num,DD]=eig(H);
Enum=diag(DD);
e0=(tR^2 - delta2^2)/(tL^2 - delta2^2);
for j=1:length(Enum)
           EE=Enum(j);
           b0=(-2*EE*tL*cos(alpha) - 2*delta2*delta1)/(tL^2 - delta2^2);
           c0=(2*cos(alpha - beta)*tL*tR + EE^2 - delta1^2 - 2*delta2^2)/(tL^2 - delta2^2);
           d0=(-2*EE*cos(beta)*tR - 2*delta2*delta1)/(tL^2 - delta2^2);
           
           pnum=[1, b0, c0, d0, e0];
             z=sort(roots(pnum),'ComparisonMethod','abs');
          
          
           z11(j)=z(1);
           z22(j)=z(2);
           z33(j)=z(3);
           z44(j)=z(4);
       
     end
      % plot(real(z11),imag(z11),'ro',MarkerSize=10);
        plot(real(z22),imag(z22),'ko',MarkerSize=5,LineWidth=1);
        plot(real(z33),imag(z33),'ro',MarkerSize=5,LineWidth=1);
        % plot(real(z44),imag(z44),'ko',MarkerSize=5);
function H = ham(tL,tR,alpha, beta,N,delta1,delta2)

a1=zeros(1,N-1);
a1(1:1:end)=tL*exp(1i*alpha);

a2=zeros(1,N-1);
a2(1:1:end)=tR*exp(1i*beta);

H1=zeros(N,N);
H1=diag(a1,1)+diag(a2,-1);

a3=zeros(1,N-1);
a3(1:1:end)=tL*exp(-1i*alpha);

a4=zeros(1,N-1);
a4(1:1:end)=tR*exp(-1i*beta);

H2=zeros(N,N);
H2=diag(a3,1)+diag(a4,-1);



a5=zeros(1,N);
a5(1:1:end)=delta1;

a6=zeros(1,N-1);
a6(1:1:end)=delta2;

H3=zeros(N,N);
H3=diag(a5)+diag(a6,-1)+diag(a6,1);

H_ori=zeros(2*N,2*N);
H_ori=[H1,H3;H3,H2];

 % 创建重排索引：交替选择基矢
    reorder_idx = zeros(1, 2*N);
    reorder_idx(1:2:end) = 1:N;           % 奇数位置放第一条链
    reorder_idx(2:2:end) = (1:N) + N;     % 偶数位置放第二条链
    
    % 重排哈密顿量
    H=zeros(2*N,2*N);
    H = H_ori(reorder_idx, reorder_idx);
end

function [Theta, Psi, Phi] = genPoly(G,tL,tR,alpha, beta,delta1,delta2,E_single,X,Y,Z)
        
        digits(G)
        eps1=10^(-10);
        

        r=((tR^2-delta2^2)/(tL^2-delta2^2))^(1/4);

        Theta= vpa(acos(X));
        Phi= vpa(acos(Y));
        Psi= vpa(acos(Z));
    
      
            EQ1=vpa(cos((1/2)*Phi)*cos((1/2)*Psi)*cos((1/2)*Theta)-((E_single*(tL*cos(alpha)*r^2+cos(beta)*tR)+delta1*delta2*(r^2+1))/(4*r^3*(tL^2-delta2^2))),20);
            EQ2=vpa(sin((1/2)*Phi)*sin((1/2)*Psi)*sin((1/2)*Theta)-(1i*(E_single*(-tL*cos(alpha)*r^2+cos(beta)*tR)+delta1*delta2*(-r^2+1))/(4*r^3*(tL^2-delta2^2))),20); 
            
            if  vpa(abs(EQ1))<eps1 && vpa(abs(EQ2))>eps1
                Theta = -Theta;
            elseif vpa(abs(EQ1))>eps1 && vpa(abs(EQ2))<eps1
                Theta = vpa(2*pi-Theta);
            elseif vpa(abs(EQ1))>eps1 && vpa(abs(EQ2))>eps1
                Theta = vpa(2*pi+Theta);
            end
            EQ1=vpa(cos((1/2)*Phi).*cos((1/2)*Psi).*cos((1/2)*Theta)-((E_single*(tL*cos(alpha)*r^2+cos(beta)*tR)+delta1*delta2*(r^2+1))/(4*r^3*(tL^2-delta2^2))),20)
            EQ2=vpa(sin((1/2)*Phi).*sin((1/2)*Psi).*sin((1/2)*Theta)-(1i*(E_single*(-tL*cos(alpha)*r^2+cos(beta)*tR)+delta1*delta2*(-r^2+1))/(4*r^3*(tL^2-delta2^2))),20)
      
       
    end
    



function [E_sorted, idx] = sort_complex_group(E, tol)

    if nargin < 2
        tol = 1e-10; % 默认容差
    end

    % 统一为列向量
    E = E(:);

    % ********** 第一步：按实部排序 **********
    [~, idx] = sort(real(E), 'ascend');
    E_sorted = E(idx);

    re = real(E_sorted);
    im = imag(E_sorted);

    N = numel(E_sorted);
    i = 1;

    % ********** 第二步：按实部近似分段 + 段内按虚部排序 **********
    while i < N
        j = i;

        % 找出 j，使 re(i:j) 在 tol 内
        while j < N && abs(re(j+1) - re(i)) <= tol
            j = j + 1;
        end

        
        if j > i
            segment = i:j;
            [~, local_order] = sort(im(segment), 'ascend');
         
            segment_sorted = segment(local_order);

            
            E_sorted(i:j) = E_sorted(segment_sorted);
            idx(i:j)      = idx(segment_sorted);
            im(i:j) = imag(E_sorted(i:j));
        end

        i = j + 1;
    end
end

