close all; clearvars;
tL =2; alphaR= pi/4; alphaL=pi/3; delta1=1;
tR =0.6;    delta2=0.4;
N =80;


Hk=zeros(2,2);
Ek=[];
K = linspace(0,2*pi,500);
for k=K 
      % D0=0;
       D0=tL*cos(alphaL)*exp(1i*k)+tR*cos(alphaR)*exp(-1i*k);
      Hk(1,1)=tL*exp(1i*k)*exp(1i*alphaL)+tR*exp(-1i*k)*exp(1i*alphaR)-D0;
      Hk(2,2)=tL*exp(1i*k)*exp(-1i*alphaL)+tR*exp(-1i*k)*exp(-1i*alphaR)-D0;
      Hk(1,2)=delta2*(exp(1i*k)+exp(-1i*k))+delta1;
      Hk(2,1)=delta2*(exp(1i*k)+exp(-1i*k))+delta1;
      Ek0=eig(Hk);
      Ek=[Ek,Ek0];
end
subplot(2,2,1);hold on; grid on; box on;

 m=linspace(0.8,1, size(Ek(1,:),2));

patch(real(Ek(1,:)),imag(Ek(1,:)), m, 'edgecolor','flat','facecolor','none');hold on
patch(real(Ek(2,:)),imag(Ek(2,:)), m, 'edgecolor','flat','facecolor','none');
 colormap(cool); 
xlabel('Re(E)',FontSize=16); ylabel('Im(E)',FontSize=16);
set (gcf, 'color', 'white');


%%numerical
 G=10;
HA = hamA(tL,tR,alphaL, alphaR,N,delta1,delta2);
[~,DD]=eig((HA));
Enum=diag(DD);
plot(real(Enum),imag(Enum),'ko','MarkerFaceColor','k','MarkerSize',3);hold on

for j=1:length(Enum)
           EE=(Enum(j));
           a0=-tL^2*cos(2*alphaL) - 2*delta2^2 + tL^2;
           b0=-4*delta1*delta2;
           c0=4*tL*tR*sin(alphaL)*sin(alphaR) + 2*EE^2 - 2*delta1^2 - 4*delta2^2;
           d0=-4*delta1*delta2;
           e0=-tR^2*cos(2*alphaR) - 2*delta2^2 + tR^2;
           pnum=[a0, b0, c0, d0, e0];
             z=sort(roots(pnum),'ComparisonMethod','abs');
           z11(j)=z(1);
           z22(j)=z(2);
           z33(j)=z(3);
           z44(j)=z(4);
       
end


subplot(2,2,2);hold on; grid on; box on;

M_full    =4000;
theta_all = (2*pi/M_full) * (1:(M_full-1));  
theta_all(M_full/2) = [];            
M         = numel(theta_all);

for j = 1:M
    theta = theta_all(j);
   
 a0 = (2*cos(alphaR)^2*tR^2 - 2*tR^2 + 2*delta2^2)*exp(2*1i*theta) - tR^2*cos(2*alphaR) + tR^2 - 2*delta2^2;

a1 = 4*exp(theta*1i)*delta1*delta2*(exp(theta*1i) - 1);  

a2 = 0;

a3 = -4*exp(2*1i*theta)*delta1*delta2*(exp(theta*1i) - 1);

a4 = -2*(cos(alphaL)^2*tL^2 + delta2^2 - tL^2)*exp((2*theta)*1i)*(exp((2*theta)*1i) - 1);



    p = [a4, a3, a2, a1, a0];

    z_roots = roots(p);    
   [z_sort,idx]=sort(z_roots,'ComparisonMethod','abs');
     z_all(:,j)=z_sort;
   
end
zk1 = z_all(1, :);
zk2 = z_all(2, :);
zk3 = z_all(3, :);
zk4 = z_all(4, :);
plot(real(zk2), imag(zk2), 'b.', MarkerSize=3); 
 plot(real(zk3), imag(zk3), 'b.',  MarkerSize=3);
  plot(real(zk4), imag(zk4), 'b.', MarkerSize=3);
 plot(real(zk1), imag(zk1), 'b.',  MarkerSize=3);

   
        plot(real(z22),imag(z22),'ro','MarkerFaceColor','r','MarkerSize',3);
        plot(real(z33),imag(z33),'ro','MarkerFaceColor','r','MarkerSize',3);
        theta1 = linspace(0, 2*pi, 400);
        unit_x = cos(theta1);
        unit_y = sin(theta1);
    plot(unit_x, unit_y, 'k--', 'LineWidth', 1.2);
xlabel('Re(\beta)',FontSize=16); ylabel('Im(\beta)',FontSize=16);
   xlim([-1.2,1.2]);
ylim([-1.2,1.2]);

%%%%% our model
clearvars;
tL =2; alphaR= pi/4; alphaL=pi/3; delta1=1;
tR =0.6;    delta2=0.4;
N =80;

Hk=zeros(2,2);
Ek=[];
K = linspace(0,2*pi,500);
for k=K 
      D0=0;
    
      Hk(1,1)=tL*exp(1i*k)*exp(1i*alphaL)+tR*exp(-1i*k)*exp(1i*alphaR)-D0;
      Hk(2,2)=tL*exp(1i*k)*exp(-1i*alphaL)+tR*exp(-1i*k)*exp(-1i*alphaR)-D0;
      Hk(1,2)=delta2*(exp(1i*k)+exp(-1i*k))+delta1;
      Hk(2,1)=delta2*(exp(1i*k)+exp(-1i*k))+delta1;
      Ek0=eig(Hk);
      Ek=[Ek,Ek0];
end
subplot(2,2,3);hold on; grid on; box on;

m=linspace(0.8,1, size(Ek(1,:),2));


patch(real(Ek(1,:)),imag(Ek(1,:)), m, 'edgecolor','flat','facecolor','none');hold on
patch(real(Ek(2,:)),imag(Ek(2,:)), m, 'edgecolor','flat','facecolor','none');
 colormap(cool); 
colormap(cool); 
xlabel('Re(E)',FontSize=16); ylabel('Im(E)',FontSize=16);

%%numerical
 
HB = hamB(tL,tR,alphaL, alphaR,N,delta1,delta2);
[~,DD]=eig((HB));
Enum=diag(DD);
plot(real(Enum),imag(Enum),'ko','MarkerFaceColor','k','MarkerSize',3);hold on
e0=(tR^2 - delta2^2)/(tL^2 - delta2^2);
Ee_edge=-delta1*tL*tR*sin(alphaL - alphaR)/(delta2*(sin(alphaL)*tL - sin(alphaR)*tR));
 plot(real(Ee_edge), imag(Ee_edge), 'r*', 'MarkerSize', 8,'LineWidth',1);

for j=1:length(Enum)
           EE=(Enum(j));
           b0=(-2*EE*tL*cos(alphaL) - 2*delta2*delta1)/(tL^2 - delta2^2);
           c0=(2*cos(alphaL - alphaR)*tL*tR + EE^2 - delta1^2 - 2*delta2^2)/(tL^2 - delta2^2);
           d0=(-2*EE*cos(alphaR)*tR - 2*delta2*delta1)/(tL^2 - delta2^2);
           
           pnum=[1, b0, c0, d0, e0];
             z=sort(roots(pnum),'ComparisonMethod','abs');
           z11(j)=z(1);
           z22(j)=z(2);
           z33(j)=z(3);
           z44(j)=z(4);
       
end
subplot(2,2,4);hold on; grid on; box on;
M_full    =4000;
theta_all = (2*pi/M_full) * (1:(M_full-1));  
theta_all(M_full/2) = [];            
M         = numel(theta_all);

for j = 1:M
    theta = theta_all(j);
   
 a0 = (exp(-1i*theta) - 1)^2 * (4*cos(alphaR)^2*exp(-1i*theta) *tR^2 - (tR^2 - delta2^2)*(exp(-1i*theta) + 1)^2) * (tR^2 - delta2^2);

a1 = 4 * delta2 * delta1 * (tR^2 - delta2^2) * (exp(-1i*theta) - 1) * (exp(-2*1i*theta) - 1);  

a2 = 4*cos(alphaL)*cos(alphaR) * tL * tR * ( exp(1i*theta) * (exp(-2*1i*theta) + 1) * (exp(-1i*theta) - 1)^2 ) * (tR^2 - delta2^2) ...
    + 4*cos(alphaR)^2 * (exp(-1i*theta) - 1)^2 * (delta1^2 + 2*delta2^2 - 2* tL * tR* cos(alphaL - alphaR))* tR^2 ...
    - 4*(exp(-1i*theta) - 1)^2 * delta2^2 * delta1^2;

a3 = 4*delta2 * delta1 * (2*1i*sin(theta) + exp(-2*1i*theta) - 1) * (-2*cos(alphaL )*cos(alphaR) * tL * tR + cos(2*alphaR) * tR^2 + delta2^2);

a4 = 16*(cos(alphaR)*cos(alphaL) * tL * tR * (2*tL* tR*cos(alphaL - alphaR) - delta1^2 - 2*delta2^2) + delta1^2*delta2^2)*(cos(theta) - 1) ...
    - 8*(tL^2*(tR^2 - delta2^2)*cos(alphaL)^2 + cos(alphaR)^2 * tR^2*(tL - delta2) * (tL + delta2))*(2*cos(theta) +1) *(cos(theta) - 1) ...
    -8*sin(theta)^2*(tR^2-delta2^2)*(tL^2-delta2^2);

a5 = 4*delta2 * delta1 *  (2*1i*sin(theta) - exp(2*1i*theta) + 1) * (2*cos(alphaL )*cos(alphaR) * tL * tR - cos(2*alphaL) * tL^2 - delta2^2);

a6 = 4*cos(alphaL)*cos(alphaR) * tL * tR* ( exp(-1i*theta) * (exp(2*1i*theta) + 1) * (exp(1i*theta) - 1)^2 )* (tL^2 - delta2^2)  ...
    + 4*cos(alphaL)^2*(exp(1i*theta) - 1)^2 * (delta1^2 + 2*delta2^2 - 2* tL * tR* cos(alphaL - alphaR))* tL^2 ...
    - 4*(exp(1i*theta) - 1)^2 * delta2^2 * delta1^2;

a7 = 4 * delta2 * delta1 * (tL^2 - delta2^2)*(exp(1i*theta) - 1) * (exp(2*1i*theta) - 1);

a8 = (exp(1i*theta) - 1)^2 * (4*cos(alphaL)^2*exp(1i*theta) *tL^2 - (tL^2 - delta2^2)*(exp(1i*theta) + 1)^2) * (tL^2 - delta2^2);

    p = [a8, a7, a6, a5, a4, a3, a2, a1, a0];

    z_roots = roots(p);    
   [z_sort,idx]=sort(z_roots,'ComparisonMethod','abs');
     z_all(:,j)=z_sort;
   
end


zk1 = z_all(1, :);
zk2 = z_all(2, :);
zk3 = z_all(3, :);
zk4 = z_all(4, :);
zk5 = z_all(5, :);
zk6 = z_all(6, :);
zk7 = z_all(7, :);
zk8 = z_all(8, :);

 plot(real(zk1), imag(zk1), 'b.',  MarkerSize=3);
   plot(real(zk2), imag(zk2), 'b.', MarkerSize=3); 
 plot(real(zk3), imag(zk3), 'b.',  MarkerSize=3);
  plot(real(zk4), imag(zk4), 'b.', MarkerSize=3);
plot(real(zk5), imag(zk5), 'b.',  MarkerSize=3);
plot(real(zk6), imag(zk6), 'b.', MarkerSize=3);
plot(real(zk7), imag(zk7),'b.',  MarkerSize=3);
 plot(real(zk8), imag(zk8),'b.',  MarkerSize=3);

        plot(real(z22),imag(z22),'ro','MarkerFaceColor','r','MarkerSize',3);
        plot(real(z33),imag(z33),'ro','MarkerFaceColor','r','MarkerSize',3);
        theta1 = linspace(0, 2*pi, 400);
        unit_x = cos(theta1);
        unit_y = sin(theta1);
    plot(unit_x, unit_y, 'k--', 'LineWidth', 1.2);
    xlabel('Re(\beta)',FontSize=16); ylabel('Im(\beta)',FontSize=16);
    xlim([-1.2,1.2]);
ylim([-1.2,1.2]);

EE=Ee_edge;
b0=(-2*EE*tL*cos(alphaL) - 2*delta2*delta1)/(tL^2 - delta2^2);
 c0=(2*cos(alphaL - alphaR)*tL*tR + EE^2 - delta1^2 - 2*delta2^2)/(tL^2 - delta2^2);
  d0=(-2*EE*cos(alphaR)*tR - 2*delta2*delta1)/(tL^2 - delta2^2);
           
pnum=[1, b0, c0, d0, e0];
 z=sort(roots(pnum),'ComparisonMethod','abs');
           
          
        plot(real(z(1)),imag(z(1)),'go', MarkerSize=5,LineWidth=1.5);
        plot(real(z(2)),imag(z(2)),'go', MarkerSize=5,LineWidth=1.5);
         plot(real(z(3)),imag(z(3)),'go', MarkerSize=5,LineWidth=1.5);
        plot(real(z(4)),imag(z(4)),'go', MarkerSize=5,LineWidth=1.5);


 function HB = hamB(tL,tR,alphaL, alphaR,N,delta1,delta2)

a1=zeros(1,N-1);
a1(1:1:end)=tL*exp(1i*alphaL);

a2=zeros(1,N-1);
a2(1:1:end)=tR*exp(1i*alphaR);

H1=zeros(N,N);
H1=diag(a1,1)+diag(a2,-1);

a3=zeros(1,N-1);
a3(1:1:end)=tL*exp(-1i*alphaL);

a4=zeros(1,N-1);
a4(1:1:end)=tR*exp(-1i*alphaR);

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


    reorder_idx = zeros(1, 2*N);
    reorder_idx(1:2:end) = 1:N;           % 奇数位置放第一条链
    reorder_idx(2:2:end) = (1:N) + N;     % 偶数位置放第二条链
    
 
    HB=zeros(2*N,2*N);
    HB = H_ori(reorder_idx, reorder_idx);
 end
  function HA = hamA(tL,tR,alphaL, alphaR,N,delta1,delta2)

a1=zeros(1,N-1);
a1(1:1:end)=1i*tL*sin(alphaL);

a2=zeros(1,N-1);
a2(1:1:end)=1i*tR*sin(alphaR);

H1=zeros(N,N);
H1=diag(a1,1)+diag(a2,-1);

a3=zeros(1,N-1);
a3(1:1:end)=-1i*tL*sin(alphaL);

a4=zeros(1,N-1);
a4(1:1:end)=-1i*tR*sin(alphaR);

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

 
    reorder_idx = zeros(1, 2*N);
    reorder_idx(1:2:end) = 1:N;           
    reorder_idx(2:2:end) = (1:N) + N;   
    
    HA=zeros(2*N,2*N);
    HA = H_ori(reorder_idx, reorder_idx);
end