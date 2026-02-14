clc; clear; close all;

 tR = 1;
tL =2;
alpha=pi/4;
beta = pi/4;






subplot(1,2,1);hold on; box on;grid on;

N=30;
 delta1 =1;    delta2=2.5;
syms E A B Gamma  
 

G=30;
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

P1=-D;

P2=-2*delta2^2 + 2*cos(alpha - beta)*tL*tR;

OBCE=vpa(MuE*P2+NuE*P1);

Etime=toc

 E_exac = vpa(root(expand(OBCE),E));   

 [E_exact, idx] = sort_complex_group(E_exac, 1e-10);

%%微扰
H = ham(tL,tR,alpha, beta,N,delta1,delta2);
alpha_ptb=pi/3;
beta_ptb=pi/6;
H(N+1,N-1)=-tL*alpha_ptb;
H(N,N-2)=tL*alpha_ptb;
H(N-1,N+1)=-tR*beta_ptb;
H(N-2,N)=tR*beta_ptb;

[V_num,DD]=eig(H);
Enum=diag(DD);


[E_num, idx] = sort_complex_group(Enum, 1e-10);
 V_num=V_num(:,idx); 



plot(real(E_num), imag(E_num), 'bo', 'MarkerSize', 8);

xlabel('Re(E)',FontSize=18);
ylabel('Im(E)',FontSize=18);
set (gcf, 'color', 'white');





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
m=linspace(0.8,1, size(Ek(1,:),2));
colormap(cool); 
 patch(real(Ek(1,:)),imag(Ek(1,:)), m, 'edgecolor','flat','facecolor','none');hold on
patch(real(Ek(2,:)),imag(Ek(2,:)), m, 'edgecolor','flat','facecolor','none');



%------eigenstates------e



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


colors = zeros(2*N, 3);
for j = 1:2*N
    psi = V_exac(:, j); 
    

    left_weight = norm(psi(1:N))^2; 
    right_weight = norm(psi(N+1:end))^2; 
   

    if left_weight > 0.6* (left_weight + right_weight) % 左局域
        color = [0.93,0.69,0.13];

    elseif right_weight > 0.6* (left_weight + right_weight) % 右局域
        color = [0.47,0.67,0.19];
    else % 扩展态或双边局域
        color = [0.15,0.15,0.15];
    end
    
   
      plot(real(E_exact(j)), imag(E_exact(j)), '.', 'Color', color, 'MarkerSize', 12);
 
  colors(j, :) = color;
 
end
Ee_edge=-delta1*tL*tR*sin(alpha - beta)/(delta2*(sin(alpha)*tL - sin(beta)*tR));
 plot(real(Ee_edge), imag(Ee_edge), 'r*', 'MarkerSize', 8,'LineWidth',1);

subplot(1,2,2);hold on; box on;grid on;
tolerance = 1e-2;
index=find(abs(E_exact) < tolerance);
N1=index(1);N2=index(2); %零模位置
h(1)=mesh(1:2*N,[N1,N2],abs(V_exac(:,[N1,N2])').^2);
set(h(1), 'FaceColor','interp', 'EdgeColor','interp', 'FaceLighting','phong','MarkerFaceColor',[1 0 0],...
 'MarkerSize',3,'Marker','o','LineStyle','none','FaceColor','none','EdgeColor',[1 0 0]);

h(2)=mesh(1:2*N,1:2*N,abs(V_num(:,1:2*N)').^2);
colormap(parula); 
 set(h(2), 'FaceColor','interp', 'EdgeColor','interp', 'FaceLighting','none'); hold on

view(3), camlight left, axis tight

xlim([1 2*N]);
xlabel("\sl n"); ylabel("Serial number"); zlabel("\sl |\langle n|\psi\rangle|^2");
set(gca,'FontName','Times New Roman','FontSize', 12 ,'LineWidth', 1);





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

        % 如果该段长度 > 1，则按虚部排序
        if j > i
            segment = i:j;
            [~, local_order] = sort(im(segment), 'ascend');

            % 排序这段的位置
            segment_sorted = segment(local_order);

            % 应用到 E_sorted 和 idx
            E_sorted(i:j) = E_sorted(segment_sorted);
            idx(i:j)      = idx(segment_sorted);

            % 更新 imag（防止下一轮使用旧的）
            im(i:j) = imag(E_sorted(i:j));
        end

        i = j + 1;
    end
end

