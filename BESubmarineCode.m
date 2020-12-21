% State representation Submarine 30 knots 
close all; clear all; clc;
a11=-0.19003;a12=4.4802;a14=0.0014673;
a21=0.0085526;a22=-0.45988;a24=-0.0056095;
a34=-15.433;

b11=-0.1855;b12=-0.57149;
b21=0.043308;b22=-0.055543;

A=[a11 a12 0 a14; a21 a22 0 a24;1 0 0 a34;0 1 0 0];
B=[b11 b12;b21 b22; 0 0; 0 0];
C=[0 0 1 0; 0 0 0 1];

sys = ss(A,B,C,zeros(2,2));
%% Analyse of the natural system
% Check observability and contorlability
Ob=obsv(A,C);
unob=length(A)-rank(Ob);%Observable

Co=ctrb(A,B);
unco=length(A)-rank(Co);%Controllable

sys = ss(A,B,C,zeros(2,2));
damp(sys);
figure 
step(sys); grid on;
figure 
impulse(sys); grid on;
figure 
pzmap(sys); grid on;
%% State feedback control 
lambda1 =  -3/100;%Pole to achieve a rise time = 100
lambda2 =  -3.85/100;
lambda3 = -3/80;%Pole to achieve a rise time = 80
lambda4 =  -3.85/80;
lambda = [lambda1 lambda2 lambda3 lambda4];

K = place(A,B,lambda)
H = -1*(C*(A-B*K)^(-1)*B)^(-1)

sys2 = ss(A-B*K,B*H,C,zeros(2,2));
G2 = tf(sys2);
figure
step(sys2); grid on;

%% Insenstive Mode placement
% For the first pole placement with
lambda1 = [-3/80 -3.85/80 -3/100 -3.85/100];

K = place(A,B,lambda1)
[AcEigVector,AEigValues]=eig(A-B*K);
Cond1=cond(AcEigVector)
% For the second pole placement with

lambda2 = [-3/80 -15/80 -3/100 -15/100 ];

K2 = place(A,B,lambda2)
[AcEigVector2,AEigValues2]=eig(A-B*K2);
Cond2=cond(AcEigVector2)

% For the third pole placement with

lambda3 = [-3/80 -15 -3/100 -15];

K3 = place(A,B,lambda3)
[AcEigVector3,AEigValues3]=eig(A-B*K3);
Cond3=cond(AcEigVector3)

% For the fourth pole placement with
lambda4 = [-3/80 -50 -3/100 -50];

K4 = place(A,B,lambda4)
[AcEigVector4,AEigValues4]=eig(A-B*K);
Cond4=cond(AcEigVector4)
% Calculate the maximum admissible perturbation
Pmaxpermited=0;

% for i=0.001:0.001:0.01
%  [maxreel1,ecartmax1]=sensible(A-B*K4,i);
%  if (maxreel1<0 && Pmaxpermited<i)
%      Pmaxpermited=i;
%  end
% end

% Sensibility of the pole placement for several pmax 
figure 
[maxreel1,ecartmax1]=sensible(A-B*K,0.5)

figure 
[maxreel2,ecartmax2]=sensible(A-B*K,0.1)

figure 
[maxreel3,ecartmax3]=sensible(A-B*K,0.01)

figure 
[maxreel4,ecartmax4]=sensible(A-B*K,0.001)

%% Degraded mode
%For U1
AU1=A; BU1=B(:,1); CU1=C(2,:);

KU1 = place(AU1,BU1,lambda1);
%HU1 = -(((CU1*((AU1-(BU1*KU1))^(-1)))*BU1)^(-1));
sys3=ss(AU1-BU1*KU1,BU1,CU1,zeros(1,1));

figure
step(sys3); grid on;
[AcEigVectorU1,AEigValuesU1]=eig(AU1-BU1*KU1);
CondU1=cond(AcEigVectorU1)
% With the same perturbations for U1
figure 
[maxreelU11,ecartmaxU11]=sensible(AU1-BU1*KU1,0.5)

figure 
[maxreelU12,ecartmaxU12]=sensible(AU1-BU1*KU1,0.1)

figure 
[maxreelU13,ecartmaxU13]=sensible(AU1-BU1*KU1,0.01)

figure 
[maxreelU14,ecartmaxU14]=sensible(AU1-BU1*KU1,0.001)

PmaxpermitedU1=0;
%Algorithme for find the Pmax (Heavy comsomption)
% for i=0.00001:0.00001:0.001
%  [maxreel1,ecartmax1]=sensible(AU1-BU1*KU1,i);
%  if (maxreel1<0 && PmaxpermitedU1<i)
%      PmaxpermitedU1=i;
%  end
% end
%%
%For U2
AU2=A; BU2=B(:,2); CU2=C(1,:);

KU2 = place(AU2,BU2,lambda1);
HU2 = -(((CU2*((AU2-(BU2*KU2))^(-1)))*BU2)^(-1));
sys4=ss(AU2-BU2*KU2,BU2*HU2,CU2,zeros(1,1));

figure
step(sys4); grid on;
[AcEigVectorU2,AEigValuesU2]=eig(AU2-BU2*KU2);
CondU2=cond(AcEigVectorU2)
% With the same perturbations for U2
figure 
[maxreelU21,ecartmaxU21]=sensible(AU2-BU2*KU2,0.5)

figure 
[maxreelU22,ecartmaxU22]=sensible(AU2-BU2*KU2,0.1)

figure 
[maxreelU23,ecartmaxU23]=sensible(AU2-BU2*KU2,0.01)

figure 
[maxreelU24,ecartmaxU24]=sensible(AU2-BU2*KU2,0.001)

PmaxpermitedU2=0;

% for i=0.00001:0.00001:0.001
%  [maxreel1,ecartmax1]=sensible(AU1-BU1*KU1,i);
%  if (maxreel1<0 && PmaxpermitedU2<i)
%      PmaxpermitedU2=i;
%  end
% end

%% Eigenstructure assignment and decoupling
% We choose the modes for the closed-loop
lambda = [ -3/80 -3.85/80 -3/100 -3.85/100];
%We obtain the sub-space for each pole
N1=-(inv(eye(4,4)*lambda(1,1)-A))*B
N2=-(inv(eye(4,4)*lambda(1,2)-A))*B
N3=-(inv(eye(4,4)*lambda(1,3)-A))*B
N4=-(inv(eye(4,4)*lambda(1,4)-A))*B
%From the decoupling constraints, we know the form of the matrix
%CV=[x x 0 0
%    0 0 x x]
%We obtain the next condition in the eigenvectors
%v14=0, v24=0, v33=0, v43=0;
%Then our matrix have the following shape
%V=[a d g j;
%   b e h k
%   c f 0 0
%   0 0 i l];
%Choosing i=j=k=l=2 we can calculate the others values
c=2;f=2;i=2;l=2;
% V=[a d g j;
%    b e h k
%    c f 0 0
%    0 0 i l];
%If we choose M= I (nxn) --> w=z
w1=inv(N1(3:4,:))*[c;0];
v1=N1*w1;

w2=inv(N2(3:4,:))*[f;0];
v2=N2*w2;

w3=inv(N3(3:4,:))*[0;i];
v3=N3*w3;

w4=inv(N4(3:4,:))*[0;l];
v4=N4*w4;

V=[v1 v2 v3 v4];
%K=W*U
K=[w1 w2 w3 w4]*inv(V)
%We obtain the matrix H from the relation UBH=Bdc, knowing that the form
%of the B in closed loop is:
%Bdc=[x 0;x 0;0 x; 0 x];
%We have to choose two values and calculate the other values
%Then we choose
%Bdc=[1 0;x 0;0 1; 0 x];
%From UBH=Bdc
UB=inv(V)*B;
UB_=[UB(1,:);UB(3,:)];
H=inv(UB_);
%We rebuild the closed loop matrix

%X*=U(A-BK)V + UBHYc
%Y=CV
Ac=inv(V)*(A-B*K)*V;
Bc=inv(V)*B*H
Cc=C*V;

SystemCL=ss(Ac,Bc,Cc,zeros(2,2));
damp(SystemCL);
figure 
step(SystemCL); grid on;

% Study of the sensitivity 
[AcEigVector,AEigValues]=eig(Ac);
Cond=cond(AcEigVector)

Pmaxpermited=0;

% for i=0.001:0.0001:0.1
%  [maxreel1,ecartmax1]=sensible(Ac,i);
%  if (maxreel1<0 && Pmaxpermited<i)
%      Pmaxpermited=i;
%  end
% end

figure 
[maxreel,ecartmax]=sensible(Ac,0.001)

%% Decouplage in static mode
Omega=  inv([A , B;
    C , zeros(2,2)]);

Om11=Omega(1:4,1:4);
Om12=Omega(1:4,5:6);
Om21=Omega(5:6,1:4);
Om22=Omega(5:6,5:6);

H_2=Om22+K*Om12;

sys_h2 = ss(A-B*K, B*H_2, C, zeros(2,2));
figure(3);
step(sys_h2);

figure(4);
impulse(sys_h2);

