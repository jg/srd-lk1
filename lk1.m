function lk1()
% Kamil Jackiewicz
% Krzysztof Koziol

global n m V0;
n = 2;              % rozmiar wektora stanu
m = 1;              % rozmiar wektora sterowania
V0 = [0; 0];   % Wartosc oczekiwana zaklocenia

xx0=[1;2;3];        % punkt pocz¹tkowy optymalizacji
xx = fmincon(@(x)fun1(x,n),xx0,[],[],[],[],[],[],@(x)con(x,n,V0))

% xx =
%     6.2839 - x1 stacjonarne
%    -1.2739 - x2 stacjonarne
%     1.3964 - sterowanie stacjonarne

xs = xx(1:n);
us = xx((n+1):(n+m));
disp('Sprawdzenie stacjonarnosci punktu:');
transf(xs, us, V0,n) - xs

% Wyniki:
%  1.0e-011 *
%         0
%   -0.7611

% Znalezienie modelu liniowego obiektu w uzyskanym punkcie pracy 
% oraz kwadratowego przyblizenia wskaznika jakosci
[A,B,G,C,R,r,Q,q,H] = model_liniowy(xs, us, V0);

% Sprowadzenie do postaci z rozdzialu 3 instrukcji:
% Znajdujemy podstawienie, w wyniku ktorego wskaznik jakosci bedzie 
% zawierac jedynie macierze Q i R
% [D,ud,xd,An,Bn,Gn,Cn,Qn,Rn] = model_liniowy_nowy(A,B,G,C,R,r,Q,q,H,n);

An = A;
Bn = B';
Cn = C;
Qn = Q;
Rn = R;
% Wyznaczanie parametrow sterowania optymalnego
[S, T] = parametry_optymalne(An, Bn, Cn, Qn, Rn, n);

D = - 0.5*inv(R)*H;
xd = inv(2*Q - 0.5*H'*inv(R)*H) * (0.5*H'*inv(R)*r - q);
ud = - 0.5*inv(R)*(r + H*xd);

% Symulacja

Amp = 0.0025;              % Amplituda zaklocen
x0 = [5.8;-2];              % Stan poczatkowy
xa = x0;
ua = [];
j = [];

for k=1:20
    x = xa(:, size(xa, 2))
    xb = x - xs;
    xp = xb - xd;
    uop = -inv(Bn'*S*Bn + Rn)*Bn'*(0.5*T'+ S*(An*xp + Cn));
    uob = uop + D * xp + ud;
    uo = uob + us;
    ua = [ua uo];
    
    v = V0 + Amp*randn(n,1);
    xk1 = transf(x,uo,v,n);
    xa = [xa xk1];
    j = [j wskjak(xk1, uo)];
end

% wykres
figure;
hold on;
grid on;
plot(1:size(xa,2), xa);
plot(1:size(ua,2), ua, 'r');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xn]=transf(x,u,v,n,m)
% x - stan x[k]
% u - sterowanie u[k]
% v - zak³ócenie v[k]
% n - rozmiar stanu
% m - rozmiar sterowania
% xn - nowy stan x[k+1]
xn=zeros(n,1);
xn(1)=1.5*x(1)-0.1*u(1)+v(1);
xn(2)=((1+v(1))^2)*x(2)+(0.7+v(1))*u(1)-x(1)*(5.2-v(1));


function [jk]=wskjak(xn,u,n,m)
% xn - stan x[k+1]
% u - sterowanie u[k]
% n - rozmiar wektora stanu
% m - rozmiar wektora sterowania
jk=xn(1)^2+xn(2)^2+0.1*u(1)^2;


function [f]=fun1(xx,n,V0)
% xx - wektor zmiennych decyzyjnych, zawiera najpierw stan, potem sterowanie
% n - rozmiar stanu
% V0 - wartosc oczekiwana zaklocenia
s=size(xx);
m=s(1)-n;
xn=xx(1:n); % wyodrêbnienie stanu
u=xx(n+1:n+m); % wyodrêbnienie sterowania
f=wskjak(xn,u,n,m);



function [con_neq,con_eq]=con(xx,n,V0)
% xx - wektor zmiennych decyzyjnych, zawiera najpierw stan, potem sterowanie
% n - rozmiar stanu
% V0 - wartosc oczekiwana zaklocenia
s=size(xx);
m=s(1)-n;
xn=xx(1:n); % wyodrêbnienie stanu
u=xx(n+1:n+m); % wyodrêbnienie sterowania
con_neq=[]; % nie mamy ograniczeñ nierównoœciowych
con_eq(1)=0.5*xn(1) - 0.1*u(1);
con_eq(2)=0.7*u(1)-xn(1)*5.2;
% transf(xn,u,V0,n,m)-;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,B,G,C,R,r,Q,q,H] = model_liniowy(xs, us, vs)
% xs - stan w punkcie stacjonarnym
% us - sterowanie w punkcie stacjonarnym
% vs - zaklocenie w punkcie stacjonarnym

A = [1.5, 0; (-5.2 + vs(1)), (1+vs(1))^2];  % wspolczynniki przy xs
B = [-0.1, (0.7+vs(1))];         % wspolczynniki przy us 
G = [1, 0; 2*(1 - vs(1))*xs(1) + us(1) - xs(1), 0;]; %wspolczynniki przy vs
C= [0;0];

R = 0.1;
r = 0;
Q = eye(2);
q = 0;
H = [0,0];

function[S, T] = parametry_optymalne(An, Bn, Cn, Qn, Rn, n)

[F, S] = dlqr(An, Bn, Qn, Rn);
T = 2*Cn'*S*(eye(n) - Bn*inv(Bn'*S*Bn + Rn)*Bn'*S)*An*inv(eye(n) - An + Bn*inv(Bn'*S*Bn + Rn)*Bn'*S*An);