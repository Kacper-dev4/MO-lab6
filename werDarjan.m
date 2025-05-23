
% Filip Bazarnicki, Kacper Wach
clear all

N = 7;                                      % Liczba sterowań
%J = @(x,u) (1/2)*sum(x(1:N).^2 + 3*u.^2);   % Podstawowy wskaźnik jakości
J = @(x,u) (1/2)*sum(x(1:N).^2 + 2*u.^2);
x0 = 20;                                    
u0 = ones(1,N);
%u0 = [2, 3, -2, -1, 4, -4, 1];
x = zeros(1,N+1);
%a = [1,4,4,11,-1,4,4,4,4,4,4,4,4,4,4,4,4];
a = [2,3,3,11,1,3,3,3,3,3,3,3,3,3,3,3,3];

alfa = 0.7;
beta = 12;
eps = 0.1;

c = 1.75;
%t = [18 ,13, 13, 18,18];
%t = [18 ,13, 13, 18,18,13,13,13,13,13,13,13,13,13,13,13,13];
r = 70;
n = 140;
t = [r,n,n,r,r,n,n,n,n,n,n,n,n,n,n,n,n];
v = a;
T(1,:) = t;

iter = 1;
while true

    u = fminsearch(@(u) Jkara(u,x0,t,a,N),u0);
    x(1) = x0;
    for i=1:N
        x(i+1) = x(i) + 2*u(i);  
    end
    X(iter,:) = x;
    U(iter,:) = u;
    wskaznikJakosci(iter) = J(x,u);
    zmodyfikowany(iter) = Jkara(u,x0,t,a,N);
    r(1) = x(N+1) - v(1);
    r(2) = max(max(0, -u(1)-v(2)));
    r(3) = max(max(0, u(1)-v(3)));
    r(4) = x(4) - v(4);
    r(5) = u(4) - v(5);
    r(6) = max(max(0, -u(2)-v(6)));
    r(7) = max(max(0, u(2)-v(7)));
    r(8) = max(max(0, -u(3)-v(8)));
    r(9) = max(max(0, u(3)-v(9)));
    r(10) = max(max(0, -u(4)-v(10)));
    r(11) = max(max(0, u(4)-v(11)));
    r(12) = max(max(0, -u(5)-v(12)));
    r(13) = max(max(0, u(5)-v(13)));
    r(14) = max(max(0, -u(6)-v(14)));
    r(15) = max(max(0, u(6)-v(15)));
    r(16) = max(max(0, -u(7)-v(16)));
    r(17) = max(max(0, u(7)-v(17)));
    
    y = norm(v + r - a);
    if y <= eps
        break;
    end

    if y<=c
        v = a - r;
        c = alfa*c;
    else
        t = beta*t;
        v = a -(1/beta)*r;
    end
    iter = iter + 1;
    T(iter,:) = t;
end

figure
plot(wskaznikJakosci,'*','LineStyle','--')
xlabel('Iteracja')
ylabel('Wartość wskaźnika jakości')


figure
plot(zmodyfikowany,'*','LineStyle','--')
xlabel('Iteracja')
ylabel('Wartość wskaźnika jakości zmodyfikowanego')

function J = Jkara(u,x0,t,a,N)
x = zeros(1,N);
x(1) = x0;
for i=1:N
    x(i+1) = x(i) + 2*u(i);  
end

% J = (1/2)*sum(x(1:N).^2 + 3*u.^2) ...
% + (1/2) *norm(x(N+1) - a(1))^2 * t(1) ...
% + ((-u(1) -a(2)) .* max(0,-u(1) - a(2))) * t(2)...
% + ((u(1) -a(3)) .* max(0,u(1) - a(3))) * t(3)...
% + ((-u(2) -a(2)) .* max(0,-u(2) - a(2))) * t(6)...
% + ((u(2) -a(3)) .* max(0,u(2) - a(3))) * t(7)...
% + ((-u(3) -a(2)) .* max(0,-u(3) - a(2))) * t(8)...
% + ((u(3) -a(3)) .* max(0,u(3) - a(3))) * t(9)...
% + ((-u(4) -a(2)) .* max(0,-u(4) - a(2))) * t(10)...
% + ((u(4) -a(3)) .* max(0,u(4) - a(3))) * t(11)...
% + ((-u(5) -a(2)) .* max(0,-u(5) - a(2))) * t(12)...
% + ((u(5) -a(3)) .* max(0,u(5) - a(3))) * t(13)...
% + ((-u(6) -a(2)) .* max(0,-u(6) - a(2))) * t(14)...
% + ((u(6) -a(3)) .* max(0,u(6) - a(3))) * t(15)...
% + ((-u(7) -a(2)) .* max(0,-u(7) - a(2))) * t(16)...
% + ((u(7) -a(3)) .* max(0,u(7) - a(3))) * t(17)...
% + norm((x(4) - a(4))^2) * t(4) ...
% + norm((u(4) - a(5))^2) * t(5);

J = (1/2)*sum(x(1:N).^2 + 2*u.^2) ...
+ (1/2) *norm(x(N+1) - a(1))^2 * t(1) ...
+ ((-u(1) -a(2)) .* max(0,-u(1) - a(2))) * t(2)...
+ ((u(1) -a(3)) .* max(0,u(1) - a(3))) * t(3)...
+ ((-u(2) -a(2)) .* max(0,-u(2) - a(2))) * t(6)...
+ ((u(2) -a(3)) .* max(0,u(2) - a(3))) * t(7)...
+ ((-u(3) -a(2)) .* max(0,-u(3) - a(2))) * t(8)...
+ ((u(3) -a(3)) .* max(0,u(3) - a(3))) * t(9)...
+ ((-u(4) -a(2)) .* max(0,-u(4) - a(2))) * t(10)...
+ ((u(4) -a(3)) .* max(0,u(4) - a(3))) * t(11)...
+ ((-u(5) -a(2)) .* max(0,-u(5) - a(2))) * t(12)...
+ ((u(5) -a(3)) .* max(0,u(5) - a(3))) * t(13)...
+ ((-u(6) -a(2)) .* max(0,-u(6) - a(2))) * t(14)...
+ ((u(6) -a(3)) .* max(0,u(6) - a(3))) * t(15)...
+ ((-u(7) -a(2)) .* max(0,-u(7) - a(2))) * t(16)...
+ ((u(7) -a(3)) .* max(0,u(7) - a(3))) * t(17)...
+ norm((x(4) - a(4))^2) * t(4) ...
+ norm((u(4) - a(5))^2) * t(5);


end