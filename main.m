clear all

N = 7;                                      % Liczba sterowań
J = @(x,u) (1/2)*sum(x(1:N).^2 + 3*u.^2);   % Podstawowy wskaźnik jakości
x0 = 20;                                    
u0 = ones(1,N);
x = zeros(1,N+1);
a = [1,11,4,4,-1];

alfa = 0.5;
beta = 2;
eps = 0.1;

c = 2;
t = [0.5 ,0.5];
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

    r(1) = x(N+1) - v(1);
    r(2) = x(4) - v(2);
    r(3) = max(max(0, -u-v(3)));
    r(4) = max(max(0, u-v(4)));
    r(5) = u(4) - v(5);
    
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
plot(wskaznikJakosci)
xlabel('Iteracja')
ylabel('Wartość wskaźnika jakości')

function J = Jkara(u,x0,t,a,N)
x = zeros(1,N);
x(1) = x0;
for i=1:N
    x(i+1) = x(i) + 2*u(i);  
end

J = (1/2)*sum(x(1:N).^2 + 3*u.^2) ...
+ (1/2) *norm(x(N+1) - a(1))^2 * t(1) ...
+ norm((x(4) - a(2))^2) * t(1) ...
+ sum((-u -a(3)) .* max(0,-u - a(3))) * t(2)...
+ sum((u -a(4)) .* max(0,u - a(4))) * t(2)...
+ norm((u(4) - a(5))^2) * t(1);
end