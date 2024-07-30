function L = EQ_1424343536(y)

x=y';

n= sqrt(58);
del= 4;
sigma = 5.16;

 e1=[1,0,0,0,0,0]';
 e2=[0,1,0,0,0,0]';
 e3=[0,0,1,0,0,0]';
 e4=[0,0,0,1,0,0]';
 e5=[0,0,0,0,1,0]';
 e6=[0,0,0,0,0,1]';

 e14=e1+e4;
 e24=e2+e4;
 e34=e3+e4;
 e35=e3+e5;
 e36=e3+e6;


 A14=e1*e4';
 A24=e2*e4';
 A34=e3*e4';
 A35=e3*e5';
 A36=e3*e6';



 mu14 = sqrt((x'*A14*x)/(x'*e14));
 mu24 = sqrt((x'*A24*x)/(x'*e24));
 mu34 = sqrt((x'*A34*x)/(x'*e34));
 mu35 = sqrt((x'*A35*x)/(x'*e35));
 mu36 = sqrt((x'*A36*x)/(x'*e36));


p14=normpdf(norminv(0.95)- (n*del/sigma)*mu14)-normpdf(-norminv(0.95)- (n*del/sigma)*mu14);
p24=normpdf(norminv(0.95)- (n*del/sigma)*mu24)-normpdf(-norminv(0.95)- (n*del/sigma)*mu24);
p34=normpdf(norminv(0.95)- (n*del/sigma)*mu34)-normpdf(-norminv(0.95)- (n*del/sigma)*mu34);
p35=normpdf(norminv(0.95)- (n*del/sigma)*mu35)-normpdf(-norminv(0.95)- (n*del/sigma)*mu35);
p36=normpdf(norminv(0.95)- (n*del/sigma)*mu36)-normpdf(-norminv(0.95)- (n*del/sigma)*mu36);

L14_1 = e14'*x*x'*A14';
L14_2 = x'*A14*(e14'*x*eye(6)-x*e14');
L14_3 = 2*sqrt(x'*A14*x)*(x'*e14)^(3/2);
L14 = (L14_1+L14_2)/L14_3;

L24_1 = e24'*x*x'*A24';
L24_2 = x'*A24*(e24'*x*eye(6)-x*e24');
L24_3 = 2*sqrt(x'*A24*x)*(x'*e24)^(3/2);
L24 = (L24_1+L24_2)/L24_3;

L34_1 = e34'*x*x'*A34';
L34_2 = x'*A34*(e34'*x*eye(6)-x*e34');
L34_3 = 2*sqrt(x'*A34*x)*(x'*e34)^(3/2);
L34 = (L34_1+L34_2)/L34_3;

L35_1 = e35'*x*x'*A35';
L35_2 = x'*A35*(e35'*x*eye(6)-x*e35');
L35_3 = 2*sqrt(x'*A35*x)*(x'*e35)^(3/2);
L35 = (L35_1+L35_2)/L35_3;

L36_1 = e36'*x*x'*A36';
L36_2 = x'*A36*(e36'*x*eye(6)-x*e36');
L36_3 = 2*sqrt(x'*A36*x)*(x'*e36)^(3/2);
L36 = (L36_1+L36_2)/L36_3;

L1 = (p14*L14)+(p24*L24)+(p34*L34)+(p35*L35)+(p36*L36);
L = abs(L1(1))+abs(L1(2))+abs(L1(3))+abs(L1(4))+abs(L1(5))+abs(L1(6));

%L = abs([L1(1),L1(2),L1(3),L1(4),x'*ones(4,1)-1]);
%L =abs(L1(1))+abs(L1(2))+abs(L1(3))+abs(L1(4))+abs(x'*[1,0,-1,0]')+abs(x'*[0,1,0,-1]')
end
