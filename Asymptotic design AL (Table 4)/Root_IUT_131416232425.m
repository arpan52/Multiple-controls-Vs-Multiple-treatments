tic
 % x = sym('x',[1 4])'
 % root_iut(x)

 fun = @EQ_131416232425
 x0 = [0.1,0.2,0.3,0.2,0.1,0.1]'
 lb=[0,0,0,0,0,0]'
 %ub=[0.5,0.5,0.5,0.5,0.5,0.5]'
 ub=[1,1,1,1,1,1]'
%x= lsqnonlin(fun,x0,lb,ub)



Aeq = [1, 1, 1, 1,1,1];
beq = 1;
A=[];
b=[];
nonlcon=[];
options = optimoptions('ga','FunctionTolerance',1e-15);
x = ga(fun,6,A,b,Aeq,beq,lb,ub,nonlcon,options)
%x = fmincon(fun,x0,[],[],[1,1,1,1],[1],lb,ub,options)
eq_value = EQ_131416232425(x)
toc