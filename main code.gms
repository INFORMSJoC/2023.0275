*----------------------------------------------------------------------------------------------------*
*----------------------------------------------------------------------------------------------------*
*Incorporating promotional effects in sales planning of retail industry using geometric programming  *
*----------------------------------------------------------------------------------------------------*
*----------------------------------------------------------------------------------------------------*

***************************************************************************************************
***************************************************************************************************
*Convex Promotion Optimization  GEOMETRIC PROGRAMMING APPROACH                                    *
***************************************************************************************************
***************************************************************************************************
*Execseed = gmillisec(jnow);
Sets
i  items or SKUs /1*25/
t  weeks in the selling season /1*5/
k /1*3/
v  promotion vehicles /1*2/
m/1*2/
ist(t) ;
alias(i,j);
alias(t,to);
set yy(i,k,t);
Option  yy(i:k,t);
option seed=129

Parameters
S(i) separation time between two consecutive promotions
K_i(i) total tumber of promotion prices
L(i) upper bound on the promotion frequnecy
M_i(i) memory parameter
C(t) upper bound on the number of promotion in the category
Ut(t) upper bound on the number of vehicles
CV(v) upper bound on the number of times of using vehicle v       /1  1,2 1/
beta(v,t) promotion vehicle effect
cost(i,t) unit cost
b0(i)      self-price elasticity
L_total Total number of promotion frequnecy
sigma(i,j)   cross-item coefficient
a(i,t) seasonality effect
b(i,m) past prices effect
objective_rr;


L_total=5;
a(i,t)=uniform(5,7);
b0(i)=uniform(2,7);
b(i,m) =uniform(1,b0(i));
beta(v,t)=0.5;
cost(i,t)=uniform(0.45,0.65);
Ut(t)=3;
c(t)=UNIFORMINT(1,10);
M_i(i)=2;
L(i)=UNIFORMINT(1,10);
s(i)=0;
K_i(i)=3;
sigma(i,j)=uniform(-2,1.5);

table q(i,k)
table q(i,k)
       1     2      3
1      1    0.95   0.9
2      1    0.95   0.9
3      1     0.9   0.85
4      1    0.95   0.9
5      1    0.95   0.9
6      1    0.95   0.9
7      1     0.9   0.85
8      1    0.95   0.9
9      1    0.95   0.9
10     1    0.95   0.9
11     1     0.9   0.85
12     1    0.95   0.9
13     1    0.95   0.9
14     1     0.9   0.85
15     1    0.95   0.9
16     1    0.95   0.9
17     1    0.95   0.9
18     1    0.95   0.9
19     1     0.9   0.85
20     1    0.95   0.9
21     1    0.95   0.9
22     1    0.95   0.9
23     1    0.95   0.85
24     1    0.95   0.9
25     1    0.95   0.9

;


free variables
vp(i,t)
yp(i,t)
zp(i,t)
hp(i,t)
pp(i,t)
miop(v,t)
rr;
binary variables
gama(i,k,t)
x(v,t);
equations
obj
cons1
cons2
cons3
cons4
cons5
cons7
cons8
cons9
cons10
cons11
cons12
cons13
cons14;

obj.. rr=e=log(sum((t,i),exp(-vp(i,t)-yp(i,t)-zp(i,t)-hp(i,t)-sum(v,miop(v,t)))));

cons1(i,t).. exp(-pp(i,t)+vp(i,t))+cost(i,t)*exp(-pp(i,t))=l=1;

cons2(i,t).. yp(i,t)-a(i,t)+b0(i)*pp(i,t)=e=0;

cons3(i,t)..zp(i,t)-sum(m $ (ord(m) le M_i(i) and b(i,m) ne 0 ), b(i,m)*pp(i,t-(ord(m))))=e=0;

cons4(i,t).. hp(i,t)-sum(j $ (ord(j) ne ord(i)),sigma(j,i)*pp(j,t))=e=0;

cons5(v,t)..miop(v,t)=e=(log(1+beta(v,t)))*x(v,t);

cons7(i,t).. pp(i,t)+1-sum(k $( ord(k) le k_i(i) and q(i,k) ne 0),q(i,k)*gama(i,k,t))=e=0;

cons8(i)..sum((t,k) $ (ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=l(i);

cons9(i,t)..sum((to,k) $(ord(to) le ord(t)+s(i) and ord(to) ge ord(t) and ord(k) ge 2 and ord(k) le k_i(i)),gama(i,k,to))=l=1;

cons10.. sum((i,t,k) $ (ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=l_total;

cons11(t)..sum((i,k) $(ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=c(t);

cons12(v).. sum(t,x(v,t))=l=Cv(v);

cons13(t)..sum(v,x(v,t))=l=Ut(t);

cons14(i,t).. sum(k $ (ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=e=1;


model promotion_optimization as /all/;
option optca=0;
option optcr=0;
option reslim=10000000;
option MINLP=baron;
solve promotion_optimization using minlp min rr;
objective_rr=exp(rr.l);
display rr.l,gama.l,pp.l,vp.l,x.l,beta,cost,a,b0,b,sigma,objective_rr;
execute_unload "results1.gdx" rr.l,gama.l,pp.l,vp.l,x.l,beta,cost,a,b0,b,sigma,q




****************************************************************************************************
****************************************************************************************************
*Non-convex Promotion Optimization                                                                 *
****************************************************************************************************
****************************************************************************************************

*Execseed = gmillisec(jnow);
Sets
i  items or SKUs /1*25/
t  weeks in the selling season /1*5/
k /1*3/
v  promotion vehicles /1*2/
m/1*2/
ist(t) ;
alias(i,j);
alias(t,to);
set yy(i,k,t);
Option  yy(i:k,t);
option seed=129

Parameters
S(i) separation time between two consecutive promotions
K_i(i) total tumber of promotion prices
L(i) upper bound on the promotion frequnecy
M_i(i) memory parameter
C(t) upper bound on the number of promotion in the category
Ut(t) upper bound on the number of vehicles
CV(v) upper bound on the number of times of using vehicle v       /1  1,2 1/
beta(v,t) promotion vehicle effect
cost(i,t) unit cost
b0(i)      self-price elasticity
L_total Total number of promotion frequnecy
sigma(i,j)   cross-item coefficient
a(i,t) seasonality effect
b(i,m) past prices effect
objective_rr;


L_total=5;
a(i,t)=uniform(5,7);
b0(i)=uniform(2,7);
b(i,m) =uniform(1,b0(i));
beta(v,t)=0.5;
cost(i,t)=uniform(0.45,0.65);
Ut(t)=3;
c(t)=UNIFORMINT(1,10);
M_i(i)=2;
L(i)=UNIFORMINT(1,10);
s(i)=0;
K_i(i)=3;
sigma(i,j)=uniform(-2,1.5);

table q(i,k)
       1     2      3
1      1    0.95   0.9
2      1    0.95   0.9
3      1     0.9   0.85
4      1    0.95   0.9
5      1    0.95   0.9
6      1    0.95   0.9
7      1     0.9   0.85
8      1    0.95   0.9
9      1    0.95   0.9
10     1    0.95   0.9
11     1     0.9   0.85
12     1    0.95   0.9
13     1    0.95   0.9
14     1     0.9   0.85
15     1    0.95   0.9
16     1    0.95   0.9
17     1    0.95   0.9
18     1    0.95   0.9
19     1     0.9   0.85
20     1    0.95   0.9
21     1    0.95   0.9
22     1    0.95   0.9
23     1    0.95   0.85
24     1    0.95   0.9
25     1    0.95   0.9
;

free variables
vp(i,t)
yp(i,t)
zp(i,t)
hp(i,t)
pp(i,t)
miop(v,t)
rr;
binary variables
gama(i,k,t)
x(v,t);
equations
obj
cons1
cons2
cons3
cons4
cons5
cons7
cons8
cons9
cons10
cons11
cons12
cons13
cons14;

obj.. rr=e=log(sum((t,i),exp(-vp(i,t)-yp(i,t)-zp(i,t)-hp(i,t)-sum(v,miop(v,t)))));

cons1(i,t).. exp(-pp(i,t)+vp(i,t))+cost(i,t)*exp(-pp(i,t))=l=1;

cons2(i,t).. yp(i,t)-a(i,t)+b0(i)*pp(i,t)=e=0;

cons3(i,t)..zp(i,t)-sum(m $ (ord(m) le M_i(i) and b(i,m) ne 0 ), b(i,m)*pp(i,t-(ord(m))))=e=0;

cons4(i,t).. hp(i,t)-sum(j $ (ord(j) ne ord(i)),sigma(j,i)*pp(j,t))=e=0;

cons5(v,t)..miop(v,t)=e=(log(1+beta(v,t)))*x(v,t);

cons7(i,t).. exp(pp(i,t))-sum(k $( ord(k) le k_i(i) and q(i,k) ne 0),q(i,k)*gama(i,k,t))=e=0;

cons8(i)..sum((t,k) $ (ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=l(i);

cons9(i,t)..sum((to,k) $(ord(to) le ord(t)+s(i) and ord(to) ge ord(t) and ord(k) ge 2 and ord(k) le k_i(i)),gama(i,k,to))=l=1;

cons10.. sum((i,t,k) $ (ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=l_total;

cons11(t)..sum((i,k) $(ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=c(t);

cons12(v).. sum(t,x(v,t))=l=Cv(v);

cons13(t)..sum(v,x(v,t))=l=Ut(t);

cons14(i,t).. sum(k $ (ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=e=1;


model promotion_optimization as /all/;
option optca=0;
option optcr=0;
option reslim=10000000;
option MINLP=baron;
solve promotion_optimization using minlp min rr;
objective_rr=exp(rr.l);
display rr.l,gama.l,pp.l,vp.l,x.l,beta,cost,a,b0,b,sigma,objective_rr;

execute_unload "results.gdx" rr.l,gama.l,pp.l,vp.l,x.l,beta,cost,a,b0,b,sigma,q



***************************************************************************************************
***************************************************************************************************
*LAGRANGIAN DECOMPOSITION ALGORITHM                                                               *
***************************************************************************************************
***************************************************************************************************

*execseed = gmillisec(jnow);
Sets
i  items or SKUs /1*25/
t  weeks in the selling season /1*3/
k /1*3/
v  promotion vehicles /1*3/
m/1*3/
ist(t) ;
alias(i,j);
alias(t,to);
set yy(i,k,t);
Option  yy(i:k,t);
option seed=1234;



Parameters
S(i) separation time between two consecutive promotions
K_i(i) total tumber of promotion prices
L(i) upper bound on the promotion frequnecy
M_i(i) memory parameter
C(t) upper bound on the number of promotion in the category
Ut(t) upper bound on the number of vehicles
CV(v) upper bound on the number of times of using vehicle v       /1  2,2 2,3 1/
beta(v,t) promotion vehicle effect
cost(i,t) unit cost
b0(i)      self-price elasticity
L_total Total number of promotion frequnecy
sigma(i,j)   cross-item coefficient
a(i,t) seasonality effect
b(i,m) past prices effect
gp_obj(t)
sum_gp ;


L_total=20;
a(i,t)=uniform(5,7);
b0(i)=uniform(2,7);
beta(v,t)=uniform(0.5,1);
cost(i,t)=uniform(0.45,0.65);
loop(i,
loop(j,
sigma(i,j)=UNIFORM(-1.5,2);
sigma(j,i)=sigma(i,j);
);
sigma(i,i)=0
);

Ut(t)=UNIFORMINT(1,10);
c(t)=UNIFORMINT(1,5);
M_i(i)=3;
L(i)=UNIFORMINT(1,10);
s(i)=0;
b(i,m) =uniform(1,b0(i));
K_i(i)=3;





table q(i,k)
       1     2      3
1      1    0.95   0.9
2      1    0.95   0.9
3      1     0.9   0.85
4      1    0.95   0.9
5      1    0.95   0.9
6      1    0.95   0.9
7      1     0.9   0.85
8      1    0.95   0.9
9      1    0.95   0.9
10     1    0.95   0.9
11     1     0.9   0.85
12     1    0.95   0.9
13     1    0.95   0.9
14     1     0.9   0.85
15     1    0.95   0.9
16     1    0.95   0.9
17     1    0.95   0.9
18     1    0.95   0.9
19     1     0.9   0.85
20     1    0.95   0.9
21     1    0.95   0.9
22     1    0.95   0.9
23     1    0.95   0.85
24     1    0.95   0.9
25     1    0.95   0.9
;



*---------------------------------------------------------------------
* RELAXED SUB-MODEL 1 "GP"
*---------------------------------------------------------------------
Parameter
   landap(v,t)   'Lagrangian multipliers'
   landapp(v,t)  'Lagrangian multipliers'
   kapap(i,t)    'Lagrangian multipliers'
   kapapp(i,t)   'Lagrangian multipliers'
   ;


Variable z1 'relaxed objective of GP';
variable gpp;
free variables
vp(i,t)
yp(i,t)
zp(i,t)
hp(i,t)
pp(i,t)
miop(v,t);

Equation
   cons1
   cons2
   cons3
   cons4
 cons55
 cons66
   cons77
   cons88
*cons99
*cons1010
   defz1     'definition of z1';

defz1..      z1=e=log(sum((t,i),exp(-vp(i,t)-yp(i,t)-zp(i,t)-hp(i,t)-sum(v,miop(v,t))))+sum((v,t),landap(v,t)*exp(miop(v,t))
                      +landapp(v,t)*exp(-miop(v,t)))+sum((i,t),kapap(i,t)*exp(pp(i,t))+kapapp(i,t)*exp(-pp(i,t))));
cons1(i,t).. exp(-pp(i,t)+vp(i,t))+cost(i,t)*exp(-pp(i,t))=l=1;
cons2(i,t).. yp(i,t)-a(i,t)+b0(i)*pp(i,t)=e=0;
cons3(i,t)..zp(i,t)-sum(m $ (ord(t) gt ord(m)), b(i,m)*pp(i,t-(ord(m))))=e=0;
cons4(i,t).. hp(i,t)-sum(j $ (ord(j) ne ord(i)),sigma(j,i)*pp(j,t))=e=0;
cons55(v,t)..miop(v,t)=l=1;
cons66(v,t)..miop(v,t)=g=0;
cons77(i,t)..pp(i,t)=l=1;
cons88(i,t)..pp(i,t)=g=-1;
*cons99(i,t)..hp(i,t)=l=2;
*cons1010(i,t)..hp(i,t)=g=0;

Model GP / cons1,cons2,cons3,cons4,cons77,cons55,cons66,cons88,defz1 /;

*---------------------------------------------------------------------
* RELAXED SUB-MODEL 2 "IP1"
*---------------------------------------------------------------------

Parameter
   landa(v,t)   'Lagrangian multipliers'
   kapa(i,t)    'Lagrangian multipliers'

Variable z2 'relaxed objective of IP1';
binary variables
gama(i,k,t)
;
Equation
   cons5
   cons6
   cons7
   cons8
   cons11
   defz2     'definition of z2';

defz2..      z2=e=-sum((i,t,k)$( ord(k) le k_i(i) and q(i,k) ne 0),kapa(i,t)*q(i,k)*gama(i,k,t));

cons5(i)..sum((t,k) $ (ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=l(i);

cons6(i,t)..sum((to,k) $(ord(to) le ord(t)+s(i) and ord(to) ge ord(t) and ord(k) ge 2 and ord(k) le k_i(i)),gama(i,k,to))=l=1;

cons7.. sum((i,t,k) $ (ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=l_total;

cons8(t)..sum((i,k) $(ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=c(t);


cons11(i,t).. sum(k $ (ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=e=1;

Model IP1 / cons5,cons6,cons7,cons8,cons11,defz2 /;

*---------------------------------------------------------------------
* RELAXED SUB-MODEL 2 "IP2"
*---------------------------------------------------------------------

Variable z3 'relaxed objective of IP2';
binary variables
x(v,t);
Equation
   cons5
   cons6
   cons7
   cons8
   cons9
   cons10
   cons11
   defz3     'definition of z2';

defz3..      z3=e=-sum((t,v),landa(v,t)*(1+beta(v,t)*x(v,t)));


cons9(v).. sum(t,x(v,t))=l=Cv(v);

cons10(t)..sum(v,x(v,t))=l=Ut(t);


Model IP2 / cons9,cons10,defz3 /;





*---------------------------------------------------------------------
* ORIGINAL MODEL
*---------------------------------------------------------------------

positive variable
mio(v,t)
p(i,t)
vi(i,t);

*---------------------------------------------------------------------
* subgradient iterations
*---------------------------------------------------------------------
set iter /iter1*iter70/;
parameter stepsize1;
parameter stepsize2;
parameter phi /0.07/;
scalar noimprovement /0/;
parameter BUB ;
scalar upperbound /+INF/;
parameter theta(i,t);
parameter pi(v,t);
scalar norm1;
scalar norm2;
scalar  bestbound /-INF/;
parameter result(iter,*);
parameter results(iter,*);
parameter norm;
parameter lowerbound;
parameter lowerb;
PARAMETER  newbub;
parameter bub1;
parameter landapervious(v,t);
parameter kapaprevious(i,t);
scalar optimal /0.007/;


*--------------------------------------------------------
* initialize LR multipliers with relaxed duals
*--------------------------------------------------------
lowerb=0;
landap(v,t)=0 ;
landapp(v,t)=0;
kapap(i,t)=0 ;
kapapp(i,t)=0 ;
landa(v,t)=0.01;
kapa(i,t)=0.001;

*--------------------------------------------------------
* ITERATIONS
*--------------------------------------------------------

loop(iter,

*---------------------------------------------------------
* solve the lagrangian dual problem
*---------------------------------------------------------

loop((i,t),
if(kapa(i,t)>=0,
kapapp(i,t)=0;
kapap(i,t)=kapa(i,t);
else
kapapp(i,t)=abs(kapa(i,t));
kapap(i,t)=0;
);
);
loop((v,t),
if(landa(v,t)>=0,
landapp(v,t)=0;
landap(v,t)=landa(v,t);
else
landapp(v,t)=abs(landa(v,t));
landap(v,t)=0;
);
*kapa(i,t)= kapapp(i,t)+ kapap(i,t);
*landa(v,t)= landapp(v,t)+ landap(v,t);
);
option optca=0;
option optcr=0;
option NLP=mosek;
solve GP minimizing z1 using nlp;
display landa,landapp,landap , kapa, kapap, kapapp, beta ;
p.fx(i,t)=exp(pp.l(i,t));
mio.fx(v,t)=exp(miop.l(v,t));





option optca=0;
option optcr=0;
option MIP= cplex;
solve IP1 minimizing z2 using mip;

option optca=0;
option optcr=0;
option MIP= cplex;
solve IP2 minimizing z3 using mip;

LOWERB=exp(z1.l)+z2.l+z3.l;
theta(i,t) =p.l(i,t)-sum(k $( ord(k) le k_i(i) and q(i,k) ne 0),q(i,k)*gama.l(i,k,t));
pi(v,t)=mio.l(v,t)-1-(beta(v,t)*x.l(v,t));
norm1 = sum((i,t),sqr(theta(i,t)));
norm2= sum((v,t),sqr(pi(v,t)));

 if(LOWERB> bestbound ,
bestbound = LOWERB;
noimprovement = 0;
else
noimprovement = noimprovement + 1;
if (noimprovement >= 2,
phi = phi/10;
noimprovement = 0;
);
);


*--------------------------------------------------------
* FEASIBILITY PAHSE (UPPER bound)
*--------------------------------------------------------
loop((i,t,v),
if(p.l(i,t) ne (sum(k $( ord(k) le k_i(i) and q(i,k) ne 0),q(i,k)*gama.l(i,k,t))),
p.l(i,t)=sum(k $( ord(k) le k_i(i) and q(i,k) ne 0),q(i,k)*gama.l(i,k,t));
pp.l(i,t)=log(p.l(i,t));
vp.l(i,t)=pp.l(i,t)-cost(i,t);
hp.l(i,t)=sum(j $ (ord(j) ne ord(i)),sigma(j,i)*pp.l(j,t));
zp.l(i,t)=sum(m $ (ord(t) gt ord(m)), b(i,m)*pp.l(i,t-(ord(m))));
yp.l(i,t)=a(i,t)+b0(i)*pp.l(i,t);
*vi.l(i,t)=exp(vp.l(i,t));
);
if(mio.l(v,t) ne (1+beta(v,t)*x.l(v,t)),
mio.l(v,t)=1+beta(v,t)*x.l(v,t);
miop.l(v,t)=log(mio.l(v,t));
);
);

BUB1= sum((i,t),vcpower((p.l(i,t)-cost(i,t)),-1)*exp(-a(i,t)) *vcpower((p.l(i,t)),b0(i))
 *prod(j $ (ord(j) ne ord(i)),rpower((p.l(j,t)),-sigma(j,i)))
 *prod(v,vcpower((mio.l(v,t)),-1))
 *prod(m $ (ord(t) gt ord(m)) ,vcpower(((p.l(i,t-ord(m)))),-b(i,m))));

BUB=(sum((t,i),exp(-vp.l(i,t)-yp.l(i,t)-zp.l(i,t)-hp.l(i,t)-sum(v,miop.l(v,t)))));
newbub=exp(bub1);


if (BUB1<upperbound,
upperbound=BUB1;
);

*------------------------------------------------------
* calculate step size
*------------------------------------------------------
result(iter,'dual obj')=LOWERB;
result(iter, 'optimal obj')=optimal;
result(iter,'ub')=upperbound;
stepsize1 = phi*(upperbound-bestbound)/norm1;
stepsize2 = phi*(upperbound-bestbound)/norm2;

*------------------------------------------------------
* update duals u
*------------------------------------------------------

kapa(i,t)=kapa(i,t)+stepsize1*theta(i,t);
landa(v,t)=landa(v,t)+ stepsize2*pi(v,t);
*------------------------------------------------------
* converged ?
*------------------------------------------------------

  if(((upperbound-bestbound) < (0.01)*bestbound),
     break;

    );
display iter;
);
display result,bestbound,upperbound,kapa,landa,cost,GAMA.L,pp.l, lowerb, z1.l,z2.l,z3.l,x.l, miop.l,mio.l;

execute_unload "results1.gdx" result,kapa,landa,GAMA.L,pp.l,upperbound, lowerb, x.l,vp.l,beta,cost,a,b0,b,sigma



***************************************************************************************************
***************************************************************************************************
*LAGRANGIAN DECOMPOSITION ALGORITHM decomposed over t to reduce the running time of sub-model 1   *
***************************************************************************************************
***************************************************************************************************

*execseed = gmillisec(jnow);
Sets
i  items or SKUs /1*25/
t  weeks in the selling season /1*10/
k /1*3/
v  promotion vehicles /1*3/
m/1*3/
ist(t) ;
alias(i,j);
alias(t,to);
set yy(i,k,t);
Option  yy(i:k,t);
*option seed=123;



Parameters
S(i) separation time between two consecutive promotions
K_i(i) total tumber of promotion prices
L(i) upper bound on the promotion frequnecy
M_i(i) memory parameter
C(t) upper bound on the number of promotion in the category
Ut(t) upper bound on the number of vehicles
CV(v) upper bound on the number of times of using vehicle v       /1  0,2 0,3 0/
beta(v,t) promotion vehicle effect
cost(i,t) unit cost
b0(i)      self-price elasticity
L_total Total number of promotion frequnecy
sigma(i,j)   cross-item coefficient
a(i,t) seasonality effect
b(i,m) past prices effect
gp_obj(t)
sum_gp ;


L_total=20;
a(i,t)=uniform(5,7);
b0(i)=uniform(2,7);
beta(v,t)=uniform(0.5,1);
cost(i,t)=uniform(0.45,0.65);
loop(i,
loop(j,
sigma(i,j)=UNIFORM(-2,1.5);
sigma(j,i)=sigma(i,j);
);
sigma(i,i)=0
);

Ut(t)=UNIFORMINT(1,10);
c(t)=UNIFORMINT(1,5);
M_i(i)=3;
L(i)=UNIFORMINT(1,10);
s(i)=0;
b(i,m) =uniform(1,b0(i));
K_i(i)=3;

table q(i,k)
       1     2      3
1      1    0.95   0.9
2      1    0.95   0.9
3      1     0.9   0.85
4      1    0.95   0.9
5      1    0.95   0.9
6      1    0.95   0.9
7      1     0.9   0.85
8      1    0.95   0.9
9      1    0.95   0.9
10     1    0.95   0.9
11     1     0.9   0.85
12     1    0.95   0.9
13     1    0.95   0.9
14     1     0.9   0.85
15     1    0.95   0.9
16     1    0.95   0.9
17     1    0.95   0.9
18     1    0.95   0.9
19     1     0.9   0.85
20     1    0.95   0.9
21     1    0.95   0.9
22     1    0.95   0.9
23     1    0.95   0.85
24     1    0.95   0.9
25     1    0.95   0.9
;




*---------------------------------------------------------------------
* RELAXED SUB-MODEL 1 "GP"
*---------------------------------------------------------------------
Parameter
   landap(v,t)   'Lagrangian multipliers'
   landapp(v,t)  'Lagrangian multipliers'
   kapap(i,t)    'Lagrangian multipliers'
   kapapp(i,t)   'Lagrangian multipliers'
   ;


Variable z1 'relaxed objective of GP';
variable gpp;
free variables
vp(i,t)
yp(i,t)
zp(i,t)
hp(i,t)
pp(i,t)
miop(v,t);

Equation
   cons1
   cons2
   cons3
   cons4
   cons55
   cons66
   cons77
   cons88
*cons99
*cons1010
   defz1     'definition of z1';

defz1(ist)..      z1=e=(sum((i),exp(-vp(i,ist)-yp(i,ist)-zp(i,ist)-hp(i,ist)-sum(v,miop(v,ist))))+sum((v),landap(v,ist)*exp(miop(v,ist))
                      +landapp(v,ist)*exp(-miop(v,ist)))+sum((i),kapap(i,ist)*exp(pp(i,ist))+kapapp(i,ist)*exp(-pp(i,ist))));
cons1(i,ist).. exp(-pp(i,ist)+vp(i,ist))+cost(i,ist)*exp(-pp(i,ist))=l=1;
cons2(i,ist).. yp(i,ist)-a(i,ist)+b0(i)*pp(i,ist)=e=0;
cons3(i,t)$ist(t).. zp(i,t)- sum(m  $ (ord(t) gt ord(m)), b(i,m)*pp(i,t-ord(m)))=e=0 ;
cons4(i,t)$ist(t).. hp(i,t)-sum(j $ (ord(j) ne ord(i)),sigma(j,i)*pp(j,t))=e=0;
cons55(v,ist)..miop(v,ist)=l=1;
cons66(v,ist)..miop(v,ist)=g=0;
cons77(i,ist)..pp(i,ist)=l=1;
cons88(i,ist)..pp(i,ist)=g=-1;



Model GP / cons1,cons2,cons3,cons4,cons55,cons66,cons77,cons88,defz1 /;
set tt(t);
tt(t)=yes;
ist(t)=no ;
*---------------------------------------------------------------------
* RELAXED SUB-MODEL 2 "IP1"
*---------------------------------------------------------------------

Parameter
   landa(v,t)   'Lagrangian multipliers'
   kapa(i,t)    'Lagrangian multipliers'

Variable z2 'relaxed objective of IP1';
binary variables
gama(i,k,t)
;
Equation
   cons5
   cons6
   cons7
   cons8
   cons11
   defz2     'definition of z2';

defz2..      z2=e=-sum((i,t,k)$( ord(k) le k_i(i) and q(i,k) ne 0),kapa(i,t)*q(i,k)*gama(i,k,t));

cons5(i)..sum((t,k) $ (ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=l(i);

cons6(i,t)..sum((to,k) $(ord(to) le ord(t)+s(i) and ord(to) ge ord(t) and ord(k) ge 2 and ord(k) le k_i(i)),gama(i,k,to))=l=1;

cons7.. sum((i,t,k) $ (ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=l_total;

cons8(t)..sum((i,k) $(ord(k) ge 2 and ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=l=c(t);


cons11(i,t).. sum(k $ (ord(k) le k_i(i)and q(i,k) ne 0),gama(i,k,t))=e=1;

Model IP1 / cons5,cons6,cons7,cons8,cons11,defz2 /;

*---------------------------------------------------------------------
* RELAXED SUB-MODEL 2 "IP2"
*---------------------------------------------------------------------

Variable z3 'relaxed objective of IP2';
binary variables
x(v,t);
Equation
   cons5
   cons6
   cons7
   cons8
   cons9
   cons10
   cons11
   defz3     'definition of z2';

defz3..      z3=e=-sum((t,v),landa(v,t)*(+1+beta(v,t)*x(v,t)));


cons9(v).. sum(t,x(v,t))=l=Cv(v);

cons10(t)..sum(v,x(v,t))=l=Ut(t);


Model IP2 / cons9,cons10,defz3 /;


*---------------------------------------------------------------------
* ORIGINAL MODEL
*---------------------------------------------------------------------

positive variable
mio(v,t)
p(i,t)
vi(i,t);

*---------------------------------------------------------------------
* subgradient iterations
*---------------------------------------------------------------------
set iter /iter1*iter30/;
parameter stepsize1;
parameter stepsize2;
parameter phi /0.4/;
scalar noimprovement /0/;
parameter BUB ;
scalar upperbound /+INF/;
parameter theta(i,t);
parameter pi(v,t);
scalar norm1;
scalar norm2;
scalar  bestbound /-INF/;
parameter result(iter,*);
parameter results(iter,*);
parameter norm;
parameter lowerbound;
parameter landapervious(v,t);
parameter kapaprevious(i,t);
scalar optimal /0.00794/;



*--------------------------------------------------------
* initialize LR multipliers with relaxed duals
*--------------------------------------------------------
gp_obj(t)=0;
landap(v,t)=0 ;
landapp(v,t)=0;
kapap(i,t)=0 ;
kapapp(i,t)=0 ;
landa(v,t)=0;
kapa(i,t)=0;

*--------------------------------------------------------
* ITERATIONS
*--------------------------------------------------------

loop(iter,

*---------------------------------------------------------
* solve the lagrangian dual problem
*---------------------------------------------------------
loop((i,t,v),
if(kapa(i,t)>=0,
kapapp(i,t)=0;
kapap(i,t)=kapa(i,t);
*kapa(i,t)= kapap(i,t)+kapapp(i,t);
else
kapapp(i,t)=abs(kapa(i,t));
kapap(i,t)=0;
*kapa(i,t)= kapap(i,t)+kapapp(i,t);
);
if(landa(v,t)>0,
landapp(v,t)=0;
landap(v,t)=landa(v,t);
*landa(v,t)=landap(v,t)+landapp(v,t);
else
landapp(v,t)=abs(landa(v,t));
landap(v,t)=0;
*landa(v,t)=landap(v,t)+landapp(v,t);
);
);

option optca=0;
option optcr=0;
option NLP=mosek;
loop(tt,
ist(tt)=yes;
solve GP minimizing z1 using nlp;
gp_obj(tt)=z1.l
loop(t$ist(t) ,
if(ord(t)=1,
pp.fx(i,"1")=pp.l(i,"1");
elseif(ord(t)=2),
pp.fx(i,"2")=pp.l(i,"2");
elseif(ord(t)=3),
pp.fx(i,"3")=pp.l(i,"3");
elseif(ord(t)=4),
pp.fx(i,"4")=pp.l(i,"4");
elseif(ord(t)=5),
pp.fx(i,"5")=pp.l(i,"5");
elseif(ord(t)=6),
pp.fx(i,"6")=pp.l(i,"6");
elseif(ord(t)=7),
pp.fx(i,"7")=pp.l(i,"7");
elseif(ord(t)=8),
pp.fx(i,"8")=pp.l(i,"8");
elseif(ord(t)=9),
pp.fx(i,"9")=pp.l(i,"9");
);
);
ist(tt)=no;
);

p.fx(i,t)=exp(pp.l(i,t));
mio.fx(v,t)=exp(miop.l(v,t));

option optca=0;
option optcr=0;
option MIP= cplex;
solve IP1 minimizing z2 using mip;

option optca=0;
option optcr=0;
option MIP= cplex;
solve IP2 minimizing z3 using mip;

sum_gp=sum(t,gp_obj(t));

lowerbound=sum_gp+z2.l+z3.l;

if (lowerbound>= bestbound,
bestbound = lowerbound;
noimprovement = 0;
else
noimprovement = noimprovement + 1;
if (noimprovement >= 2,
phi = phi/10;
noimprovement = 0;
);
);


theta(i,t) = p.l(i,t)-sum(k $( ord(k) le k_i(i) and q(i,k) ne 0),q(i,k)*gama.l(i,k,t));
pi(v,t)=mio.l(v,t)-1-(beta(v,t)*x.l(v,t));
norm1 = sum((i,t),sqr(theta(i,t)));
norm2= sum((v,t),sqr(pi(v,t)));


*--------------------------------------------------------
* FEASIBILITY PAHSE (UPPER bound)
*--------------------------------------------------------
loop((i,t,v),
if(p.l(i,t) ne (sum(k $( ord(k) le k_i(i) and q(i,k) ne 0),q(i,k)*gama.l(i,k,t))),
p.l(i,t)=sum(k $( ord(k) le k_i(i) and q(i,k) ne 0),q(i,k)*gama.l(i,k,t));
*vi.l(i,t)=exp(vp.l(i,t));
);
if(mio.l(v,t) ne (1+beta(v,t)*x.l(v,t)),
mio.l(v,t)=(1+(beta(v,t)*x.l(v,t)));
);
);


BUB= sum((i,t),vcpower((p.l(i,t)-cost(i,t)),-1)*exp(-a(i,t)) *vcpower((p.l(i,t)),b0(i))
 *prod(j $ (ord(j) ne ord(i)),rpower((p.l(j,t)),-sigma(j,i)))
 *prod(v,vcpower((mio.l(v,t)),-1))
  *prod(m $ (ord(t) gt ord(m)) ,vcpower(((p.l(i,t-ord(m)))),-b(i,m))));



if (BUB<upperbound,
upperbound=BUB;
);

result(iter,'dual obj')=bestbound;
result(iter, 'optimal obj')=optimal;
result(iter,'ub')=upperbound;



*------------------------------------------------------
* calculate step size
*------------------------------------------------------

stepsize1 = phi*(upperbound-bestbound)/norm1;
stepsize2 = phi*(upperbound-bestbound)/norm2;

*------------------------------------------------------
* update duals u
*------------------------------------------------------

kapa(i,t)=kapa(i,t)+stepsize1*theta(i,t);
landa(v,t)=landa(v,t)+ stepsize2*pi(v,t);
*------------------------------------------------------
* converged ?
*------------------------------------------------------

 if(((upperbound-bestbound) < (1e-2)*bestbound),
    break;

);
display iter;
);

display result,bestbound,upperbound,kapa,landa,cost,GAMA.L,pp.l,sigma, z1.l,z2.l,z3.l,x.l, miop.l,mio.l,gp_obj,sum_gp;

execute_unload "results1.gdx" GAMA.L,pp.l,upperbound, x.l,vp.l,beta,cost,a,b0,b,sigma,Ut,c,L,bestbound
