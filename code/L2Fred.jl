using Jacobi
using SpecialFunctions
using GR
using XLSX
using LinearAlgebra

function GL_Source(x,D,GL_nodes_N,gam,sigma)
    a=D[1];b=D[2];

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x1=zglj(GL_nodes_N,-gam,0)#JacobiGL(-gam,0,GL_nodes_N);
    w1=wglj(x1,-gam,0);
    Z1=(x-a)/2*x1.+(x+a)/2;
    U1=excat_u(Z1,sigma);
    I1=((x-a)/2)^(1-gam)*U1'*w1;
    
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    x2=zglj(GL_nodes_N,0,-gam);
    w2=wglj(x2,0,-gam);
    Z2=(b-x)/2*x2.+(b+x)/2;
    U2=excat_u(Z2,sigma);
    I2=((b-x)/2)^(1-gam)*U2'*w2;
    
    II=I1+I2;
    return II
end

function excat_u(x,sigma)

    ut=sin.(pi*x).*exp.(x)#;(-(x.-1)).^(2.1).*
    ###################### Continuous %%%%%%%%%%%%%%%%%%%%%%%%%%
    #ut=(x.^2.*(-(x.-1))).^2);#.*(x.*(-(x.-1))).^(sigma/2)x.*.^(sigma/2)(sin.(pi*x)).*
    # ut=exp.(x).*(x.*(-(x.-1))).^(sigma);
    ###################### Space Singularity Ⅰ %%%%%%%%%%%%%%%%%%%%%%%%%%
    #=
    alpha=sigma;
    CU=2^(-alpha)*gamma(1/2)/(gamma(1+alpha/2)*gamma((1+alpha)/2));
        ut=CU*(x.*(-(x.-1))).^(alpha/2);
         ut=exp(t).*ut;
    =#
    ######################## Space Singularity Ⅱ %%%%%%%%%%%%%%%%%%%%%%%%%%
     #ut=exp(t).*(x.^(1/2)+(-(x.-1))).^(1/2)-1);
    
    return ut
end

function ff(x,gam,sigma,F1)
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      f=lam(x,gam).*excat_u(x,sigma)-F1;
      #f=excat_u(x,t,sigma).+f; 
      return f
end

function lam(x,gam)
    global D
    l=1/(1-gam)*((x.-D[1]).^(1-gam)+(-(x.-D[2])).^(1-gam));
    #l=l.+1;
    return l
end

function Gradmesh(D,N,GM_type,r)
    #global r;
    P=Array{BigFloat}(undef,1,N+1);
    a=D[1];b=D[2];
    
################     Left side
    if GM_type=="L" 
        for j=1:N+1
            P[j]=a+(BigFloat(b-a))*((BigFloat(j-1)/N))^r;
        end
################     Right side
    elseif GM_type=="R"
        for j=1:N+1
            P[j]=b-(BigFloat(b-a))*(1-(BigFloat(j-1)/N))^r;
        end
######################################   Double side 
    elseif GM_type=="D"
        for j=1:(N)/2
            j=Int64(j);
            P[j]=((BigFloat(b-a))/2)*(2*(BigFloat(j-1))/(N))^r.+a;
        end
        for j=(N)/2+1:N+1
            j=Int64(j);
            P[j]=-((BigFloat(b-a))/2)*(-2*(BigFloat(j-1))/(N).+2)^r.+b;
        end
    end
    return P
end

function Assemble_Space(N,P,D,GL_nodes_N,gam,sigma)
    A=zeros(BigFloat,N-1,N-1);F1=zeros(BigFloat,N-1,1);
   
    C_g=1/(1-gam)/(2-gam);
    for i=2:N
        xi=P[i];
        #%hi=P[i]-P[i-1];hi1=P(i+1)-P(i);

        for j=2:N
            xj=P[j];xj_1=P[j-1];xj1=P[j+1];
            hj=xj-xj_1; hj1=xj1-xj;
            Cj=[1/hj -1/hj-1/hj1 1/hj1];

            Xj=[P[j-1] P[j] P[j+1]];
            Dji=(abs.(Xj'.-xi)).^(2-gam);
            A[i-1,j-1]=(-C_g*Cj*Dji)[1];
            
        end

        A[i-1,i-1]=lam(xi,gam) +A[i-1,i-1];

        F1[i-1,1]= GL_Source(xi,D,GL_nodes_N,gam,sigma);

    end

return [A,F1]
end
function L2error(u,P,N,Gn,sigma)
    uu=[0;u;0];
    z = zglj(Gn, 0.0, 0.0);
    w = wglj(z, 0.0, 0.0);
    Sum=0;
    for n=1:N
        Pl=P[n];Pr=P[n+1];
        #########################
        uhL=uu[n];uhR=uu[n+1];
        z1=(z.+1)./2;
        f1=(uhR-uhL).*z1.+uhL;
        #########################
        z2=(Pr-Pl).*z1.+Pl;
        uut=excat_u(z2,sigma);
        #########################
        ef=(f1-uut).^2;
        S=sum(w.*ef);
        Sum=Sum+S*(Pr-Pl)/2;
    end
error=Sum^0.5;
    return error
end

    
function Nonlocal_Temporal(D,N,GM_type,r,gam,sigma,GL_nodes_N)

    P=Gradmesh(D,N,GM_type,r);
    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P_in=P[2:N];
    AM=Assemble_Space(N,P,D,GL_nodes_N,gam,sigma);
    A=AM[1];F1=AM[2]; #AA=Float64.(A);a,b=eigen(AA);maximum(a)/minimum(a)
##################################################################
    f=ff(P_in,gam,sigma,F1);
    U=A\f;

    err=L2error(U,P,N,5,sigma)
    #Ut=excat_u(P_in,sigma);

    #error=abs.(U-Ut);
    #err=maximum(error);

    return err
end

#setprecision(300)
global D,Result_e
D=[0 1];

GM_type="D";
GL_nodes_N=500;
m=3;
nn0=25
index=[0.3  0.7];
#index=0.5;0.5
#Result_e=[r (exp.(log(2)*((1:m).+nn0)))'];
r=0.1
Result_e=[r (exp.(log(2)*((1:m))).*nn0)'];
for mm=1:length(index)

    global Result_e
gam=index[mm];
sigma=0.5#0.5;#%4;
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 #1;
 #r=2/(1+sigma);#6/(2+sigma-2*gam);#%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
error=zeros(1,m);
for nn=1:m
    N=2^(nn)*nn0;#N=2^(nn+nn0);
error[nn]=Nonlocal_Temporal(D,N,GM_type,r,gam,sigma,GL_nodes_N);
end
rho=log.(error[1:m-1]./error[2:m])./log(2);
Result_e=[Result_e;gam error;0 0 rho'];
end
Result_e
#P=Result_e(1,2:end);

#=
XLSX.openxlsx("NResult_sinE_r15.xlsx", mode="w") do xf
#XLSX.openxlsx("NResult_Ex1x03_r4.xlsx", mode="w") do xf
    sheet=xf[1]
    sheet["A2"]=Result_e;
    #sheet["A1"]="Uniform"
    ##sheet["A1"]="r=2/(1+sigma)";
    sheet["A1"]="r=1.5";#"r=2/(1+sigma-gam)";#
    sheet["B1"]="GL_nodes_N=500";
    sheet["C1"]="ut=sin(pi*x1)exp(x)";#"ut=EX(x.*(-(x.-1))).^(sigma)";#
    #sheet["D1"]="Sigma=0.3";
end
=#