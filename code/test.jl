using SpecialFunctions
using LinearAlgebra
using GR
using XLSX

function ff(x,s)
    f=1;
    #=a=2*s;
    f=1/(a*(1-a))*x.^(1-a)-2/(a*(1-a)*(2-a))*x.^(2-a)+1/(a*(1-a))*(1-x).^(1-a)-2/(a*(1-a)*(2-a))*(1-x).^(2-a)=#
    return f
end

function excat_u(x,s)
    C=(2*s*gamma(1/2))/(2*pi^(1/2)*gamma(1-s)*gamma(1+s));
    u=C*(-x.*(x.-1)).^s;
    #u=-x.*(x.-1);
    return u
end

function Gradmesh(D,s,N)
    P=Array{BigFloat}(undef,1,N+1);
    #P=Array{Float64}(undef,1,N+2);
    a=D[1];b=D[2];
    #alpha=3*(1-s)/s;
    #alpha=2*(1-s)/s;
    alpha=(1-s)/s;
    #alpha=(1-s)/(2*s);
    #alpha=1;
   for j=1:(N)/2
        j=Int64(j);
        P[j]=((BigFloat(b-a))/2)*(2*(BigFloat(j-1))/(N))^alpha.+a;
    end
    for j=(N)/2+1:N+1
        j=Int64(j);
        P[j]=-((BigFloat(b-a))/2)*(-2*(BigFloat(j-1))/(N).+2)^alpha.+b;
    end
    return P
end

function collection1D(D,N,s)
   h=1/N;
    P=Gradmesh(D,s,N);
    
    A=zeros(BigFloat,N-1,N-1);f=zeros(BigFloat,N-1,1);
    
    if s!=0.5
    C_s=1/(2*s*(1-2*s));
    end
    
    for i=2:N
        x_i=P[i];
        h_i=P[i]-P[i-1];h_i1=P[i+1]-P[i];
        
        eta_i=(h_i^(2-2*s)+h_i1^(2-2*s))/(1-2*s);
        for j=2:N
            h_j=P[j]-P[j-1];
            h_j1=P[j+1]-P[j];
            c_j=1/h_j;
            c_j1=1/h_j1;
            C_j=[c_j -c_j-c_j1 c_j1];            
            X_j=[P[j-1] P[j] P[j+1]];
            D_j_i=(abs.(X_j'.-x_i)).^(1-2*s);           
            if s>=0.5
                if j==i-1
                    D_j_i[3]=0;
                elseif j==i
                    D_j_i[2]=0;
                elseif j==i+1
                    D_j_i[1]=0;
                end
            end
            
            if j==i-1
                D_j_i[3]=0;
            elseif j==i
                D_j_i[2]=0;
            elseif j==i+1
                D_j_i[1]=0;
        end                      
        A[i-1,j-1]=(C_s*C_j*D_j_i)[1];
           #= if s>0.5
                if j==i-1
                    A[i-1,j-1]=A[i-1,j-1]+eta_i/(h_i*(h_i+h_i1));
                elseif j==i
                    A[i-1,j-1]=A[i-1,j-1]-eta_i/(h_i*(h_i1));
                elseif j==i+1
                    A[i-1,j-1]=A[i-1,j-1]+eta_i/(h_i1*(h_i+h_i1));
                end
        end=#
        end
        f[i-1]=ff(x_i,s);    
    end
    u=A\f;
    
    t=P[2:N]';
    ut=excat_u(t,s);
    err=abs.(u-ut')
    error=maximum(err)
    #plot(t,u,t,ut)
    return err[1],err[Int64(N/4)],err[Int64(N/2)],error
end

global Result,Result11,Result44
setprecision(600)
index=[0.05 0.15 0.25 0.35 0.45];
#index=[0.05 0.25 0.45]
D=[0 1];
m=6

Result=zeros(BigFloat,1,m+2);E=zeros(BigFloat,1,m+1);
Result11=zeros(BigFloat,1,m+1);Result44=zeros(BigFloat,1,m+1);
#nn=1
for nn=1:size(index,2)
    global Result,E,Result11,Result44
    s=index[nn];
    error1=zeros(BigFloat,1,m)
    errorN4=zeros(BigFloat,1,m)
    errorN=zeros(BigFloat,1,m)
    error=zeros(BigFloat,1,m)

    for n=1:m                                                                                                                                                                                              m
    N=2^(4+n);
    err=collection1D(D,N,s);
    error1[n]=err[1];errorN4[n]=err[2];
    errorN[n]=err[3];
    error[n]=err[4];
    end
    rho1=log.(error1[1:m-1]./error1[2:m])./log(2)
    rhoN4=log.(errorN4[1:m-1]./errorN4[2:m])./log(2)
    rhoN=log.(errorN[1:m-1]./errorN[2:m])./log(2)
    rho=log.(error[1:m-1]./error[2:m])./log(2)
    Result1=[2*s "error1" error1;0 0 0 rho1';0 "errorN" errorN;0 0 0 rhoN';0 "errorN_8" errorN4;0 0 0 rhoN4';0 "error" error;0 0 0 rho'];
    Result=[Result;Result1];
    Result11=[Result11;2*s error1;0 0 rho1'];
    Result44=[Result44;2*s errorN4;0 0 rhoN4'];
    E=[E;2*s error;"x1" error1;"N/8" errorN4;"N/2" errorN]
end

#=
t=1:m;
XLSX.openxlsx("Result_t2.xlsx", mode="w") do xf
    sheet=xf[1]
    sheet["A2"]=Result;
    sheet["A2"]=["alpha" (2*ones(Int64,1,m)).^(t.+4)']
    #sheet["A1"]="Uniform"
    sheet["A1"]="r=(1-s)/s" 
end
=#