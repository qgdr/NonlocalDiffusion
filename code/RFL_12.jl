using SpecialFunctions
using LinearAlgebra
using GR
using XLSX

function function_f(x,alpha)
    f=1;
    return f
end

function excat_u(x,alpha)
CC=2^(-alpha)*gamma(1/2)/(gamma(1+alpha/2)*gamma((1+alpha)/2));
ut=CC*(-x.*(x.-1)).^(alpha/2);
return ut
end


function Gradmesh(D,alpha,N,GM_type)
    global r
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

function Collocation_RFL(D,N,alpha)

        CR=1/(2*cos(pi*alpha/2)*gamma(2-alpha));

        A=zeros(BigFloat,N-1,N-1);f=zeros(BigFloat,N-1,1);
        P=Gradmesh(D,alpha,N,"D");

        C_a=1/(2-alpha)/(3-alpha);
        for i=2:N
            xi=P[i];xi_1=P[i-1];xi1=P[i+1];
            hi=xi-xi_1;hi1=xi1-xi;Hi=2/(hi+hi1);
            Xi=[P[i-1] P[i] P[i+1]];
            Ci=[1/hi -1/hi-1/hi1 1/hi1];

            for j=2:N
            xj=P[j];xj_1=P[j-1];xj1=P[j+1];
            hj=xj-xj_1;hj1=xj1-xj;
            Xj=[P[j-1] P[j] P[j+1]];
            Cj=[1/hj -1/hj-1/hj1 1/hj1];
            ###########################################
            Dij=abs.(Xj.-Xi').^(3-alpha);
            A[i-1,j-1]=(CR*C_a*Hi*Ci*Dij*Cj')[1];
            end
            f[i-1]=function_f(xi,alpha);
        end
        u=A\f#inv(A)*f;
        t=P[2:N];
        ut=excat_u(t,alpha);
        error=maximum(abs.(u-ut));

    return error
end

#setprecision(300)
global r,Result_e
D=[0,1]

m=6;
index=[1.1 1.3 1.5 1.7 1.9];
Result_e=[0 (exp.(log(2).*((1:m).+4)))']#zeros(BigFloat,1,m)#;
for mm=1:length(index)#,2
    global Result_e,r
alpha=index[mm];

r=1;
#r=4/alpha
#r=6/alpha;

error=zeros(BigFloat,1,m);
for nn=1:m
error[nn]=Collocation_RFL(D,2^(nn+4),alpha)

end
rho=log.(error[1:m-1]./error[2:m])./log(2);
Result_e=[Result_e;alpha error;0 0 rho'];
end
Result_e


# XLSX.openxlsx("Result_U.xlsx", mode="w") do xf
#     sheet=xf[1]
#     sheet["A2"]=Result_e;
#     #sheet["A2"]=["alpha" 0 (2*ones(Int64,1,m)).^(t.+4)']
#     sheet["A1"]="Uniform"
#     #sheet["A1"]="r=4/alpha" 
# end
#==#