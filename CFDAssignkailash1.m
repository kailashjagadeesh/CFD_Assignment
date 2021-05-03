%CFD Course Assignment - kailash jagadeesh
%{
dz-element length, H- half height of the channel between 2 plates
Nz- total number of mesh points
z- z coordinate in space of the mesh points
den-density of the fluid, eps- permittivity, E- electric field strength
eta-flow consistency, sigw- wall electric potential,rot-rotational velocity
%}
Kset=[10,20,30,40,60];nset=[0.7,0.8,1,1.2,1.5];
eta=0.9*10^-3; sigw=-0.25; rot=100; eps=709*10^-12; 
E=10^4;dz=2.5*10^-6;H=10^-4;z=[0:dz:H];

[velx1,vely1]= mysolver(eta,-.1,Kset(2),0.0,nset(4),eps,E);
[velx2,vely2]= mysolver(eta,-.1,Kset(2),200,nset(4),eps,E);
[velx3,vely3]= mysolver(eta,-.1,Kset(2),400,nset(4),eps,E);
[velx4,vely4]= mysolver(eta,-.1,Kset(2),600,nset(4),eps,E);
[velx5,vely5]= mysolver(eta,-.1,Kset(2),800,nset(4),eps,E);

figure(1);
plot(z/H,velx1);
hold on;
plot(z/H,velx2);
hold on;
plot(z/H,velx3);
hold on;
plot(z/H,velx4);
hold on;
plot(z/H,velx5);
hold on;
legend('omega=0','omega=200','omega=400','omega=700','omega=800')
figure(2);
plot(z/H,vely1);
hold on;
plot(z/H,vely2);
hold on;
plot(z/H,vely3);
hold on;
plot(z/H,vely4);
hold on;
plot(z/H,vely5);
hold on;
legend('omega=0','omega=200','omega=400','omega=700','omega=800')

function [unorm,vnorm]= mysolver(eta,sigw,K,rot,n,eps,E)
dz=2.5*10^-6;H=10^-4;Nz=41; 
%initialising velocities to 0
u=zeros(Nz,1);
v=zeros(Nz,1);
z=[0:dz:H];alp=-(3.9571*sigw).*z;den=1.06*10^3;
ueo=(-eps*E*sigw)/eta;
k=K/H;
%stability condition
dt=0.01*(den/eta)*(dz^(n+1));
%calculation of zeta potential through every mesh point
for j=1:Nz
    if j~=1
    sig(j)=((4*sigw)/alp(j))*atanh(tanh(0.25*alp(j))*exp(k*(z(j)-H)));
    else 
        sig(j)=sigw*exp(-k*H);
    end
end
a=1;
%main loop to iterate through time, a is a counter variable
  while a~=0
        uold=u;
        vold=v;
        if n==1
           vis=eta.*ones(Nz,1);
           
        else
        %calculation of viscosity through space using central difference
        %scheme
        for j=2:Nz-1
           if (uold(j+1)-uold(j-1))^2+(vold(j+1)-vold(j-1))^2~=0
            vis(j)=eta*((2*dz)^(-n+1))*(((uold(j+1)-uold(j-1))^2+(vold(j+1)-vold(j-1))^2)^((n-1)/2));
           %if gradient dissappears then its a newtonian fluid behaviour       
           else 
               vis(j)=eta;
           end
           %using backward difference for the end element
           if ((-uold(Nz-1))^2+(-vold(Nz-1))^2)~=0
            vis(Nz)=eta*((dz)^(-n+1))*(((-uold(Nz-1))^2+(-vold(Nz-1))^2)^((n-1)/2));
           else 
               vis(Nz)=eta;
           end
        end
        end
        %loop to iterate for velocities through space
        for j=2:Nz-1
            u(j)=uold(j)+2*rot*dt*vold(j)+(dt/den)*(dz^-2)*(vis(j+1)*(uold(j+1)-uold(j))-vis(j)*(uold(j)-uold(j-1)))...
                -(eps*k^2*E*sig(j))*(dt/den);
            v(j)=vold(j)-2*rot*dt*uold(j)+(dt/den)*(dz^-2)*(vis(j+1)*(vold(j+1)-vold(j))-vis(j)*(vold(j)-vold(j-1)));
        end
        %boundary conditions
            u(1)=u(2);
            v(1)=v(2);
            a=0;
        %checking for convergence
            for j=1:Nz-1
            if abs((u(j)-uold(j))/u(j))>10^-15||abs((v(j)-vold(j))/v(j))>10^-15
                
                a=a+1;
            end
            end
  end
    
    unorm=u./ueo;
    vnorm=v./ueo;
end
    

            
    


