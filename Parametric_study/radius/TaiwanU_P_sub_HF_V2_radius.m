clear;clc
q=1.6022e-19;
Na=1.25e21;
ni=1.18e16;
phiB=0.026*log(Na/ni)-0.0425;
dielectSi=11.9*8.85e-12;
dielectOx=4*8.85e-12;
R=10e-6;
Tox=0.1e-6;
h=100e-6
% charge in Oxide
Qox=0;
% range of fais
range=1050;
fais=zeros(range,1);
electricfield1=zeros(range,7);
electricfield2=zeros(range,7);
%tsv_lf=zeros(range,4);
%tsv_hf=zeros(range,4);
cox=2*pi*dielectOx*h/log(R/(R-Tox));
%%% Parametric study-Radius
R=[20e-6,10e-6,5e-6,2e-6,1e-6]
cox=2*pi*dielectOx*h./log(bsxfun(@rdivide, R, bsxfun(@minus,R,Tox)));

for par=1:5    

for i=2:range
% fais range -0.15 to 0.8
fais(i,1)=-0.23+0.001*i;
F=(0.026*1.6022e-19*Na/dielectSi)*((exp(-fais(i)/0.026)+fais(i)/0.026-1)+(ni^2/Na^2)*(exp(fais(i)/0.026)-fais(i)/0.026+1));
E=fais(i)/R(par)+sqrt((fais(i)/R(par))^2+2*F);
%positive direction for surface voltage
electricfield1(i,2)=fais(i);
if electricfield1(i,2)>0
    electricfield1(i,1)=-E;
else
     electricfield1(i,1)=E;
end
%charge amount
electricfield1(i,3)=2*pi*R(par)*h*dielectSi*electricfield1(i,1);
%voltage on liner
electricfield1(i,4)=-electricfield1(i,3)/cox(par);
% voltage of TSV
electricfield1(i,5)=electricfield1(i,4)+electricfield1(i,2);
%Cdep
electricfield1(i,6)=abs(electricfield1(i,3)-electricfield1(i-1,3))/0.001;
%Ctsv
electricfield1(i,7)=1/((1/electricfield1(i,6))+(1/cox(par)));
end
figure(1)
plot(electricfield1(:,5),electricfield1(:,7))

%%FOR TSV2
for i=2:range
    %charge on TSV2
    electricfield2(i,3)=Qox-electricfield1(i,3);
    %electric field
    E=electricfield2(i,3)/(2*pi*R(par)*h*dielectSi)
    electricfield2(i,1)=E;
    % voltage of liner
    electricfield2(i,4)=-electricfield2(i,3)/cox(par);
      for j=2:range-1
        if ((electricfield1(j,1)-E)*(electricfield1(j+1,1)-E)<0)
            %fais=electricfield1(j,2);
            %fais on TSV2
            electricfield2(i,2)=electricfield1(j,2)+(electricfield2(i,3)-electricfield1(j,3))*0.002/(electricfield1(j+1,3)-electricfield1(j-1,3));
            break;
        %else
        %errordlg('error')
        end
      end
% voltage of TSV
electricfield2(i,5)=electricfield2(i,4)+electricfield2(i,2);
%Cdep
electricfield2(i,6)=abs(electricfield2(i,3)-electricfield2(i-1,3))/abs(electricfield2(i,2)-electricfield2(i-1,2));
%Ctsv
electricfield2(i,7)=1/((1/electricfield2(i,6))+(1/cox(par)));
end
figure(2)
plot(electricfield2(:,5),electricfield2(:,7))

%% CV curve for TSV pair
for i=1:range
tsv_total_cap(i,1)=1/((1/electricfield1(i,7))+(1/electricfield2(i,7)));
tsv_total_voltage(i,1)=electricfield1(i,5)-electricfield2(i,5);
end
electricfield1(:,8)=smooth(electricfield1(:,7),10);
electricfield2(:,8)=smooth(electricfield2(:,7),10);
tsv_total_cap(:,2)=smooth(tsv_total_cap(:,1),10);
figure(3)
plot(tsv_total_voltage,tsv_total_cap)
figure(4)
plot(tsv_total_voltage,tsv_total_cap(:,2),tsv_total_voltage,electricfield1(:,7),tsv_total_voltage,electricfield2(:,8))
%% symerical transform
for i=1:range
    if tsv_total_voltage(i)<0
    tsv_lf(i,1)=tsv_total_voltage(i);
    tsv_lf(i,2)=electricfield1(i,7);
    tsv_lf(i,3)=electricfield2(i,8);
    tsv_lf(i,4)=tsv_total_cap(i,2);
    end
end
mat=size(tsv_lf)
for i=mat(1)+1:mat(1)*2
    tsv_lf(i,1)=-tsv_lf(2*mat(1)+1-i,1);
    tsv_lf(i,3)=tsv_lf(2*mat(1)+1-i,2);
    tsv_lf(i,2)=tsv_lf(2*mat(1)+1-i,3);
    tsv_lf(i,4)=tsv_lf(2*mat(1)+1-i,4);
end
plot(tsv_lf(:,1),tsv_lf(:,2),tsv_lf(:,1),tsv_lf(:,3),tsv_lf(:,1),tsv_lf(:,4))

%% HF for TSV1
for i=2:range-1
    if (electricfield1(i,2)-(2*phiB))*(electricfield1(i+1,2)-(2*phiB))<=0
        thresh=i
        break
    end
end
for i=2:range
    if electricfield1(i,2)>electricfield1(thresh,2)
        electricfield1(i,6)=electricfield1(thresh,6);
        electricfield1(i,7)=1/((1/electricfield1(i,6))+(1/cox(par)));
    end
end
figure(1)
plot(electricfield1(:,5),electricfield1(:,7))

%% HF for TSV2
for i=2:range-1
    if (electricfield2(i-1,2)>0)&((electricfield2(i,2)-(2*phiB))*(electricfield2(i+1,2)-(2*phiB))<=0)
        thresh=i
        break
    end
end
for i=2:range
    if electricfield2(i,2)>electricfield2(thresh,2)
        electricfield2(i,6)=electricfield2(thresh,6);
        electricfield2(i,7)=1/((1/electricfield2(i,6))+(1/cox(par)));
    end
end
bench_cap(par)=1/((1/electricfield2(thresh,6))+(1/cox(par))+(1/cox(par)))
figure(2)
plot(electricfield2(:,5),electricfield2(:,7))
%% CV curve for TSV pair
electricfield1(:,8)=smooth(electricfield1(:,7),10);
electricfield2(:,8)=smooth(electricfield2(:,7),10);
for i=1:range
tsv_total_cap(i,1)=1/((1/electricfield1(i,7))+(1/electricfield2(i,7)));
tsv_total_voltage(i,1)=electricfield1(i,5)-electricfield2(i,5);
end
tsv_total_cap(:,2)=smooth(tsv_total_cap(:,1),10);
figure(3)
plot(tsv_total_voltage,tsv_total_cap)
%% symerical transform
for i=1:range
    if tsv_total_voltage(i)<0
    tsv_hf(i,1)=tsv_total_voltage(i);
    tsv_hf(i,2)=electricfield1(i,7);%take not smoothed data of tsv1 
    tsv_hf(i,3)=electricfield2(i,8);%take smoothed data of tsv2
    tsv_hf(i,4)=tsv_total_cap(i,2);
    end
end
mat=size(tsv_hf)
for i=mat(1)+1:mat(1)*2
    tsv_hf(i,1)=-tsv_hf(2*mat(1)+1-i,1);
    tsv_hf(i,3)=tsv_hf(2*mat(1)+1-i,2);%exchange tsv1 for tsv2 for a continus curve
    tsv_hf(i,2)=tsv_hf(2*mat(1)+1-i,3);
    tsv_hf(i,4)=tsv_hf(2*mat(1)+1-i,4);
end
plot(tsv_hf(:,1),tsv_hf(:,2),tsv_hf(:,1),tsv_hf(:,3),tsv_hf(:,1),tsv_hf(:,4))

para_radius_vol(:,par)=tsv_hf(:,1)
para_radius_cap(:,par)=tsv_hf(:,4)/bench_cap(par)
clear tsv_hf
end
plot(para_radius_vol(:,1),para_radius_cap(:,1),para_radius_vol(:,2),para_radius_cap(:,2),para_radius_vol(:,3),para_radius_cap(:,3),para_radius_vol(:,4),para_radius_cap(:,4),para_radius_vol(:,5),para_radius_cap(:,5))


    

