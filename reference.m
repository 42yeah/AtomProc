
clc;clear;
file = 'packing-PBT-101_227repeated-units-PBT-Cu-big-voids--longFilling-1800atm-packing-800-Cu-0K-2-5ps-2-No-Cu-H_O_1.data'; % file path 
Natoms           = 114408   ; % atom total number  
Bond_number      =  122544  ;
Angle_number     =  147024  ;
Dihedrals_number =  171504  ;
Impropers_number =  16344 ;
delete_type=3;
a=929.8119352918654-669.5345893141348;
b=714.2359190958651-453.95857311814075;
c=567.5667939802589-341.2362060197376;   %Box size

%---------------------------------------------------------------%
try
    dump = fopen(file,'r'); 
catch
    error('Dumpfile not found!');
end                    
%---------------------------------------------------------------%
atom_data_molecule      = zeros(Natoms,10);          % atom-ID molecule-ID atom-type q x y z
atom_bond      = zeros(Bond_number,4);     % id type atom1 atom2
atom_angle     = zeros(Angle_number,5);    % id type atom1 atom2 atom3
atom_dihedrals = zeros(Dihedrals_number,6);% id type atom1 atom2 atom3 atom4 atom5
atom_impropers = zeros(Impropers_number,6);% id type atom1 atom2 atom3 atom4
atom_masses    = zeros(Natoms,2);
atom_masses(:,1)=1:Natoms;
atom_masses(:,2)=12;
%---------------------------------------------------------------%
while feof(dump) == 0
  id = fgetl(dump);
  %----------------------------------------------%
  % notice format 
   if (strncmpi(id,'Atoms',numel('Atoms')))
  %----------------------------------------------%
      if (~isempty(id))
         disp(fgetl(dump));
         for j = 1 : 1: Natoms
            atom_data_molecule(j,:) = sscanf(fgetl(dump),'%d %d %d %f %f %f %f %d %d %d',[1 10]);
         end
      end
  %----------------------------------------------%
  % notice format 
      elseif (strncmpi(id,'Bonds',numel('Bonds')))
  %----------------------------------------------%
      if (~isempty(id))
         disp(fgetl(dump));
         for j = 1 : 1: Bond_number
            atom_bond(j,:) = sscanf(fgetl(dump),'%d %d %f %f',[1 4]);
         end
      end
  %----------------------------------------------%
  % notice format 
      elseif (strncmpi(id,'Angles',numel(' Angles')))
  %----------------------------------------------%
      if (~isempty(id))
         disp(fgetl(dump));
         for j = 1 : 1: Angle_number
            atom_angle(j,:) = sscanf(fgetl(dump),'%d %d %f %f %f',[1 5]);
         end
      end
  %----------------------------------------------%
  % notice format 
      elseif (strncmpi(id,'Impropers',numel('Impropers')))
  %----------------------------------------------%
      if (~isempty(id))
         disp(fgetl(dump));
         for j = 1 : 1: Impropers_number
            atom_impropers(j,:) = sscanf(fgetl(dump),'%d %d %f %f %f %f',[1 6]);
         end
      end
  %----------------------------------------------%
  % notice format 
      elseif (strncmpi(id,'Dihedrals',numel('Dihedrals')))
  %----------------------------------------------%
      if (~isempty(id))
         disp(fgetl(dump));
         for j = 1 : 1: Dihedrals_number
            atom_dihedrals(j,:) = sscanf(fgetl(dump),'%d %d %f %f %f',[1 6]);
         end
      end
   end     
end

position=find(atom_data_molecule(:,3)==delete_type);   
ring_carbon=atom_data_molecule(position,:);

%%%%%%%%%%%%%%%%%%% 将在一个苯环上的芳香碳原子放在连续六行 %%%%%%%%%%%%%
ring_carbon_new=zeros(length(ring_carbon),13);
mol_number=max(ring_carbon(:,2));
for i=1:mol_number;    % i 表示分子编号
    i
   mol=find(ring_carbon(:,2)==i);   
   ring_carbon_mol=ring_carbon(mol,:);
    for j=1:length(ring_carbon_mol)
        atom_id=zeros(3,1);
        atom_1=find(atom_bond(:,3)==ring_carbon_mol(j,1));
        atom_2=find(atom_bond(:,4)==ring_carbon_mol(j,1));
        atom_1_id=atom_bond(atom_1,4);
        atom_2_id=atom_bond(atom_2,3);
        
        if (length(atom_1_id)+length(atom_2_id))==2
        atom_id(1:2)=[atom_1_id; atom_2_id];
        else
        atom_id=[atom_1_id; atom_2_id];
        end
        
        ring_carbon_new((i-1)*length(ring_carbon_mol)+j,:)=[ring_carbon_mol(j,:),atom_id'];
    end
end


ring_carbon_find=zeros(length(ring_carbon_new),19);
for k=1:length(ring_carbon_new)
  k;
if ring_carbon_new(k,13)==0;
    ring_carbon_find(k,1:13)=ring_carbon_new(k,1:13);
    ring_carbon_find(k,14)=ring_carbon_new(k,1);
    atom_1_row=find(ring_carbon_new(:,1)==ring_carbon_new(k,11));
    atom_2_row=find(ring_carbon_new(:,1)==ring_carbon_new(k,12));
    if     ring_carbon_new(atom_1_row,13)>0;
      ring_carbon_find(k,15)=ring_carbon_new(atom_1_row,1);
    end
    if ring_carbon_new(atom_2_row,13)>0;
      ring_carbon_find(k,15)=ring_carbon_new(atom_2_row,1);
    end
    
    
    %if ring_carbon_find(k,15)>0
    %    k
    %else
    %   continue 
    %end
    
    atom_3_row=find(ring_carbon_new(:,1)==ring_carbon_find(k,15));
    sel_atom_3=ring_carbon_new(atom_3_row,11:13);
    sel_atom_3(sel_atom_3==ring_carbon_new(k,1))=[];
    
    if ismember(sel_atom_3(1),ring_carbon_new(:,1))==1
    sel_atom_3(2)=[];
    else
    sel_atom_3(1)=[];
    end
    ring_carbon_find(k,16)=sel_atom_3;
    
    sel_atom_4=ring_carbon_new(ring_carbon_new(:,1)==ring_carbon_find(k,16),11:12);
    sel_atom_4(sel_atom_4==ring_carbon_find(k,15))=[];
    ring_carbon_find(k,17)=sel_atom_4;

    sel_atom_5=ring_carbon_new(ring_carbon_new(:,1)==ring_carbon_find(k,17),11:12);
    sel_atom_5(sel_atom_5==ring_carbon_find(k,16))=[];
    ring_carbon_find(k,18)=sel_atom_5;

    sel_atom_6=ring_carbon_new(ring_carbon_new(:,1)==ring_carbon_find(k,18),11:13);
    sel_atom_6(sel_atom_6==ring_carbon_find(k,17))=[];
    if ismember(sel_atom_6(1),ring_carbon_new(:,1))==1
    sel_atom_6(2)=[];
    else
    sel_atom_6(1)=[];
    end
    ring_carbon_find(k,19)=sel_atom_6;

end
end
ring_carbon_find(ring_carbon_find(:,1)==0,:)=[];
save a.mat

order=ring_carbon_find(:,14:19);
ring=zeros(length(order),6);
    for i=1:length(order)
    ring(i,:)=sort(order(i,:));
    end
 ring=sortrows(ring);
 ring_sort=ring(1:4:length(ring),:)';
 ring_carbon_order=ring_carbon;
 for   i=1:length(ring_carbon)
     ring_carbon_order(i,:)=ring_carbon(ring_carbon(:,1)==ring_sort(i),:);
 end
    
 
 
 ring_carbon=ring_carbon_order;
%%%%%%%%%%%%%search atom_bond,统计每个芳香碳的键数%%%%%%%%%%%
N_bond=zeros(length(ring_carbon),1);
for i=1:length(ring_carbon)
   N_bond(i,1)=length(find(atom_bond(:,3:4)==ring_carbon(i,1)));
end
ring_carbon=[ring_carbon N_bond];

%%%%%%%%%%%%%%%%%%筛选要删除的原子%%%%%%%%%%%
ring_carbon_sort=sortrows(ring_carbon,2);
[m,n]=find(ring_carbon_sort(:,11)==3);   
ring_carbon_sort(m,:)=[];

%%%%%%%%%%%计算四个原子间距离%%%%%
distance=zeros(length(ring_carbon_sort),1);
for i=1:length(ring_carbon_sort)/4
   atom_1=ring_carbon_sort(4*i-3,5:7);
   for j=2:4
   deta_d=ring_carbon_sort(4*i-4+j,5:7)-atom_1;
    if abs(deta_d(1))>7
        deta_d(1)=abs(deta_d(1))-a;
    end
    if abs(deta_d(2))>7
        deta_d(2)=abs(deta_d(2))-b;
    end 
    if abs(deta_d(3))>7
        deta_d(3)=abs(deta_d(3))-c;
    end
   distance(4*i-4+j,1)=sqrt(deta_d*deta_d');
   end
end
max_distace=max(distance)     %输出最大值进行检验
ring_carbon_sort=[ring_carbon_sort distance];

%%%%%%%%%%%%%对距离最小的2个原子删除%%%%%%%%%
delete_atomid=zeros(length(ring_carbon_sort)/2,1);

for i=1:length(ring_carbon_sort)/4
    delete_atomid(2*i-1)=ring_carbon_sort(4*i-3,1);
    [p,q]=min((distance(4*i-2:4*i)));
    delete_atomid(2*i)=ring_carbon_sort(4*i-3+q,1);
end

%%%%%%%%%%%%%%%删除原子，bond,angle,dihedral,improper%%%%%%%%

atom_data_molecule_new=atom_data_molecule;

for i=1:length(delete_atomid)
atom_data_molecule_new(atom_data_molecule_new(:,1)==delete_atomid(i),:)=[];
end

atom_bond_new=atom_bond;
for  r=1:length(delete_atomid)
    
    [a1,a2]=find(atom_bond_new(:,3:4)==delete_atomid(r));
   atom_bond_new(a1,:)=[];
 
end

atom_angle_new=atom_angle;
for  r=1:length(delete_atomid)
    
    [a1,a2]=find(atom_angle_new(:,3:5)==delete_atomid(r));
   atom_angle_new(a1,:)=[];
 
end

atom_dihedrals_new=atom_dihedrals;
for  r=1:length(delete_atomid)
    
    [a1,a2]=find(atom_dihedrals_new(:,3:6)==delete_atomid(r));
   atom_dihedrals_new(a1,:)=[];
 
end


atom_impropers_new=atom_impropers;
for  r=1:length(delete_atomid)
    
    [a1,a2]=find(atom_impropers_new(:,3:6)==delete_atomid(r));
   atom_impropers_new(a1,:)=[];
 
end
%%%%%%%%%%%%%%%%%带后缀new的分别对应data文件的各个部分%%%%%%%%%%%%%






%%%%%%%%%%%%%% 读取原子Velocities %%%%%%%%%%%%%%%%%%%%
%{
}fidin = fopen('16-16.data'); % 文件名
atom_v_16=zeros(112000,4);
i=0;
n=0;
while ~feof(fidin)  % 判断是否为文件末尾
 tline = fgetl(fidin);
 i=i+1;
        if i>=112301
         a=str2num(tline);
         n=n+1;
        atom_v_16(n,:)= a;
        end
        if n==112000
            break
        end
            
end
fclose(fidin);
atom_v_16_new=atom_v_16;
atom_v_16_new(delete_atomid,:)=[];
%}
















































