function [newformulas] = quick13C_KL_1(formulas,NeutralMass,showProgress)  
%function [newformulas] = quick13C_KL_1(formulas,NeutralMass,showProgress)  
% attempt at finding 13C pairs in MAB0905 dataset
%   Written by E. Kujawinski-Behn 
%   Edited by M. Bhatia 06.30.09 to look for peaks at every mass so that
%   the code could be used for multiply charged molecules 
% Note that this will require the MAT file atom_masses.mat
% Krista Longnecker 9/1/2011 - changed the error calculation to spit out the absolute
% error to be consistent with findformula
%%contact information:
%
% Krista Longnecker
% Department of Marine Chemistry and Geochemistry
% Woods Hole Oceanographic Institution
% klongnecker@whoi.edu

NomIUPAC = round(NeutralMass);
range = [min(NomIUPAC):1:max(NomIUPAC)];
% EvenOdd = rem(range(1),2);  % MB commented out lines 10-15
% if EvenOdd == 0
%     a = 1;
% else
%     a = 2;
% end
load atom_masses.mat;
elements = [C; H; O; N; C13; S ; P ; Na];
t=0; % MB set counter

% for i=a:2:(length(range)-1)   % MB commented out line 20, inserted line 21
for i=1:(length(range)-1)
%     i
%     range(i)
    tempParent = find(round(NeutralMass)==range(i));
    tempDaughter = find(round(NeutralMass)==range(i+1));
    if isempty(tempParent) == 0 & isempty(tempDaughter)==0
        for x = 1:length(tempParent)
            diff = NeutralMass(tempDaughter) - NeutralMass(tempParent(x));
            iso = find(diff<1.0035 & diff>1.003);
%             diff(iso);
            if length(iso)>1
                tp = find(min(abs(1.0034-diff(iso))));
                iso = iso(tp);
            end
            if isempty(iso)==0 % MB instered counter (lines 30-33)
                t=t+1;
             c13pairs(t,:) = [tempParent(x) tempDaughter(iso)];  % MB added to keep track of indices of C13 pairs
            end
            if isempty(iso) == 0 & formulas(tempParent(x),1)>0
                ParentForm = [formulas(tempParent(x),1:8),0];
                newDForm = [ParentForm(1)-1, ParentForm(2:4), 1, ParentForm(6:8)];  
                oldDErr = formulas(tempDaughter(iso),9);
                newDMass = newDForm*elements;
                newDErr = abs((NeutralMass(tempDaughter(iso))-newDMass)/newDMass*1e6);
                if abs(newDErr) < 0.5
                    formulas(tempDaughter(iso),:)=[newDForm,newDErr];
                    if showProgress %KL put this here 9/1/2011 bc don't want to see this
                        fprintf('Parent #=%1.0f Daughter # = %1.0f\n',tempParent(x),tempDaughter(iso))
                        fprintf('Old Error = %f New Error = %f\n', oldDErr,newDErr)
                    end
                end
            end
            clear diff iso
        end
    end
    clear tempParent tempDaughter   
end

% totC13pairs = t   % MB included new output total C13 pairs found 
newformulas = formulas;   