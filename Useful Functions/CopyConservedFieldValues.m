

function [Struct, bSetSomething, bModifiedSomething] = CopyConservedFieldValues(Struct, StructCopyFromEdit, bAllowNewFieldCreation, StructCopyFromOrig)
%function [Struct, bSetSomething, bModifiedSomething] = CopyConservedFieldValues(Struct, StructCopyFromEdit, bAllowNewFieldCreation, StructCopyFromOrig)
%
%Inputs and outputs are all structures.
%This function recurses the structure fields.  Whenever StructCopyFromEdit has an identical
%field as Struct, the value of this field is copied to Struct.  
%
%if bAllowNewFieldCreation==0,  No fields are added or removed.
%if bAllowNewFieldCreation==1,  fields present in StructCopyFromEdit that are not present in Struct will be created
%
%if StructCopyFromOrig has a field that no longer exists in StrcutCopyFromEdit
%   then this field will be removed from Struct.  To disable this feature, just set StructCopyFromOrig=[]
%
%
%
%NOTE: structures are NEVER copied, only terminal fields.
%
%Note: if any structures (or substructures) are vectors, then, two possible cases:
%  1. if StructCopyFromEdit is a single item, then the value of the single item in StructCopyFromEdit
%     will be copied to every element in the Vector of Struct
%  2. if StructCopyFromEdit is more than one item, then each element of
%     StructCopyFromEdit is copied to Struct.  
%     (So if Struct is shorter, then it's length will be increased)
%     (and if Struct is longer, then it's last few elements be unchanged)

%Copyright (c) 2004 by Michael Tadross. All rights reserved. 
%You may not without express written permission:
%  1. Redistribute this source code, fraction of source code, or modified version of source code.
%  2. Copy or use any source code, fraction of source code, or modified version of source code for any commercial purpose;
%  3. Remove or change any copyright, trademark or other intellectual property right notices contained in the original material copies or printed off from these materials. 
%You may:
%  1. Use this code for non-commercial personal or educational research use only.
%  2. Collaborate with the creator of this source code to develop added functionality.  All such collaborations will be subject to the copyright disclaimer herein.
%
%Code is provided "as is."  Use of this code is at your sole risk
% $Revision: 1.0 $  $Date: 2004/04/08 05:29:34 $


if nargout==3
    bNeedModifiedInfo = 1;
else
    bNeedModifiedInfo = 0; %for speed
end
[Struct, bSetSomething, bModifiedSomething] = CopyConservedFieldValues2(Struct, StructCopyFromEdit, bAllowNewFieldCreation, StructCopyFromOrig, bNeedModifiedInfo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Struct, bSetSomething, bModifiedSomething] = CopyConservedFieldValues2(Struct, StructCopyFromEdit, bAllowNewFieldCreation, StructCopyFromOrig, bNeedModifiedInfo)

%set flag if something was modified
bSetSomething = 0;
bModifiedSomething = 0;

%The only special case is if neither is a struct.  Then just copy the value.
if ~isstruct(Struct) && ~isstruct(StructCopyFromEdit)
    %Note: we only get here if all the cel stuff is OK . . see code prior to recursive call.
    bSetSomething = 1;
    if bNeedModifiedInfo && ~StructuresAreIdentical(Struct, StructCopyFromEdit)
        bModifiedSomething = 1;
    end
    Struct = StructCopyFromEdit;
    return;
end

%make sure StructCopyFromOrig is a structure 
if ~isstruct(StructCopyFromOrig)
    StructCopyFromOrig = struct; %make it a struct with no fields
end


%Otherwise we recurse:
FieldName = unique([GetFieldNames(StructCopyFromEdit); GetFieldNames(StructCopyFromOrig)]);  %get the field names in the CopyFrom struct. (unique from both the Edit and Orig)
for k=1:length(FieldName)
    if ~isfield(StructCopyFromEdit, FieldName{k}) 
        %this field was either deleted or renamed! . . . 
        %remove this field, it was deleted in StructCopyFromEdit, yet it exists in StructCopyFromOrig
        if isfield(Struct, FieldName{k})
            Struct = rmfield(Struct, FieldName{k});
            bSetSomething = 1;
            bModifiedSomething = 1;
        end
        
    elseif bAllowNewFieldCreation ||  isfield(Struct, FieldName{k})
        %the field exists in StructCopyFromEdit, and it either exists in Struct
        %OR we are allowing the creation of new fields
        
        %NOTE: length(Struct) may be == 0, so this is actually necessary
        if length(StructCopyFromEdit)==1 &&  length(Struct)>1
            sizeStruct = size(Struct);
            nMax = length(Struct(:));  %copy single item to every item of Struct
        else
            sizeStruct = size(StructCopyFromEdit);
            nMax = length(StructCopyFromEdit(:)); %only go up to length of StructCopyFromEdit, grow Struct if necessary
        end

        for n=1:nMax  %linear indexing
            %get the CopyFrom field value
            if length(StructCopyFromEdit) == 1
                j = 1;
            else
                j=n;
            end
            %edit value always exists
            FromValueEdit = StructCopyFromEdit(j).(FieldName{k});
            %orig value may not exist . .  so do a try/catch
            try
                FromValueOrig = StructCopyFromOrig(j).(FieldName{k});
            catch
                FromValueOrig = []; %this is equivalent to not having any info
            end

            %now, we have some rules  . .
            %There's a lot going on here.
            %possible values of "FromValueEdit"
            %   1.  [12 23 24]  or  'haha'              --  field exists in all,          with single value
            %   2.  {[12 23 24]}  or  {'haha'}          --  field does not exist in all,  single value in those that have it
            %   3.  {[1 34]  [3]}  or {'one' 'two'}     --  field exists in all,          but different values
            %   4.  {{[1 34]  [3]}}  or {{'one' 'two'}} --  field does not exist in all,  different value in those that have it

            %we only copy a value from FromValueEdit to ToValue if:
            %Case1 --  field exists in Struct OR  bAllowNewFieldCreation
            %Case2 --  only if field already exists in Struct,  
            %Case3/4 -- never!
            [FromValueEdit, nTypeFromValueEdit] = GetCellValue(FromValueEdit);
            if any(nTypeFromValueEdit==[3 4])
                continue;  %don't copy this one
            end
        
            %get the CopyTo value
            if ~isfield(Struct, FieldName{k})  ||  length(Struct(:))<n  
                %field doesn't exist in Struct . . OR Struct is 0x0 (or just shorter than StructCopyFromEdit)
                %only continue if bAllowNewFieldCreation AND nTypeFromValueEdit==1
                if nTypeFromValueEdit==1 && bAllowNewFieldCreation
                    if isstruct(FromValueEdit) && ~isempty(FromValueEdit)
                        ToValue = struct; %a 1x1 struct with no fields
                    else
                        ToValue = []; %could either end up being a terminal field or a 0x0 struct with fields
                    end
                else
                    %skip this time
                    continue;
                end
            else
                %field exists in Struct . . no need to create a new field
                ToValue = Struct(n).(FieldName{k});
            end
            
            %recurse on this field and set it.
            [NewStructValue, bSet, bModified] =  CopyConservedFieldValues2(ToValue, FromValueEdit, bAllowNewFieldCreation, FromValueOrig, bNeedModifiedInfo);
            if bSet
                %only set it if it has at least one terminal field somewhere . . .
                %determine indices
                switch length(sizeStruct)
                    case 1
                        %vector
                        Struct(n).(FieldName{k}) = NewStructValue;
                    case 2
                        %matrix
                        [I1, I2] = ind2sub(sizeStruct, n);
                        Struct(I1,I2).(FieldName{k}) = NewStructValue;
                    case 3
                        %3-D matrix
                        [I1, I2, I3] = ind2sub(sizeStruct, n);
                        Struct(I1,I2,I3).(FieldName{k}) = NewStructValue;
                    otherwise
                        error('Error in CopyConservedFieldValues.  Higher than 3-D struct matrix not yet supported');
                end
            end
            bSetSomething = bSetSomething || bSet;
            bModifiedSomething = bModifiedSomething || bModified;
        end
    end
end

%Finally, make the order of fields as much like StructCopyFromEdit as possible
%NOTE: Not all the fields in Struct are guaranteed to be in StructCopyFromEdit (nor vice versa)
%      so we will only change the order of fields that are in both
FromEditFields = GetFieldNames(StructCopyFromEdit);
StructFields   = GetFieldNames(Struct);
if ~isempty(StructFields)
    Order = zeros(length(StructFields), 1); %start with all zeros
    for k=1:length(StructFields)
        n = find(strcmp(StructFields{k}, FromEditFields));
        if ~isempty(n)
            %the k'th field of Struct exists in FromEditFields
            %store it's field position
            Order(k) = n;
        end
    end
    %This code is confusing, but the end result is it gives us the
    %correct order of the fields in Struct
    n = find(Order);           %indices of field names that exist in both
    [dum, I] = sort(Order(n)); %sort order to get these names into the same order as in FromEditFields
    [dum, II] = sort(I);       %inverse map
    Order = [1:length(Order)]; %start with 1 2 3 4 . .  n
    Order(n) = n(II);          %re-order the ones that are in both
    [dum, Order] = sort(Order); % . ..  yup, it's tricky
    %now order the fields
    Struct = orderfields(Struct, Order);
end

%That's all.



function Names = GetFieldNames(Struct)
%function Names = GetFieldNames(Struct)
%
%very similar to the function fieldnames(Struct) except
%it also supports Struct = []

%Copyright (c) 2004 by Michael Tadross. All rights reserved. 
%You may not without express written permission:
%  1. Redistribute this source code, fraction of source code, or modified version of source code.
%  2. Copy or use any source code, fraction of source code, or modified version of source code for any commercial purpose;
%  3. Remove or change any copyright, trademark or other intellectual property right notices contained in the original material copies or printed off from these materials. 
%You may:
%  1. Use this code for non-commercial personal or educational research use only.
%  2. Collaborate with the creator of this source code to develop added functionality.  All such collaborations will be subject to the copyright disclaimer herein.
%
%Code is provided "as is."  Use of this code is at your sole risk
% $Revision: 1.0 $  $Date: 2004/04/08 05:29:34 $

if ~isstruct(Struct)
    Names = {};
else
    Names = fieldnames(Struct);
end




function [Value, nType] = GetCellValue(Value)
%function [Value, nType] = GetCellValue(Value)
%
%There's a lot going on here.
%possible forms of Value
%   1.  [12 23 24]  or  'haha'              --  field exists in all,          with single value
%   2.  {[12 23 24]}  or  {'haha'}          --  field does not exist in all,  single value in those that have it
%   3.  {[1 34]  [3]}  or {'one' 'two'}     --  field exists in all,          but different values
%   4.  {{[1 34]  [3]}}  or {{'one' 'two'}} --  field does not exist in all,  different value in those that have it
%
%Value is like 1/3 . . . 
%nType is 1,2,3,4 as above

bExist = 1;
if iscell(Value)
    %case 2,3, or 4
    if length(Value)==1
        %case 2 or 4
        Value = Value{1};
        if iscell(Value)
            nType = 4;
        else
            nType = 2;
        end
    else
        nType = 3;
    end
else
    nType = 1;
end



