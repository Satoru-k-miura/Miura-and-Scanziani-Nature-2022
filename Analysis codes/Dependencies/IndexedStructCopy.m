function T = IndexedStructCopy(S, Condition, FieldList)
if nargin == 2
   FieldList = fieldnames(S)';
end 
FieldList{2,1} = {};
T = struct(FieldList{:});
numDat = sum(Condition);
T(1:numDat) = S(Condition);
end
