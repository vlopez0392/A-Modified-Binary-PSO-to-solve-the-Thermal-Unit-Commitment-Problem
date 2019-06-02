function [I_C_SORT_EXS,IDX_EXS,I_C_SORT_INS,IDX_INS] = AFLC(ai,bi,ci,PGI_MAX)
%Returns a sorted array in a descending fashion according to AFLC criteria.
%Returns the unit (idx) in the sorted array 

    I_C_sort=((2*(ai.*ci).^(1/2)+bi)./PGI_MAX);

    [I_C_SORT_EXS,IDX_EXS]=sort(I_C_sort,'descend');

    [I_C_SORT_INS,IDX_INS]=sort(I_C_sort,'ascend');
end

