function idx = val2idx(val_arr,val_inp)
%VAL2IDX   Finds the index into an array in which the indexed array
%          value is closest to an input value.
%
%          IDX = VAL2IDX(VAL_ARR,VAL_INP) Given an array of values,
%          VAL_ARR, and a scalar input value, VAL_INP, the function
%          returns the index into the array in which the indexed array
%          value is closest to the input value.
%
%          NOTES:  None.
%
%         20-Mar-2024 * Mack Gardner-Morse
%

%#######################################################################
%
% Check Inputs
%
if (nargin<2)
  error(' *** ERROR in val2idx:  Two inputs are required!');
end
%
% Find Index to Value Closest to the Input Value
%
dv = val_arr-val_inp;
[~,idx] = min(abs(dv));
%
return
