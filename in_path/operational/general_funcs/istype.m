function istrue = istype(datatype, field)

% Returns istrue = true if the field is true. If no output is specified,
% then an error message is generated.

% datatype is a struct, field is a string



% Check if the field exists. If it does not, assume that it is false.
istrue = isfield(datatype, field);

% Now check (and overwrite) the contents of the field.
if istrue
    eval(['foi = datatype.' field])
    istrue = foi;
end

% So now deal with the case where if there is no output argument
if nargout < 1
    if ~istrue
        error('The field specified does not exist or is false')
    end
    clear istrue
end