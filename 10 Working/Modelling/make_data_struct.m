
function [ data_struct ] = make_data_struct( varargin )

% make_data_struct(data, FS, spl)
% make_data_struct(data, FS, spl, audiogram_struct, 'prescription_type')


if length(varargin) < 2
    disp('Error, not enough inputs');
    data_struct = struct([]);
    return;
elseif (length(varargin) == 5) && isnumeric(varargin{1})
    x = varargin{1};
    FS = varargin{2};
    spl = varargin{3};
    calc_details = varargin{4};
    synaptopathy = varargin{5};
    audi_struct = [];
    pres_type = [];
    
elseif (length(varargin) == 7) && isnumeric(varargin{1}) && isstruct(varargin{4}) && ischar(varargin{5})
    x = varargin{1};
    FS = varargin{2};
    spl = varargin{3};
    audi_struct = varargin{4};
    pres_type = varargin{5};
    calc_details = varargin{6};
    synaptopathy = varargin{7};
else
    disp('Check your inputs!')
    data_struct = struct([]);
    return;
end

data_struct.data = set_spl(x(:)',spl);
data_struct.FS = FS;
data_struct.SPL = spl;
data_struct.SYNAPTOPATHY = synaptopathy;

if strncmpi(calc_details, 'detailed', 6)
    %data_struct.fundamental = get_fund(x,FS);
    %data_struct.formants = get_form(x,FS);
    [data_struct.F, data_struct.DATA, data_struct.PHASE ] = quickfft( data_struct.data, data_struct.FS );
    %data_struct.approx_formants = round(data_struct.formants/data_struct.fundamental)*data_struct.fundamental;
    data_struct.calc_details = 'detailed';
else
    data_struct.calc_details = 'simple';
end

if isstruct(audi_struct) && ischar(pres_type)
    data_struct = ampl_pres( data_struct, audi_struct, pres_type, calc_details );
end


