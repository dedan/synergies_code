
%clc
diary('log.txt');
for monk = {'vega', 'chalva', 'darma'}
    runscript4emgdat(char(monk))
end
diary('off');