
%clc
diary('log.txt');
for monk = {'vega'} %{'vega', 'chalva', 'darma'}
    runscript4emgdat(char(monk))
end
diary('off');