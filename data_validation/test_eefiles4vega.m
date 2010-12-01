function test_eefiles4vega( )

fid = fopen( 'report.txt', 'w');
hdir = ['i:\vega\data\'];
eedir = 'i:\vega\data\UpdateEdfilesMergedFinal\';

dirlist = dir([hdir '\v*07']);
dirlist = sortdirs( dirlist );

for i=1:length(dirlist),
    curdir = dirlist{i};
    if isdir([hdir curdir '\info']),
        disp(curdir);
        fprintf(fid,'%s\n', curdir);
        dinfo = load([hdir curdir '\info\' curdir '_param.mat']);
        ID = dinfo.DDFparam.ID;
        IDstr = padzero(ID);
        if isdir([hdir curdir '\MAT']),
            for subindx=1:length(dinfo.SESSparam.SubSess),
                files = dinfo.SESSparam.SubSess(subindx).Files;
                files = files(1):files(end);
                istr = padzero(subindx);
                for findx = files,
                    str = padzero2mat(findx);
                    bhvname = [hdir curdir '\MAT\' curdir str '_bhv.mat'];
                    ename = [hdir curdir '\edfiles\v' IDstr istr 'e.' num2str(findx-files(1)+1) '.mat'];
                    eename = [eedir 'v' IDstr istr 'ee.' num2str(findx-files(1)+1) '.mat'];
                    if ~exist(bhvname,'file'),
                        disp([bhvname '--> BHV file is missing']);
                        fprintf(fid,'%s --> BHV file is missing\n',bhvname);
                    elseif ~exist(ename,'file'),
                        disp([ename '--> E file is missing']);
                        fprintf(fid,'%s --> E file is missing\n',ename);
                    elseif ~exist(eename,'file'),
                        disp([eename '--> EE file is missing']);
                        fprintf(fid,'%s --> EE file is missing\n',eename);

                    else
                        bdat = load(bhvname, 'TrqFE*');
                        edat = load(ename);
                        eedat = load(eename);
                        T1 = length(bdat.TrqFE)/bdat.TrqFE_KHz;
                        T2 = max(edat.analog_time);
                        T3 = max(eedat.analog_time);
                        if T2 ~= T3 || abs(T1 - T2) > 100,
                            disp('Mismatch times');
                            fprintf(fid,'-----------------------------\n');
                            fprintf(fid,'%s\t%d\n', bhvname, T1);
                            fprintf(fid,'%s\t%d\n', ename, T2);
                            fprintf(fid,'%s\t%d\n', eename, T3);
                            fprintf(fid,'-----------------------------\n');

                        end
 
                        
                        
                    end;
                end
                
            end
        else
            disp('No mat files');
            fprintf(fid,'No mat files available\n');
        end
    end
end
            
fclose(fid);
            