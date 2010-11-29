function [edname, id] = extract_edname( monk,  subss, sessname, findx , e2add)

id      = -1;
i       = 1;
edname  = '';
while i < length(subss),
    if strcmp( char(subss(i).Session), sessname ), % this is the right session,
        if ismember( findx, subss(i).Files(1):subss(i).Files(end)),
            ID = subss(i).ID;
            id = subss(i).ID;
            subID = subss(i).SubSess;
            ednum = find(  subss(i).Files(1):subss(i).Files(end) == findx);
            if ID < 10,
                str1 = ['0' num2str(ID)];
            else
                str1 = num2str(ID);
            end
            if subID < 10,
                str2 = ['0' num2str(subID)];
            else
                str2 = num2str(subID);
            end
            edname = [monk str1 str2 'e' e2add '.' num2str(ednum) '.mat'];
            return;
        end
    end
    i = i+1;
end
disp( 'EDname not found' );