function ctx_coord = read_cortical_data( PRM)

% get information about the cortical stimulation for a session
% information contains which quadrant, session id, number of electrodes
% used for stimulation, threshold and the response to this stimulation.


ctx_coord(1).id = PRM.ID;

for i=1:2,
    ctx_coord(i).qd = [];
    ctx_coord(i).x = [];
    ctx_coord(i).y = [];
    ctx_coord(i).x = [];
    ctx_coord(i).thr = [];
    ctx_coord(i).eff = [];
end;
if isfield(  PRM, 'Cortex1'), % there are two electrodes
    if PRM.Cortex1.Flag
        ctx_coord(1).qd = PRM.Cortex1.Quad;
        ctx_coord(1).x = PRM.Cortex1.Y;
        ctx_coord(1).y = PRM.Cortex1.X;
        ctx_coord(1).thr = PRM.CTXmap(1).Thr;
        ctx_coord(1).eff = quantify_effect( PRM.CTXmap(1).Resp);
    end
    if PRM.Cortex2.Flag,
        ctx_coord(2).qd = PRM.Cortex2.Quad;
        ctx_coord(2).x = PRM.Cortex2.Y;
        ctx_coord(2).y = PRM.Cortex2.X;
        ctx_coord(2).thr = PRM.CTXmap(2).Thr;
        ctx_coord(2).eff = quantify_effect( PRM.CTXmap(2).Resp);
        
    else
        return;
    end;
else
    if PRM.Cortex.Flag,
        disp('Single electrode');
        ctx_coord(1).qd = PRM.Cortex.Quad;
        ctx_coord(1).x = PRM.Cortex.Y;
        ctx_coord(1).y = PRM.Cortex.X;
        ctx_coord(1).thr = PRM.CTXmap(1).Thr;
        ctx_coord(1).eff = quantify_effect( PRM.CTXmap(1).Resp);
    else
        return;
    end
end