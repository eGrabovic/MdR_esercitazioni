function updatePlot_localPOE(Glocal, q, transforms)
%
%
%

Glocalnum = cell(1, length(q));
[Glocalnum{:}] = Glocal(q);
Glocalnum = cellfun(@full, Glocalnum, 'UniformOutput',false);

for i = 1:length(q)
   
    transforms{i}.Matrix = Glocalnum{i};
    
end



end