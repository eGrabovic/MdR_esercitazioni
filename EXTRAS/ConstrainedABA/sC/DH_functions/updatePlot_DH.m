function updatePlot_DH(Tj_fun, q_num, T0, Te, transforms)
%
%
%

Tj_num = cell(length(q_num), 1);
[Tj_num{:}] = Tj_fun(q_num);

T = T0;
for j = 1:length(q_num)
    T = T*full(Tj_num{j});
    transforms{j}.Matrix = T;
end
transforms{j+1}.Matrix = Te;

end