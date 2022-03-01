x = [7000 1000 17 0];
count = 1;
distr = [];
for id_1 = 1 : length(x)
    for id_2 = 1 : x(id_1)
        distr(count) = id_1 - 1;
        count = count + 1;
    end
end
pd = fitdist(distr','gamma')