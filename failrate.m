br96 = brRTS96;
mpc19 = loadcase('RTS_GMLC.m');
gen = mpc19.gen;
[rownum,~] = size(gen);
lamdagen = zeros(rownum,1);
for i = 1:rownum
    pmax = gen(i,9);
    switch pmax
        case 20
            lamdagen(i) = 19.4467;
        case 76
            lamdagen(i) = 4.46940;
        case 100
            lamdagen(i) = 7.3;
        case 197
            lamdagen(i) = 9.2211;
        case 12
            lamdagen(i) = 2.9796;
        case 155
            lamdagen(i) = 9.125;
        case 400
            lamdagen(i) = 7.9636;
        case {50,55}
            lamdagen(i) = 4.4242;
        case {350,355}
            lamdagen(i) = 7.6174;
        otherwise
            lamdagen(i) = 0.01;
    end
end

lamdabus = br96(:,10);

save("lamda.mat", "lamdagen", "lamdabus");
