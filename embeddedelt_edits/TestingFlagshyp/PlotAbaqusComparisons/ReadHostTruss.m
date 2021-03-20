%Read Abaqus Files
%Read One Host with One Embedded Truss
function [AbqEHost, AbqETruss, AbqE] = ReadHostTruss()
file = fopen('C:/Users/Valerie/Documents/GitHub/flagshyp/embeddedelt_edits/TestingFlagshyp/OneHostOneTrussResults.txt','r');
for i=1:4
    tline = fgetl(file);
end

formatSpec = '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f';
sizeA = [20 inf ];
% A = fscanf(file,[repmat(formatSpec,18)],sizeA);
A = fscanf(file,formatSpec,sizeA); A=A';
fclose(file);

           


AbqEHost.time = A(:,1);
AbqEHost.Acceleration = A(:,2);
AbqETruss.Acceleration = A(:,3);
AbqEHost.Strain = A(:,[9 11]);
AbqETruss.Strain = A(:,10);
AbqEHost.Force = A(:,12);
AbqETruss.Force = A(:,13);
AbqEHost.Stress = A(:,[14 16]);
AbqETruss.Stress = A(:,15);
AbqEHost.Displacement = A(:,[17 19]);
AbqETruss.Displacement = A(:,[18 20]);

AbqE.IE = A(:,4);
AbqE.KE = A(:,5);
AbqE.SE = A(:,6);
AbqE.WK = A(:,7);
AbqE.ETOTAL = A(:,8);

end