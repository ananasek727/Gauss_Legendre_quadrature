% Comparison of real and calculated values of nodes and coefficients
% For n >= 3 & n <=101 & n mod 2 == 1
format longE
% Read the text file
data = fileread('Gauss-Legendre.txt');
% Split the data based on 'n = '
parts = strsplit(data, 'n = ');

coefficientsCell = cell(1,numel(parts));
nodesCell = cell(1,numel(parts));
% Loop through each part and extract coefficients and nodes (for odd n)
for i = 3:2:numel(parts)
    % Extract the value of n
    n = str2double(parts{i}(1));
        % Find the starting index of Coefficients and Nodes
        startCoefficients = strfind(parts{i}, 'Coefficients:') + length('Coefficients:');
        startNodes = strfind(parts{i}, 'Nodes:') + length('Nodes:');

        % Extract Coefficients and Nodes
        coefficients = str2num(parts{i}(startCoefficients:startNodes - length('Nodes:')- 1));
        nodes = str2num(parts{i}(startNodes:end));

        % Store coefficients and nodes in cell arrays
        coefficientsCell{i} = coefficients;
        nodesCell{i} = nodes;
end
% Calculating mean error for every n
error_coeff = zeros(1,numel(parts));
error_nodes = zeros(1, numel(parts));
for i = 3:2:numel(parts)
    [c, n] = coef_nodes(i);
    error_coeff(i) = mean(abs(c' - coefficientsCell{i}));
    % n is increasing but nodesCell are decreasing so we flip
    error_nodes(i) = mean(abs(flip(n',1) - nodesCell{i}));
end
% 
error_coeff_format = error_coeff(1, 3:2:numel(parts));
error_nodes_format = error_nodes(1,3:2:numel(parts));
n_format = 3:2:numel(parts);
T = table(n_format',error_coeff_format',error_nodes_format', 'VariableNames' ...
    ,{'n', 'error_coefficient', 'error_node'});
writetable(T,'Legendre_error.csv')