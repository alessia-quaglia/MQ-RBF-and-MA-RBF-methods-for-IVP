function table_latex(met_name, N, met_err, met_ord)
%
% Usage:        table(met_name, N, met_err, met_ord)
% Purpose:      it generates a LaTeX table showing the global errors and 
%               the convergence orders versus different values of N
%               for the analyzed methods
% Input: met_name = each element is a string containing the names of the methods
%               N = vector containing the number of steps considered
%         met_err = each element is a vector of global errors
%         met_ord = each element is a vector of convergence orders
% Output:       this function does not return values; it prints the LaTeX 
%               code directly to the command window
%
fprintf('\\begin{table}[ht]\n\\centering\n');
fprintf('\\begin{tabular}{llcc}\n\\hline\n');
fprintf('Method & $N$ & Global error & Order \\\\ \\hline\n');
for m = 1:length(met_name)
    current_err  = met_err{m};
    current_ord  = met_ord{m};
    for i = 1:length(N)
        if i == 1
           name_str = met_name{m};
        else
           name_str = '';
        end
        if i == 1
           ord_str = '-';
        else
           ord_str = sprintf('%.4f', current_ord(i));
        end
        fprintf('%-20s & %d & %.2E & %s \\\\ \n', ...
        name_str, N(i), current_err(i), ord_str);
    end
    if m < length(met_name)
       fprintf('\\hline \n');
    end
end
fprintf('\\hline\n\\end{tabular}\n');


