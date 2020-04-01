function result = ncVarInFile(file,var)

vars = nclistvars(file);

result = sum(ismember(vars,var));