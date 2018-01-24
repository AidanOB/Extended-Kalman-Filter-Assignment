function [n, assoc, ind, indMap] = dataAssociation(xxO, yyO, xxG, yyG)
    
    [C, indMap] = min(pdist2(single([xxG, yyG]), single([xxO', yyO'])));
    large = find(C > 0.5);
    assoc = indMap;
    assoc(large) = 0;
    indMap(large) = [];
    ind = find(assoc);
    n = length(indMap);
    assoc = assoc';
    indMap = indMap';
    ind = ind';
  
end