function nkzfs8, lambda
x = lambda
n = SQRT( 1 + 1.62693651*x^2/(x^2-0.010880863) + 0.24369876*x^2/(x^2-0.0494207753) + 1.62007141*x^2/(x^2-131.009163) )
return, n
end
