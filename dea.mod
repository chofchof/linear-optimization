set INPUT;
set OUTPUT;

param N >= 1;
param u{1..N, INPUT}; # input values
param v{1..N, OUTPUT}; # output values

param K; # the DMU to assess
param eps > 0;

var x{i in INPUT} >= eps;
var y{j in OUTPUT} >= eps;

maximize objective:
  sum {j in OUTPUT} y[j] * v[K, j];

subject to this_dmu:
  sum {i in INPUT} x[i] * u[K, i] = 1;

subject to other_dmus{k in 1..N}:
  sum {j in OUTPUT} y[j] * v[k, j] <= sum {i in INPUT} x[i] * u[k, i];

data;

param eps := 0.00001;
param K := 1;
param N := 10;

set INPUT := realestate wages;
set OUTPUT := economics business mathematics;

param u: realestate wages :=
  1 72 81
  2 73 82
  3 70 59
  4 87 83
  5 53 64
  6 71 85
  7 65 68
  8 59 62
  9 134 186
  10 134 140;

param v: economics business mathematics :=
  1 77 73 78
  2 73 70 69
  3 72 67 80
  4 69 74 84
  5 57 65 65
  6 78 72 73
  7 81 71 69
  8 64 66 56
  9 150 168 172
  10 134 130 130;

end;
