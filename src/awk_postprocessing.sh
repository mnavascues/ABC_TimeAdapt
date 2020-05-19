#!/bin/sh

awk '
BEGIN { srand(seed.awk);
  missing[10]=0.2;
  missing[11]=0.3;
  missing[12]=0.0005;
  missing[13]=0.0029;
  missing[14]=0.001;
  missing[15]=0.0082;
  missing[16]=0.02;
  missing[17]=0.00002;
  missing[18]=0.002;
  missing[19]=0.0002;
  missing[20]=0.0012;
  missing[21]=0.005;
}
/#/  { print $0 }
!/#/ {
  for (i=10; i<=11; i++) {
    split($i,haplo,"|",seps);
    x=haplo[int(1+2*rand())];
    $i = x"/"x
  };
  for (i=10; i<=21; i++) {
    if (rand()<missing[i]){ $i = "./."}
  };
  print $0
}
' results/genotypes.vcf > results/genotypes2.vcf