awk_args <- paste("-v OFS='\t' '",
                  "BEGIN { srand(",seed.awk,");")
for (samp in seq_len(total_sample_size)){
  awk_args <- paste(awk_args,
                    "missing[",9+samp,"]=",missing[samp],";")
}
awk_args <- paste(awk_args,
                  "}",
                  "/#/ { print $0 }",
                  "!/#/ {",
                  "  for (i=",9+total_sample_size-na,"; i<=",9+total_sample_size,"; i++) {",
                  "    split($i,haplo,\"|\",seps);",
                  "    x=haplo[int(1+2*rand())];",
                  "    $i = x\"/\"x",
                  "  };",
                  "  for (i=10; i<=",9+total_sample_size,"; i++) {",
                  "    if (rand()<missing[i]){ $i = \"./.\"}",
                  "  };",
                  "  print $0",
                  "}",
                  "' results/genotypes.vcf > results/genotypes2.vcf")
