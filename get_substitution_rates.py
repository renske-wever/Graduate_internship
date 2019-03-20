import sys

input_file = open(sys.argv[1], 'r')
output_file = open(sys.argv[2], 'w')

total = 0
AT = 0
AG = 0
AC = 0
TA = 0
TG = 0
TC = 0
GA = 0
GT = 0
GC = 0
CA = 0
CT = 0
CG = 0
counter = 0 
long_region = 0
double = 0

next(input_file)
for line in input_file:
    line = line.strip()
    line = line.split("\t")
    origin = line[2]
    alt = line[3]
    if len(origin) == 1 and len(alt) > 1 and "," not in alt:
        long_region += 1
        pass
    if len(origin) > 1:
        long_region += 1
        pass
    if len(origin) == 1 and "," in alt:
        double +=1
        alterations = alt.split(",")
        for i in alterations:
            if origin == "A":
                if i == "T":
                    AT += 1
                    total += 1
                if i == "C":
                    AC += 1
                    total += 1
                if i == "G":
                    AG += 1
                    total += 1
            if origin == "T":
                if i == "A":
                    TA += 1
                    total += 1
                if i == "C":
                    TC += 1
                    total += 1
                if i == "G":
                    TG += 1
                    total += 1
            if origin == "G":
                if i == "T":
                    GT += 1
                    total += 1
                if i == "C":
                    GC += 1
                    total += 1
                if i == "A":
                    GA += 1
                    total += 1
            if origin == "C":
                if i == "T":
                    CT += 1
                    total += 1
                if i == "A":
                    CA += 1
                    total += 1
                if i == "G":
                    CG += 1
                    total += 1
                   
    if len(origin) == 1 and len(alt) == 1 and "," not in alt:
        counter +=1
        if origin == "A":
            if alt == "T":
                AT += 1
                total += 1
            if alt == "C":
                AC += 1
                total += 1
            if alt == "G":
                AG += 1
                total += 1
        if origin == "T":
            if alt == "A":
                TA += 1
                total += 1
            if alt == "C":
                TC += 1
                total += 1
            if alt == "G":
                TG += 1
                total += 1
        if origin == "G":
            if alt == "T":
                GT += 1
                total += 1
            if alt == "C":
                GC += 1
                total += 1
            if alt == "A":
                GA += 1
                total += 1
        if origin == "C":
            if alt == "T":
                CT += 1
                total += 1
            if alt == "A":
                CA += 1
                total += 1
            if alt == "G":
                CG += 1
                total += 1
        
output_file.write(
"total number of snp's"+"\t"+str(total)+"\n"
"long region subsitutions"+"\t"+str(long_region)+"\n"
"double mutations "+"\t"+str(double)+"\n"
"A > T"+"\t"+str(AT)+"\n"
"A > G"+"\t"+str(AG)+"\n"
"A > C"+"\t"+str(AC)+"\n"
"T > A"+"\t"+str(TA)+"\n"
"T > G"+"\t"+str(TG)+"\n"
"T > C"+"\t"+str(TC)+"\n"
"G > A"+"\t"+str(GA)+"\n"
"G > T"+"\t"+str(GT)+"\n"
"G > C"+"\t"+str(GC)+"\n"
"C > A"+"\t"+str(CA)+"\n"
"C > T"+"\t"+str(CT)+"\n"
"C > G"+"\t"+str(CG)+"\n"
)

input_file.close()
output_file.close()

