import argparse
import random

def arguments():
        parser = argparse.ArgumentParser(description="misum on ms")
        parser.add_argument("-i", "--input", help="path to ms file", required=True)
        parser.add_argument('-o', '--output', help='output path with file name omitted. file name will be the input file name.', required=True)
        parser.add_argument('-n', "--haploids", help='total number of haploids simulated, including simulated reference', required=True, type=int)
        args = parser.parse_args()
        return(args)

def divergence(x,y):

        return(sum ( x[i] != y[i] for i in range(len(x)) ))

def simsum(lines, id, outpath, haploids):
        outfile = outpath + id + '.simsum'
        outfile = open(outfile, 'w')
        outfile.write('minsum' + '\t' + 'maxsum' + '\t' + 'segsites' + '\t' + 'SNPs' + '\n')
        pairs=[]
        i=haploids+20
        for line in lines:
                if '//' in line:
                        i=haploids+3
                # Put all of the haplotypes into a list and calc minsum on pairs
                # Pairs function as diploids

                if i==1:
                        r=line
                        x=True
                        n=0
                        while x:
                                g1=pairs[n]
                                g2=pairs[n+1]
                                hap1 = divergence(g1, r)
                                hap2 = divergence(g2, r)
                                if hap1>hap2:
                                        max = hap1
                                        min = hap2
                                elif hap2>hap1:
                                        max = hap2
                                        min = hap1
                                elif hap2==hap1:
                                        choice=random.randint(1,2)
                                        if choice==1:
                                                max = hap2
                                                min = hap1
                                        elif choice==2:
                                                max = hap1
                                                min = hap2
                                snps = divergence(g1,g2)
                                outfile.write(str(min) + '\t' + str(max) + '\t' + str(len(g1)-1) + '\t' + str(snps) + '\n')
                                if g2==pairs[-1]:
                                        x=False
                                        pairs=[]
                                g1=[]
                                g2=[]
                                hap1=[]
                                hap2=[]
                                n+=2
                                i=haploids+20
                        r=[]
                        pairs=[]

                elif i<=haploids and i>1:
                        pairs.append(line)
                i-=1
        outfile.close()

##############################################
args = arguments()

# get infile info
input = args.input
infile = open(input, 'r')
lines = infile.readlines()

# get outpath info
outpath = args.output
haploids=args.haploids

#parse out name of indiv for outfile name
id = input.split('/')[-1]

print 'running simsum'
simsum(lines, id, outpath, haploids)
print 'simsum finished'

infile.close()
