import argparse
# import os
​
parser = argparse.ArgumentParser("C:/users/no/experimentsbyid/cp121/CP121Draft1.py")
parser.add_argument('-inp', dest="BEDInputFile", help='Path to your .BED input file.')
parser.add_argument('-dir', dest="FileDir", help='Path to desired output folder. E.g., "C:/Users/Documents/" ')
parser.add_argument('-len', dest="PromoterLength", default=3000, type=int, help='Length of promoter sequence. Defaults to 3000bp.')
parser.add_argument('-tm', dest="TmTarget", default=62, type=int, help='Ideal primer melting temperature. Defaults to 62C.')
parser.add_argument('-tmlow', dest="TmLow", default=60, type=int, help="Lowest acceptable primer temperature. Defaults to 60C.")
parser.add_argument('-tmhi', dest="TmHigh", default=70, type=int, help="Highest acceptable primer temperature. Defaults to 70C.")
parser.add_argument('-lenmin', dest="MinPrimerLength", default=18, type=int, help="Minimum primer length. Defaults to 18nt.")
parser.add_argument('-buf', dest="BufferID", default="q5-0", help="Sets buffer/polymerase ID, defaults to Q5. Options are not defined yet. See manual at https://tmapi.neb.com/#usage for more info.")
#Add prodcode arguments.
parser.add_argument('-conc', dest="PrimerConc", default=.5, type=int, help="Sets primer concentration in nM. Default is .5 uM.")
args = vars(parser.parse_args())
​
#Primermaker:
#The input file should be in .BED format (i.e., with columns: ["","chrom","chromStart","chromEnd","gene_name","id_number","strand"]).
​
from Bio import Entrez, SeqIO
import csv
​
# This makes the conversion dictionary for chromosomes to accession number:
​
Keys = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","X","Y"]
Values = ['NC_000001.11', 'NC_000002.12', 'NC_000003.12', 'NC_000004.12', 'NC_000005.10', 'NC_000006.12',
          'NC_000007.14', 'NC_000008.11', 'NC_000009.12', 'NC_000010.11', 'NC_000011.10', 'NC_000012.12',
          'NC_000013.11', 'NC_000014.9', 'NC_000015.10', 'NC_000016.10', 'NC_000017.11', 'NC_000018.10', 'NC_000019.10',
          'NC_000020.11', 'NC_000021.9', 'NC_000022.11', 'NC_000023.11', 'NC_000024.10']
​
ReferenceDict = dict(zip(Keys, Values))
​
#Entering email for NCBI reference:
​
Entrez.email = "penneyc@email.chop.edu"
​
#This function takes an Inputfile (.BED format) path as well as a desired sequence length as an integer, and outputs "gene_name"+","+"SEQUENCE"+\n to a given output file:
​
def SequenceGrabber(InputFile, Length, FileDir):
    with open(InputFile, 'r') as inputfile, open(FileDir + "Output.txt", 'w') as outputfile:
        readCSV = csv.reader(inputfile, delimiter=',')
        next(inputfile)
        for row in readCSV:
            handle = Entrez.efetch(db="nucleotide",
                                   id=ReferenceDict[row[1]],
                                   rettype="fasta",
                                   strand=row[6],
                                   seq_start=row[2],
                                   seq_stop=(int(row[2]) - Length))
            record = SeqIO.read(handle, "fasta")
            handle.close()
            outputfile.write(str(row[4]) + "," + str(record.seq) + "\n")
​
import requests
import json
import csv
​
def GatherCharacteristics(Sequences):
    url = 'https://tmapi.neb.com/tm/batch'
    seqpairs = Sequences
    input = {"seqpairs": seqpairs, 'conc': args['PrimerConc'], 'prodcode': args['BufferID']}
    headers = {'content-type': 'application/json'}
    res = requests.post(url, data=json.dumps(input), headers=headers)
    r = json.loads(res.content)
    if r['success']:
        for row in r['data']:
            #print('Seq1: {}  Tm1: {}  Seq2: {} Tm2: {}  Ta: {}'.format(row['seq1'], row['tm1'], row['seq2'], row['tm2'], row['ta']))
            return (row['tm1'], row['tm2'], row['ta'])
    else:
        print('request failed')
        print(r['error'][0])
​
def PrimerGenerator(InputFile, OutputFile, TmTarget, TmLow, TmHigh, MinLength):
    with open(InputFile, 'r') as csvfile, open(OutputFile, 'w') as output:
        readCSV = csv.reader(csvfile, delimiter=',')
        for row in readCSV:
            Gene_Name = row[0]
            Promoter_Sequence = row[1]
            Sequences = [[Promoter_Sequence[0:35], Promoter_Sequence[-35:-1]]]
            print(Sequences)
            ErrorAnneal = ""
            ErrorTm1 = ""
            ErrorTm2 = ""
            A = False
            while A == False:
                if (GatherCharacteristics(Sequences)[0] != TmTarget):
                    if len(Sequences[0][0]) > MinLength:
                        Sequences[0][0] = Sequences[0][0][0:(len(Sequences[0][0]) - int(1))]
                        print(len(Sequences[0][0]), Gene_Name)
                    elif len(Sequences[0][0]) == MinLength:
                        if TmLow < (GatherCharacteristics(Sequences)[0]) < TmHigh:
                            A = True
                            Tm1=GatherCharacteristics(Sequences)[0]
                            Length1=len(Sequences[0][0])
                            print(Tm1, Length1, Gene_Name)
                        else:
                            A = True
                            ErrorTm1 = "Primer 1 Tm Outside Bounds."
                            Tm1 = GatherCharacteristics(Sequences)[0]
                            Length1 = len(Sequences[0][0])
                            print(ErrorTm1, Gene_Name)
                            print(Tm1)
                elif (GatherCharacteristics(Sequences)[0] == int(TmTarget)):
                    if len(Sequences[0][0]) > int(MinLength):
                        A = True
                        Tm1 = GatherCharacteristics(Sequences)[0]
                        Length1 = len(Sequences[0][0])
                        print(Tm1, Length1, Gene_Name)
                    else:
                        A = True
                        Tm1 = GatherCharacteristics(Sequences)[0]
                        Length1 = len(Sequences[0][0])
                        print("Error Making Primer1" + Gene_Name)
            B = False
            while B == False:
                if (GatherCharacteristics(Sequences)[1] != TmTarget):
                    if len(Sequences[0][1]) > MinLength:
                        Sequences[0][1] = Sequences[0][1][0:(len(Sequences[0][1]) - int(1))]
                        print(len(Sequences[0][1]), Gene_Name)
                        Tm2 = GatherCharacteristics(Sequences)[1]
                        Length2 = len(Sequences[0][1])
                        Ta = GatherCharacteristics(Sequences)[2]
                    elif len(Sequences[0][1]) == MinLength:
                        if TmLow <= (GatherCharacteristics(Sequences)[1]) <= TmHigh:
                            if (Tm1 - 5) <= GatherCharacteristics(Sequences)[1] <= (Tm1 + 5):
                                Tm2 = GatherCharacteristics(Sequences)[1]
                                Length2 = len(Sequences[0][1])
                                Ta = GatherCharacteristics(Sequences)[2]
                                print(Tm2, Length2, Ta, Gene_Name)
                                B = True
                            else:
                                B = True
                                ErrorTm2 = "Primer 2 " + Gene_Name + " Tm Outside Bounds."
                                Tm2 = GatherCharacteristics(Sequences)[1]
                                Length2 = len(Sequences[0][1])
                                Ta = GatherCharacteristics(Sequences)[2]
                                print("Reverse Primer not within 5C of first." + Gene_Name , GatherCharacteristics(Sequences)[1])
                        else:
                            B=True
                            print("Error: primer2 out of Tm range")
                elif (GatherCharacteristics(Sequences)[1] == TmTarget):
                    if len(Sequences[0][1]) > MinLength:
                        B = True
                        Tm2 = GatherCharacteristics(Sequences)[1]
                        Length2 = len(Sequences[0][1])
                        Ta = GatherCharacteristics(Sequences)[2]
                        print(Tm2, Length2, Ta)
                    else:
                        B = True
                        print("Error Making Primer2")
            if Ta > Tm1 or Ta > Tm2:
                AnnealError = True
                Error = "Annealing Temperature Higher than Melt Temp"
            output.write(Gene_Name + "," + Sequences[0][0] + "," + str(Tm1) + "," + Sequences[0][1] + "," + str(Tm2) + "," + str(Ta) + "," + ErrorAnneal + "," + ErrorTm1 + "," + ErrorTm2 + "\n")
​
SequenceGrabber(InputFile = args['BEDInputFile'], Length=args['PromoterLength'], FileDir = args['FileDir'])
PrimerGenerator(InputFile=(args['FileDir'] + "Output.txt"), OutputFile=(args['FileDir'] + "PrimersOutput.txt"), TmTarget=args['TmTarget'], TmLow = args['TmLow'], TmHigh = args['TmHigh'], MinLength= args['MinPrimerLength'])
