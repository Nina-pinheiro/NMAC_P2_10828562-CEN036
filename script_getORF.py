# To add a new cell, type '# %%'
# To add a new markdown cell, type '# %% [markdown]'
# %% [markdown]
# ### Testando importanção de arquivos do módulo sys

# %%
import sys
import pandas as pd

caminho = "sequencia.fasta"

arquivo = open(caminho, 'r')

sequencia = arquivo.readlines()

print(sequencia)

# %% [markdown]
# ### Fazendo o exercício

# %%
#Lendo o arquivo e gerando a sequência complementar pelo

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

with open("sequencia.fasta", "r") as seq_file:
        for record in SeqIO.parse(seq_file, "fasta"):
            seq_dna = record.seq
            complementar = seq_dna.reverse_complement()

print(seq_dna)
print(complementar)


# %%
#Função para gerar os ORF's
def getORFs(dna, frame):
    #frame_count = 0
    for i in range(frame, len(dna), 3):
        codon1 = dna[i:i+3]
        if codon1 == 'atg':
            start = i
            for j in range(start, len(dna), 3):
                codon2 = dna[j:j+3]
                if codon2 in ['taa', 'tag', 'tga']:
                    end = j
                    orflength = end-start+3
                    # frame += 1
                    # print(frame) 
                    #100 - 0+3 = 97
                    orf = dna[start:end+3]
                    yield [orflength, orf, start, end+3] #frame_count]

#1 - Criando a lista com o ORF's
sequencias = [seq_dna, complementar]
orfs = []
frames = [0,1,2]

for a in sequencias:
    for i in frames:
        print(i)
        for j in getORFs(a, i):
            print(j)
            orfs.append(j)
print(orfs)


# %%
#2 - Obtendo o ORF de maior comprimento

maxLengthORF = [] 
from operator import itemgetter

def getMaxLengthORF(orfList, i):
    return max(enumerate(map(itemgetter(i), orfList)), key=itemgetter(1))

maxLengthORF.append(getMaxLengthORF(orfs, 0))
print(maxLengthORF)


# %%
#3 - Traduzindo a sequência do pepitídio do ORF's com o maior comprimento

maxLengthORFIndex = maxLengthORF[0][0]
novo_maxLengthORFIndex = maxLengthORF[0][1]
orfToTranslate = orfs[1][maxLengthORFIndex]
print(orfToTranslate.translate(stop_symbol=""))


# %%
#4,5 e 6 - Salvar arquivo 

#Gerando o nome do arquivo para salvar no seu id
s = 'id_sequencia01_%s_%i_%i' % (1, orfs[1][2], orfs[1][3])

#Salvando o arquivo ORF.fna
ORF_fna = pd.DataFrame({'%s' %(s): [orfs[1][1]]})
ORF_fna.to_csv('ORF.fna')

#Salvando o arquivo ORF.faa
ORF_faa = pd.DataFrame({'%s' %(s): [orfs[1][1].translate(stop_symbol="")]})
ORF_faa.to_csv('ORF.faa')


# %%
#7 - Carregar arquivos no repositório


# %%
#8 - verificar integridade dos arquivos

import os

def index(directory):
    stack = [directory]
    files = []
    while stack:
        directory = stack.pop()
        for file in os.listdir(directory):
            fullname = os.path.join(directory, file)
            if fullname.endswith('mp3'):
                files.append(fullname)
            if os.path.isdir(fullname) and not os.path.islink(fullname):
                stack.append(fullname)
    return files

def check(directory):
    files = index(directory)
    hvalues = []
    for x in files:
        cmd = 'md5sum' + ' ' + x
        fp = os.popen(cmd)
        res = fp.readline()
        hvalues.append(res)
        stat = fp.close() 
    return hvalues

def check(directory):
    files = index(directory)
    hvalues = []
    for f in files:
        cmd = ['md5sum', f]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        hvalues.append(proc.stdout.readline())
        proc.stdout.close()
        stat = os.waitpid(proc.pid, 0)
    return hvalues


# %%
index('ORF.fna')
index('ORF.faa')
check('ORF.fna')
check('ORF.faa')


# %%



