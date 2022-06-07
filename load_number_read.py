import os
import subprocess
import re
import sys


def check_disease(filename):

    #tind=['01','02','03','04','05','06','07','08','09']
    #nind=['10','11','12','13','14','15','16','17','18','19']
    tmp=filename.split('/')
    tmpx=tmp[len(tmp)-1]
    #print(tmpx)
    tmpxx=tmpx.split(".")
    tmpxxx=tmpxx[1]
    tmpxxxx=tmpx.split("-")
    tmpx5=tmpxxxx[3]
    idx=re.findall('\d+', tmpx5)
    idx=''.join(idx)
    #print(tmpx5)
    #print(idx)
    #print(type(idx))
    if idx == '01':
        return "tumor"
    elif idx == '10':
        return "blood"
    elif idx =='11':
        return 'normal'
    else:
        return 'NA'


def load_read_number(workdir,outdir,samp,out):

    readfile=outdir+samp+"/coutread_"+samp
    lsfile=workdir+samp+"/countread_"+samp+".lsf"

    fx=open(out,'a')

    f=open(lsfile,'r')

    bam=[]
    for line in f:
        line=line.strip()
        if "samtools view -c -F 4" in line:
            tmp=line.split()
            bam.append(tmp[5])
    f.close()
    bam=list(set(bam))   
    
    f=open(readfile,'r')
    c=0
    read=[]
    nb=0
    for line in f:
        line=line.splitlines()
        read.append(''.join(line))
        c+=1
        if c==4:
            c=0
            fx.write("%s\t" % samp)
            fx.write("%s\t" % bam[nb])
            fx.write("%s\t" % check_disease(bam[nb]))
            for i in range(3):
                fx.write("%s\t" % read[i])
            fx.write("%s\n" % read[3])
            read=[]
            nb=nb+1

    fx.close()

def load_all_sample(samplelist,workdir,outdir,out):
    
    f=open(samplelist,'r')

    samp=[]                                     #load all 412 samples
    ns=0
    for line in f:
        if ns >0:
            tmp=line.split(",");
            samp.append(tmp[0])
        ns=ns+1

    f.close()
    ns=ns-1

    for i in range(ns):
        samp[i]=samp[i].strip('\"')
        print(i)
        load_read_number(workdir,outdir,samp[i],out) 
        
if __name__=="__main__":
    workdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/workdir/countread/"    
    outdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/output/countread/"
    samplelist="/sc/arion/projects/zhuj05a/Wenhui/bladder/tcga_blca_clinical/GDCdata/tcga_blca_clincial_x.csv"
    outfile="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/bam_hla_readcount"
    f=open(outfile,'w')
    f.write("sample\tfile\ttissue\tall_read\thla_a\thla_b\thla_c\n")
    f.close()
    load_all_sample(samplelist,workdir,outdir,outfile)
