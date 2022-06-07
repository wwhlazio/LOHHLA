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


def check_bam(samplelist,datadir,nbam_file):

    f=open(samplelist,'r')

    samp=[]
    ns=0
    for line in f:
        if ns >0:
            tmp=line.split(",");
            samp.append(tmp[0])
        ns=ns+1

    f.close()
    ns=ns-1

    print("number of sampls%d\n" % ns)
    fx=open(nbam_file,"w")
    
    for i in range(ns) :
        samp[i]=samp[i].strip("\"")
        cmd="find %s -name \"*%s*.bam\"" % (datadir,samp[i])
	print(cmd)
        a=subprocess.check_output(cmd,shell=True).decode("utf-8")
        b=a.split("\n")
        c='|'.join(b)
        d=list()
        for j in range(len(b)-1):
            d.append(check_disease(b[j]))
        dd='|'.join(d)
        fx.write("%s\t%d\t%s\t%s\n" % (samp[i],len(b)-1,c,dd))
    fx.close()


if __name__=="__main__":
    
    datadir="/sc/arion/projects/zhuj05a/Wenhui/bladder/tcga_blca/"
    samplelist="/sc/arion/projects/zhuj05a/Wenhui/bladder/tcga_blca_clinical/GDCdata/tcga_blca_clincial_x.csv"
    nbam_file="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/tcga_blca_bam"
    check_bam(samplelist,datadir,nbam_file) 
