#!/usr/bin/python
#play with samples having tumor and solid normal

#for tumor and normal, update from verion v1
import os
import subprocess
import re
import sys


def merge_count(count_dir,dep_dir,samplelist,workdir,errorinputfile):
    
    acc = "premium"
    mem1 = "20000"
    time1 = "01:00"

    f=open(samplelist,'r')

    samp=[]					#load all the 412 samples
    ns=0
    for line in f:
        if ns >0:
            tmp=line.split(",");
            samp.append(tmp[0])
        ns=ns+1

    f.close()
    ns=ns-1
    
    print("number of sampls%d\n" % ns)
    dis=["tumor","normal"]
    
    sh_file=workdir+"run_merge_count.sh"
    fx=open(sh_file,'w')
    fxx=open(errorinputfile,'w')
    for i in range(ns):
        samp[i]=samp[i].strip('\"')
        cmd="find %s -name \"*%s*.bam\"" % (datadir,samp[i])		#find samples with bam file
        a=subprocess.check_output(cmd,shell=True).decode("utf-8")
        b=a.split("\n")
        if len(b) > 4:							#filter samples with 3 or more bams
            print("Error number of bam files!")
            print(a)
            continue


        dis=list()
        for j in range(len(b)-1):
            tmp=check_disease2(b[j])
            if tmp =='tumor' or tmp == 'normal':			#find samples with one tumor and one normal
                dis.append(check_disease2(b[j]))
        if len(dis)<2:
            continue
        
        for j in range(2):
            samp_countdir=count_dir+samp[i]+"_"+dis[j]+"/"
            #print(samp_countdir)
            if not os.path.exists(samp_countdir):
                fxx.write("%s\n" % samp[i]) 
            else:
                lsf=workdir+samp[i]+"_"+dis[j]+"_mergecount.lsf"
                outputfile=dep_dir+samp[i]+"_"+dis[j]
                #inputdir=count_dir+samp[i]+"_"+dis[j]+"/"
                f=open(lsf,'w')
                f.write("#!/bin/bash\n")
                f.write("#BSUB -P acc_BD2K\n")
                f.write("#BSUB -q %s\n" % acc)
                f.write("#BSUB -n 1\n")
                f.write("#BSUB -J mergec_%s\n" % i)
                f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
                f.write("#BSUB -W %s\n" % time1)
                f.write("#BSUB -o %s/mergec_%s_%s.stdout\n" % (workdir,samp[i],dis[j]))
                f.write("#BSUB -eo %s/mergec_%s_%s.stderr\n" % (workdir,samp[i],dis[j]))
                f.write("#BSUB -L /bin/bash\n")
                f.write("module load R\n")
                f.write("Rscript -e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/code/merge_pos.R\\\")\" -e \"merge_pos(\\\"%s\\\",\\\"%s\\\")\"\n" % (samp_countdir,outputfile))
                f.close()
                fx.write("bsub < %s\n" % lsf)
    fx.close()
    fxx.close()
		
def run_gen_input(dep_dir,act_dir,samplelist,workdir,cov_cut):
    acc = "premium"
    mem1 = "20000"
    time1 = "01:00"

    f=open(samplelist,'r')

    samp=[]				#load all samples
    ns=0
    for line in f:
        if ns >0:
            tmp=line.split(",");
            samp.append(tmp[0])
        ns=ns+1

    f.close()
    ns=ns-1

    print("number of sampls%d\n" % ns)
    sh_file=workdir+"run_gen_act.sh"
    fx=open(sh_file,'w')
    for i in range(ns):
        samp[i]=samp[i].strip('\"')
        norm_input=dep_dir+samp[i]+"_normal"
        tumor_input=dep_dir+samp[i]+"_tumor"
        print(norm_input)
        print(tumor_input)
        if os.path.isfile(norm_input) and os.path.isfile(tumor_input):
            lsf=workdir+samp[i]+"_genact.lsf"
            outputfile=act_dir+samp[i]
            f=open(lsf,'w')
            f.write("#!/bin/bash\n")
            f.write("#BSUB -P acc_BD2K\n")
            f.write("#BSUB -q %s\n" % acc)
            f.write("#BSUB -n 1\n")
            f.write("#BSUB -J genact_%s\n" % i)
            f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
            f.write("#BSUB -W %s\n" % time1)
            f.write("#BSUB -o %s/geneact_%s.stdout\n" % (workdir,samp[i]))
            f.write("#BSUB -eo %s/geneact_%s.stderr\n" % (workdir,samp[i]))
            f.write("#BSUB -L /bin/bash\n")
            f.write("module load R\n")
            f.write("Rscript -e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/code/generate_ascat_input.R\\\")\" -e \"generate_ascat_input(\\\"%s\\\",\\\"%s\\\",\\\"%s\\\",\\\"%s\\\")\"\n" % (tumor_input,norm_input,outputfile,cov_cut))
            f.close()
            fx.write("bsub < %s\n" % lsf)
    fx.close()



    
    

def check_disease(filename):

    tind=['01','02','03','04','05','06','07','08','09']
    nind=['10','11','12','13','14','15','16','17','18','19']
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
    if idx in tind:
        return "tumor"
    elif idx in nind:
        return "normal"
    else:
        print("Error!")
        return
        

def check_disease2(filename):

    tind='01'
    nind='11'
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
    if idx == tind:
        return "tumor"
    elif idx == nind:
        return "normal"
    else:
        return "NA"

     

def run_alleleCount(datadir,workdir,outdir,bindir,samplelist,poslist):
    
    acc = "premium"
    mem1 = "10000"
    time1 = "00:20"


    f=open(samplelist,'r')
    
    samp=[]					#load all 412 samples
    ns=0
    for line in f:
        if ns >0:
            tmp=line.split(",");
            samp.append(tmp[0])
        ns=ns+1

    f.close()
    ns=ns-1
    
    f=open(poslist,'r')
     
    pos=[]					#load position files
    np=0
    possuf=[]
    for line in f:
        tmp=line.split('\n')
        pos.append(tmp[0])
        tmpx=tmp[0].split('.')
        possuf.append(tmpx[len(tmpx)-1])
        np=np+1

    


    print("number of sampls%d\n" % ns)
    print("number of poslist%d\n" % np)

    
    sh_file = workdir + "run_allelecount.sh"
    fxx=open(sh_file,'w')
    for i in range(ns):
       
        samp[i]=samp[i].strip('\"')
        cmd="find %s -name \"*%s*.bam\"" % (datadir,samp[i])
        a=subprocess.check_output(cmd,shell=True).decode("utf-8")
        b=a.split("\n")
        if len(b) > 4:							#filter samples with 3 or more bam files
            print("Error number of bam files!") 
            print(a)
            continue
        
        
        dis=list()							#only keep samples with one normal and one tumor
        disfile=list()
        for j in range(len(b)-1):
            tmp=check_disease2(b[j])
            if tmp =='tumor' or tmp == 'normal':
                dis.append(check_disease2(b[j]))
                disfile.append(b[j])
        if len(dis)<2:
            continue
        
        for j in range(2):
            samp_workdir=workdir+samp[i]+"_"+dis[j]
            if not os.path.exists(samp_workdir):
                os.makedirs(samp_workdir)
            samp_outdir=outdir+samp[i]+"_"+dis[j]
            if not os.path.exists(samp_outdir):
                os.makedirs(samp_outdir)
            samp_sh=samp_workdir+"/run_allelecount_"+samp[i]+"_"+dis[j]+".sh" 
            #print(samp_workdir)
            if not os.path.isdir(samp_workdir): 
                print("Error of samp_workdir\n")
                sys.exit()
            #print(samp_outdir)
            #if not os.path.exists(samp_outdir):
            #    print("Error of samp_outdir\n")
            fx=open(samp_sh,"w")
            for k in range(np):
                samp_lsf=samp_workdir+"/allelec_"+samp[i]+"_"+dis[j]+"_"+possuf[k]+".lsf"
                f=open(samp_lsf,'w')
                f.write("#!/bin/bash\n")
                f.write("#BSUB -P acc_BD2K\n")
                f.write("#BSUB -q %s\n" % acc)
                f.write("#BSUB -n 1\n")
                f.write("#BSUB -J allelec_%s\n" % i)
                f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
                f.write("#BSUB -W %s\n" % time1)
                f.write("#BSUB -o %s/allelec_%s_%s_%s.stdout\n" % (samp_workdir,samp[i],dis[j],possuf[k]))
                f.write("#BSUB -eo %s/allelec_%s_%s_%s.stderr\n" % (samp_workdir,samp[i],dis[j],possuf[k]))
                f.write("#BSUB -L /bin/bash\n")
                f.write("export PERL5LIB=$PERL5LIB:/hpc/packages/minerva-common/bioperl/1.6.619/lib/perl5\n")
                f.write("export PERLLIB=$PERLLIB:/hpc/packages/minerva-common/bioperl/1.6.619/lib/perl5\n")
                f.write("export PERL5LIB=$PERL5LIB:/hpc/users/wangm08/perl5/lib/perl5\n")
                f.write("export PERLLIB=$PERLLIB:/hpc/users/wangm08/perl5/lib/perl5\n")
                f.write("export PATH=$PATH:/sc/arion/projects/zhuj05a/Wenhui/packages/alleleCount/bin\n")
                f.write("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/sc/arion/projects/zhuj05a/Wenhui/packages/alleleCount/lib\n")
                f.write("export PERL5LIB=$PERL5LIB:/sc/arion/projects/zhuj05a/Wenhui/packages/alleleCount/lib/perl5/\n")
                f.write("export PERLLIB=$PERLLIB:/sc/arion/projects/zhuj05a/Wenhui/packages/alleleCount/lib/perl5/\n")
                f.write("module load htslib\n")
                f.write("cd %s\n" % bindir)
                f.write("./alleleCounter -l %s -b %s -o %s/%s_%s_%s\n" % (pos[k],disfile[j],samp_outdir,samp[i],dis[j],possuf[k]))
                f.close()
                fx.write("bsub < %s\n" % samp_lsf)
            fx.close()
            fxx.write("chmod +x %s\n" % samp_sh)
            fxx.write("%s\n" % samp_sh)
    
    fxx.close()

def run_polysolver(datadir,samplist,workdir,outdir):

    acc = "premium"
    mem1 = "10000"
    time1 = "02:00"


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
    
    print(ns)

    sh_file = workdir + "run_polysolver.sh"
    fxx=open(sh_file,'w')

    for i in range(ns):

        samp[i]=samp[i].strip('\"')
        cmd="find %s -name \"*%s*.bam\"" % (datadir,samp[i])
        a=subprocess.check_output(cmd,shell=True).decode("utf-8")
        b=a.split("\n")
        if len(b) > 4:
            print("Error number of bam files!")
            print(a)
            continue

        dis=list()
        disfile=list()
        for j in range(len(b)-1):
            tmp=check_disease2(b[j])
            if tmp =='tumor' or tmp == 'normal':
                dis.append(check_disease2(b[j]))
                disfile.append(b[j])
        #print(dis)
        if len(dis)<2:
            continue
        

        for j in range(len(dis)):
            if dis[j] == "normal":
                ind=j
        #print(ind)
        samp_workdir=workdir+samp[i]        
        samp_lsf=samp_workdir+"/polysolver_"+samp[i]+".lsf"
        samp_outdir=outdir+samp[i]


        if not os.path.exists(samp_workdir):
            os.makedirs(samp_workdir)
        if not os.path.exists(samp_outdir):
            os.makedirs(samp_outdir)

        
        f=open(samp_lsf,'w') 
        f.write("#!/bin/bash\n")
        f.write("#BSUB -P acc_BD2K\n")
        f.write("#BSUB -q %s\n" % acc)
        f.write("#BSUB -n 1\n")
        f.write("#BSUB -J polysolver_%s\n" % i)
        f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
        f.write("#BSUB -W %s\n" % time1)
        f.write("#BSUB -o %s/polysolver_%s.stdout\n" % (samp_workdir,samp[i]))
        f.write("#BSUB -eo %s/polysolver_%s.stderr\n" % (samp_workdir,samp[i]))
        f.write("#BSUB -L /bin/bash\n")
        f.write("export PERL5LIB=$PERL5LIB:/hpc/packages/minerva-common/bioperl/1.6.619/lib/perl5\n")
        f.write("export PERLLIB=$PERLLIB:/hpc/packages/minerva-common/bioperl/1.6.619/lib/perl5\n")
        f.write("export PERL5LIB=$PERL5LIB:/hpc/users/wangm08/perl5/lib/perl5\n")
        f.write("export PERLLIB=$PERLLIB:/hpc/users/wangm08/perl5/lib/perl5\n")
        f.write("source /sc/arion/projects/zhuj05a/Wenhui/packages/polysolver/polysolver/scripts/config.bash\n")
        f.write("cd %s\n" % samp_workdir)
        f.write("/sc/arion/projects/zhuj05a/Wenhui/packages/polysolver/polysolver/scripts/shell_call_hla_type %s Unknown 1 hg38 STDFQ 0 %s >& log1\n" % (disfile[ind],samp_outdir))
        f.close()
        fxx.write("bsub < %s\n" % samp_lsf)
    fxx.close()

def run_ascat(samplelist,datadir,workdir,outdir):
    
    acc = "premium"
    mem1 = "10000"
    time1 = "02:00"


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

    sh_file = workdir + "run_ascat.sh"
    fxx=open(sh_file,'w')

    for i in range(ns):
         samp[i]=samp[i].strip('\"')
         asc_input=datadir+samp[i]
         
         if os.path.isfile(asc_input):
            workdir_samp=workdir+samp[i]+"/"
            if not os.path.exists(workdir_samp):
                os.makedirs(workdir_samp)   
            lsf=workdir_samp+samp[i]+"_ascat.lsf"
            outputfile=outdir+samp[i]
            f=open(lsf,'w')
            f.write("#!/bin/bash\n")
            f.write("#BSUB -P acc_BD2K\n")
            f.write("#BSUB -q %s\n" % acc)
            f.write("#BSUB -n 1\n")
            f.write("#BSUB -J ascat_%d\n" % i)
            f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
            f.write("#BSUB -W %s\n" % time1)
            f.write("#BSUB -o %s/ascat_%s.stdout\n" % (workdir_samp,samp[i]))
            f.write("#BSUB -eo %s/ascat_%s.stderr\n" % (workdir_samp,samp[i]))
            f.write("#BSUB -L /bin/bash\n")
            f.write("module load R\n")
            f.write("Rscript -e \"source(\\\"/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/code/run_ascat.R\\\")\" -e \"a=read.table(\\\"%s\\\",header=T,as.is=T)\" -e \"run_ascat(a,\\\"%s\\\",\\\"%s\\\",\\\"%s\\\")\"\n" % (asc_input,workdir_samp,outputfile,samp[i]))
            f.close()
            fxx.write("bsub < %s\n" % lsf)
    fxx.close()

def run_lohhla(samplelist,data_bam,data_cnv,data_allele,workdir,outdir):
    
    acc = "premium"
    mem1 = "10000"
    time1 = "02:00"


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

    sh_file = workdir + "run_lohhla.sh"
    fxx=open(sh_file,'w')
    for i in range(ns):
        samp[i]=samp[i].strip('\"')
        samp_allele=data_allele+samp[i]+"/winners.hla.txt"        
        samp_cnv=data_cnv+samp[i]
        
        if os.path.isfile(samp_allele) and  os.path.isfile(samp_cnv):
            samp_work=workdir+samp[i]+"/"
            samp_out=outdir+samp[i]+"/"
            if not os.path.exists(samp_work):
                os.makedirs(samp_work)
            if not os.path.exists(samp_out):
                os.makedirs(samp_out)
            lsf=samp_work+samp[i]+"_lohhla.lsf"
            f=open(lsf,'w')
            f.write("#!/bin/bash\n")
            f.write("#BSUB -P acc_BD2K\n")
            f.write("#BSUB -q %s\n" % acc)
            f.write("#BSUB -n 1\n")
            f.write("#BSUB -J lohhla_%d\n" % i)
            f.write("#BSUB -R \"rusage[mem=%s]\"\n" % mem1)
            f.write("#BSUB -W %s\n" % time1)
            f.write("#BSUB -o %s/lohhla_%s.stdout\n" % (samp_work,samp[i]))
            f.write("#BSUB -eo %s/lohhla_%s.stderr\n" % (samp_work,samp[i]))
            f.write("#BSUB -L /bin/bash\n")
            f.write("module load BEDTools/2.26.0\n")
            f.write("module load samtools/1.3\n")
            f.write("module load R/3.6.0\n")
            f.write("module load jellyfish/2.2.6\n")
            f.write("export PATH=$PATH:\"/sc/arion/projects/zhuj05a/Wenhui/packages/novoalign/novocraft/\"\n")
            
            cmd="find %s -name \"*%s*.bam\"" % (data_bam,samp[i])           
            #print(cmd)
            a=subprocess.check_output(cmd,shell=True).decode("utf-8")
            b=a.split("\n")
            if len(b) > 4:
                print("Error number of bam files!")
                print(a)
                continue
            
            print(b) 
            dis=list()
            disfile=list()
            for j in range(len(b)-1):
                tmp=check_disease2(b[j])
                if tmp =='tumor' or tmp == 'normal':
                    dis.append(check_disease2(b[j]))
                    disfile.append(b[j])
                if len(dis)<2:
                    continue    
            ind1=0
            ind2=0
            print(dis)
            for j in range(len(dis)):
                print(dis[j])
                if dis[j] == "normal":
                    ind1=j
                elif dis[j] == "tumor":
                    ind2=j

            f.write("cd %s\n" % samp_work)
            f.write("ln -s %s %s%s_normal.bam\n" % (disfile[ind1],samp_work,samp[i]))
            f.write("ln -s %s %s%s_tumor.bam\n" % (disfile[ind2],samp_work,samp[i]))
            
            tmp1=disfile[ind1]
            tmp1=tmp1[:-1]
            tmp1=tmp1+"i"
            tmp2=disfile[ind2]
            tmp2=tmp2[:-1]
            tmp2=tmp2+"i"
            
            f.write("ln -s %s %s%s_normal.bai\n" % (tmp1,samp_work,samp[i]))
            f.write("ln -s %s %s%s_tumor.bai\n" % (tmp2,samp_work,samp[i]))
            	    	            
            f.write("Rscript /sc/arion/projects/zhuj05a/Wenhui/packages/LOHHLA/lohhla/LOHHLAscript.R \\\n")
            f.write("--patientId %s \\\n" % samp[i])
            f.write("--outputDir %s \\\n" % samp_out) 
            f.write("--normalBAMfile %s%s_normal.bam \\\n" % (samp_work,samp[i]))
            f.write("--BAMDir %s \\\n" % samp_work)
            f.write("--hlaPath %s \\\n" % samp_allele)
            f.write("--HLAfastaLoc /sc/arion/projects/zhuj05a/Wenhui/packages/polysolver/polysolver/data/abc_complete.fasta \\\n")
            f.write("--CopyNumLoc %s \\\n" % samp_cnv)
            f.write("--mappingStep TRUE \\\n")
            f.write("--minCoverageFilter 10 \\\n")
            f.write("--fishingStep TRUE \\\n")
            f.write("--cleanUp FALSE \\\n")
            f.write("--gatkDir /hpc/packages/minerva-common/picard/1.93/bin/ \\\n")
            f.write("--novoDir /sc/arion/projects/zhuj05a/Wenhui/packages/novoalign/novocraft/ \\\n")
            f.write("--HLAexonLoc /sc/arion/projects/zhuj05a/Wenhui/packages/LOHHLA/lohhla/data/hla.dat >& log\n")
            f.close()
            fxx.write("bsub < %s\n" % lsf)
    fxx.close()
        
         
 
if __name__=="__main__":
    
    datadir="/sc/arion/projects/zhuj05a/Wenhui/bladder/tcga_blca/"
    workdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/workdir/alleleCount/normal_solid/"
    outdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/output/alleleCount/normal_solid/"
    bindir="/sc/arion/projects/zhuj05a/Wenhui/packages/alleleCount/bin/"
    samplelist="/sc/arion/projects/zhuj05a/Wenhui/bladder/tcga_blca_clinical/GDCdata/tcga_blca_clincial_x.csv"          
    poslist="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/pos_dir/poslist" 
    #run_alleleCount(datadir,workdir,outdir,bindir,samplelist,poslist)
    count_dir=outdir
    dep_dir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/output/depth_count/normal_solid/"
    workdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/workdir/merge_count/normal_solid/"
    errorinputfile="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/samp_with_3bam_normal_solid"
    #merge_count(count_dir,dep_dir,samplelist,workdir,errorinputfile)
       
    act_dir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/output/ascat_input/normal_solid/"
    workdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/workdir/act_dir/normal_solid/"
    cov_cut=30
    #run_gen_input(dep_dir,act_dir,samplelist,workdir,cov_cut)                
    
    datadir=act_dir
    workdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/workdir/run_asact/normal_solid/"
    data_cnv="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/output/asact_out/normal_solid/"
    #run_ascat(samplelist,datadir,workdir,data_cnv) 
    
    datadir="/sc/arion/projects/zhuj05a/Wenhui/bladder/tcga_blca/"
    workdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/workdir/polysolver/normal_solid/"
    data_allele="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/output/polysolver/normal_solid/" 
    #run_polysolver(datadir,samplelist,workdir,data_allele)

    datadir="/sc/arion/projects/zhuj05a/Wenhui/bladder/tcga_blca/"
    workdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/workdir/run_lohhla/normal_solid/"
    outdir="/sc/arion/projects/zhuj05a/Wenhui/bladder/lohhla/output/run_lohhla/normal_solid/"
    run_lohhla(samplelist,datadir,data_cnv,data_allele,workdir,outdir)
