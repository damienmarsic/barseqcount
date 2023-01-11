#!/usr/bin/env python
__version__='0.1.1'
last_update='2023-01-12'
author='Damien Marsic, damien.marsic@aliyun.com'
license='GNU General Public v3 (GPLv3)'

import dmbiolib as dbl
import argparse,sys,os,itertools,regex,math
from glob import glob
import numpy as np
from collections import defaultdict
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

script=os.path.basename(__file__).split('.')[0]

def main():
    parser=argparse.ArgumentParser(description="Analysis of DNA barcode sequencing experiments. For full documentation, visit: https://"+script+".readthedocs.io")
    parser.add_argument('-v','--version',nargs=0,action=override(version),help="Display version")
    subparser=parser.add_subparsers(dest='command',required=True)
    parser_a=subparser.add_parser('count',help="Count barcodes from read files")
    parser_a.add_argument('-c','--configuration_file',default=script+'_count.conf',type=str,help='Configuration file for the '+script+' count program (default: '+script+'_count.conf), will be created if absent')
    parser_a.add_argument('-n','--new',default=False,action='store_true',help="Create new configuration file and rename existing one")
    parser_c=subparser.add_parser('analyze',help="Analyze data")
    parser_c.add_argument('-c','--configuration_file',default=script+'_analyze.conf',type=str,help="Configuration file for the "+script+" analyze program (default: "+script+"_analyze.conf), will be created if absent")
    parser_c.add_argument('-n','--new',default=False,action='store_true',help="Create new configuration file and rename existing one")
    parser_c.add_argument('-f','--file_format',type=str,default='Single multipage pdf',help="Save each figure in separate file with choice of format instead of the default single multipage pdf file. Choices: svg, png, jpg, pdf, ps, eps, pgf, raw, rgba, tif")
    args=parser.parse_args()
    if args.command=='count':
        count(args)
    if args.command=='analyze':
        analyze(args)

def count(args):
    fname=args.configuration_file
    if args.new:
        dbl.rename(fname)
    if args.new or not dbl.check_file(fname,False):
        countconf(fname,args)
        return
    read=''
    proj=''
    rfiles=[]
    templ=''
    tfile=False
    bc={}
    defn={}
    i=0
    probe=5
    fail=''
    f=open(fname,'r')
    for line in f:
        ln=line.strip()
        if ln[:9]=='# PROJECT':
            read='proj'
        if ln[:11]=='# READ FILE':
            read='rfiles'
        if ln[:10]=='# TEMPLATE':
            read='templ'
        if ln[:18]=='# PRIMERS/BARCODES':
            read='bc'
        if (not ln or ln[:2]=='# ') and read=='bc' and i in [bc[k][1] for k in bc]:
            i+=1
        if ln[:13]=='# DEFINITIONS':
            read='defn'
        if ln[:7]=='# PROBE':
            read='probe'
        if ln[:3]=='===' or ln[:2]=='# ' or ln[:13]=='Instructions:' or not ln:
            continue
        if read=='proj':
            proj=ln
        if read=='rfiles':
            x=ln.split()
            fail+=dbl.check_read_file(x[-1])
            if len(x)>2:
                fail+='\n  Too many items per line under READ FILES! Each line must contain a prefix followed by a single file name (merged reads if paired-end sequencing), separated by space or tab!'
            if len(x)==1:
                z=dbl.prefix([x[0]])[0]
                x.insert(0,x[0][:z])
            if x[0] in [k[0] for k in rfiles]:
                fail+='\n  Duplicate prefix '+x[0]+' found under READ FILES! Each line must contain a different prefix!'
            else:
                rfiles.append(x)
        if read=='templ':
            if not templ:
                x,y=dbl.getfasta(ln,'atgc'+dbl.ambiguous,'atgcn',False)
                if 'could not' in y:
                    templ=ln
                else:
                    templ=list(x.values())[0]
                    fail+=y
                    tfile=True
            elif not tfile and '.fa' not in ln:
                templ+=ln
            else:
                fail+='\n  A single item (either file name or nucleotide sequence) must be entered in the TEMPLATE SEQUENCE section !'
        if read=='bc':
            x=ln.split()
            if (len(x)==2 or (len(x)==3 and x[2] in ('+','-'))) and dbl.check_seq(x[1],'atgc','atgc')==(1,1):
                x[1]=x[1].lower()
                if x[-1]=='-':
                    x[1]=dbl.revcomp(x[1])
                if x[0] in bc:
                    fail+='\n  Duplicate barcode name: '+ln
                    continue
                if x[1] in [bc[k][0] for k in bc if k[1]==i]:
                    fail+='\n  Duplicate barcode: '+ln
                    continue
                bc[x[0]]=[x[1],i]
            else:
                fail+='\n  Incorrect barcode definition: '+ln
        if read=='defn':
            x=ln.split()
            if x[0] in defn:
                fail+='\n  Duplicate definition name: '+ln
                continue
            if x[1:] in defn.values():
                fail+='\n  Duplicate definition: '+ln
                continue
            defn[x[0]]=x[1:]
        if read=='probe':
            if not ln.isdigit() or int(ln)<1 or int(ln)>50:
                fail+='\n  Probe length must be an integer between 1 and 50'
            else:
                probe=int(ln)
    f.close()
    rname=proj+'_count_report.txt'
    dbl.rename(rname)
    print('\n  Checking configuration... ',end='')
    if not proj:
        fail+='\n  Project name is missing!'
    if not rfiles:
        fail+='\n  Read files are missing!'
    if templ:
        t,req=dbl.check_seq(templ,'atgcn','atgc')
        if not t:
            fail+='\n  Template sequence contains invalid characters!'
        if not req:
            fail+='\n  Template sequence does not contain any nucleotide!'
    else:
        fail+='\n  Template sequence is missing!'
    if not bc:
        print('  Warning! Primers/barcodes are missing. No error correction can be performed. Barcode sequences will be dispplayed in results.\n')
    if bc and not defn:
        print('  Warning! Sample definitions are missing. Barcode names will be displayed in results.\n')
    if defn:
        x=set([k for n in defn for k in defn[n]])
        for n in sorted(bc):
            if n not in x:
                del bc[n]
        for n in x:
            if n not in bc:
                fail+='\n  Unknown barcode: '+n
    for n in bc:
        if len(bc[n][0])<10:
            continue
        l0=bc[n][0]
        l1=dbl.revcomp(bc[n][0])
        a0,x0,b0,y0=maxmatch(l0,templ,10)
        a1,x1,b1,y1=maxmatch(l1,templ,10)
        if a0+b0==a1+b1==0:
            fail+='\n  Primer does not match the template: '+n+' '+bc[n][0]
            continue
        elif a0+b0>a1+b1:
            l,a,x,b,y=l0,a0,x0,b0,y0
        else:
            l,a,x,b,y=l1,a1,x1,b1,y1
        bc[n][0]=l
        if a and not b:
            templ=templ[x:]
            bc[n]=[bc[n][0],-1,a,0]
        if b and not a:
            if y+b<len(templ):
                templ=templ[:y+b]
            bc[n]=[bc[n][0],-1,0,b]
        if a and b:
            bc[n]=[bc[n][0],-1,a,b]
    w=[k for k in bc if len(bc[k][0])>=10]
    y=min([bc[k][-2] for k in w if bc[k][-2]>0])
    z=min([bc[k][-1] for k in w if bc[k][-1]>0])
    x0=set([bc[k][0][-y:] for k in w if y])
    x1=set([bc[k][0][:z] for k in w if z])
    for n in x0:
        a=min([bc[k][2] for k in w if bc[k][0].endswith(n)])
        for m in [k for k in w if bc[k][0].endswith(n)]:
            bc[m][2]=a
    for n in x1:
        b=min([bc[k][3] for k in w if bc[k][0].startswith(n)])
        for m in [k for k in w if bc[k][0].startswith(n)]:
            bc[m][3]=b
    x=set([(bc[k][0][-bc[k][2]:],bc[k][0][:bc[k][3]],len(bc[k][0])-bc[k][2]-bc[k][3]) for k in w if bc[k][-1] and bc[k][-2]])
    for n in x:
        a=templ.find(n[0])
        b=templ.find(n[1])+len(n[1])
        templ=templ[:b]+'n'*n[2]+templ[a:]
    w=[k for k in bc if not bc[k][-1] and len(bc[k][0])>=10]
    for i in range(1,len(bc[w[0]])):
        if len(set([bc[k][0][:i] for k in w]))>1:
            i-=1
            break
    x=bc[w[0]]
    templ=x[0][:i]+'n'*(len(x[0])-i-x[-2])+templ
    for n in w:
        bc[n][1]=bc[n][-1]=i
    w=[k for k in bc if not bc[k][-2] and len(bc[k][0])>=10]
    for i in range(1,len(bc[w[0]])):
        if len(set([bc[k][0][-i:] for k in w]))>1:
            i-=1
            break
    x=bc[w[0]]
    templ+='n'*(len(x[0])-x[-1]-i)+x[0][-i:]
    for n in w:
        bc[n][-2]=i
        bc[n][1]=templ.find(bc[n][0][:bc[n][-1]])+bc[n][-1]
    for n in [k for k in bc if bc[k][1]==-1]:
        bc[n][1]=templ.find(bc[n][0][:bc[n][-1]])+bc[n][-1]
    bcr=dbl.find_ambiguous(templ)
    x=[k for k in bc if len(bc[k][0])<10]
    for n in x:
        bc[n].append(len(bc[n][0]))
    y=sorted(set([(bc[k][1],bc[k][2]) for k in x]))
    z=[k for k in bcr if k not in [n[1] for n in bc if len(bc[n])==4]]
    if list(bcr.values())==[k[1] for k in y]:
        for i in range(len(bcr)):
            for n in [k for k in x if bc[k][-2:]==list(y[i])]:
                bc[n]=[bc[n][0],list(bcr.keys())[i]]
    elif len(y)==1==len(z) and bcr[z[0]]==y[0][1]:
        for n in x:
            bc[n]=[bc[n][0],z[0]]
    x=[k for k in bc if len(bc[k])==3]
    if x:
        fail+='\n  Some barcodes could not be assigned to the template sequence: '+', '.join(x)
    else:
        for n in [k for k in bc if len(bc[k])==4]:
            bc[n]=[bc[n][0][bc[n][-1]:-bc[n][-2]],bc[n][1]]
    BC={}
    for n in bcr:
        if not bc:
            break
        BC[n]={bc[k][0]:k for k in bc if bc[k][1]==n}
        x=[k[0] for k in bc.values() if k[1]==n]
        if not x:
            bcr[n]=[bcr[n],False,False]
            continue
        if len(set(x))!=len(x):
            fail+='\n  Barcode sequence duplicate found at position '+str(n)+'!'
        bcr[n]=[bcr[n],dbl.diff(x)>2]
        y=set([len(dbl.compress(k)) for k in x])
        z=set([len(k) for k in x])
        if len(y)==1 and y==z and (not n or templ[n-1] not in [k[0] for k in x]) and (n+bcr[n][0]==len(templ) or templ[n+bcr[n][0]] not in [[k[-1] for k in x]]):
            bcr[n].append(True)
        else:
            bcr[n].append(False)
    dr=sorted(set([tuple(set([bc[n][1] for n in k])) for k in defn.values()]))
    for n in bcr:
        if len([m for k in dr for m in k if m==n])>1:
            fail+='\n  Inconsistency in sample/variant definitions!'
    DEF={n:{} for n in dr}
    for n in dr:
        x=[]
        for m in n:
            x.append(list(BC[m].values()))
        for m in defn:
            y=[]
            for a in x:
                y.append([])
                for b in defn[m]:
                    if b in a:
                        y[-1].append(b)
            if not y[-1]:
                continue
            z=list(itertools.product(*y))
            for a in z:
                DEF[n][a]=m
    if fail:
        print('\n'+fail+'\n')
        sys.exit()
    fail=''
    print('OK\n\n  Checking read files...    ',end='')
    for j in range(len(rfiles)):
        nr,fail=dbl.readcount(rfiles[j][1],fail)
        rfiles[j].append(nr)
        if nr<100:
            fail+='\n  Number of reads in '+rfiles[j][0]+' is too low!'
        f,step=dbl.initreadfile(rfiles[j][1])
        c=0
        a=0
        b=0
        lpmet=dbl.revcomp(templ)
        while c<200:
            l,f,c,_=dbl.getread(f,step,c)
            d=min(30,len(l))
            if dbl.match(l[:d],templ[:d]):
                a+=1
                continue
            if dbl.match(l[:d],lpmet[:d]):
                b+=1
        f.close()
        x='+'
        if b and not a:
            x='-'
        elif b and a:
            x='Both'
        elif not a and not b:
            x='None'
            fail+='\n  Read file '+rfiles[j][1]+' does not contain sequences matching the template sequence!'
        rfiles[j].append(x)
    r=open(rname,'w')
    if fail:
        dbl.pr2(r,'Problems found!\n'+fail+'\n')
    else:
        print('OK\n')
    dbl.pr2(r,'  Read file prefix         Read file                     Number of reads    Read orientation')
    for n in rfiles:
        dbl.pr2(r,'  '+n[0].ljust(25)+n[1].ljust(30)+f'{n[2]:,}'.rjust(15)+n[3].center(24))
    dbl.pr2(r,'')
    if fail:
        r.close()
        sys.exit()
    dbl.pr2(r,'  Template sequence:\n'+dbl.format_dna(templ,2,80,10))
    if bc:
        x={k:[(n,bc[n][0]) for n in bc if bc[n][1]==k] for k in bcr}
        z={n:max([len(str(k)) for k in x[n]])-5 for n in x}
        w=max([len(x[n]) for n in x])
        dbl.pr2(r,'  Barcodes:\n  '+''.join([str(k).ljust(z[k]) for k in x]))
        for i in range(w):
            for n in x:
                if len(x[n])<w:
                    x[n].append(('',''))
            dbl.pr2(r,'  '+''.join([(x[k][i][0]+' '+x[k][i][1]).ljust(z[k]) for k in x]))   #####  add blank spaces or .ljust !!!!
        dbl.pr2(r,'\n  Error correction:\n  Position  Single substitutions  Indels within homopolymer')
        for n in bcr:
            dbl.pr2(r,str(n).rjust(10)+str(bcr[n][1]).center(22)+str(bcr[n][2]).center(27))
    if defn:
        x={**{k:str(k[0]) for k in DEF if len(k)==1},**{k:str(k) for k in DEF if len(k)>1}}
        z=[max(11,max([len(k) for k in x.values()]))]
        z.append(max(12,max([len(', '.join(k))+2 for n in DEF for k in DEF[n]])))
        z.append(max(12,max([len(k)+2 for n in DEF for k in DEF[n].values()])))
        dbl.pr2(r,'\n  Definitions:\n  '+'Position(s)'.center(z[0])+' Barcode(s) '.center(z[1])+' Definition '.center(z[2]))
        for n in DEF:
            for m in DEF[n]:
                dbl.pr2(r,'  '+x[n].center(z[0])+(', '.join(m)).center(z[1])+DEF[n][m].center(z[2]))
    for i in bcr:
        bcr[i].append(max(0,i-probe))
        bcr[i].append(i+bcr[i][0])
        bcr[i].append(min(len(templ),i+bcr[i][0]+probe))
    ctempl=dbl.compress(templ)
    cbcr=dbl.find_ambiguous(ctempl)
    for i in cbcr:
        cbcr[i]=[cbcr[i]]
        cbcr[i].append(max(0,i-probe))
        cbcr[i].append(i+cbcr[i][0])
        cbcr[i].append(min(len(ctempl),i+cbcr[i][0]+probe))
    counts=defaultdict(int)
    ec=[0,0,0,0,0]  # alternate position, compressed mode, compressed mode + alternate position, barcode single substitution, corrected reads
    C=0
    for rfile in rfiles:
        print()
        f,step=dbl.initreadfile(rfile[1])
        t='Processing reads from '+rfile[1]+'...'
        show=dbl.progress_start(rfile[2],t)
        c=0
        compr=True in [bcr[k][2] for k in bcr]
        while True:
            X1=X2=None
            p1=p2=0
            l,f,c,_=dbl.getread(f,step,c)
            if not l:
                break
            cl=''
            if compr:
                cl=dbl.compress(l)
            if rfile[3]!='-':
                X1,p1,ec1=find_bc(l,templ,bcr,cl,ctempl,cbcr)
            if rfile[3]!='+':
                X2,p2,ec2=find_bc(dbl.revcomp(l),templ,bcr,cl,ctempl,cbcr)
            if X1 and (not X2 or (X2 and p2>p1)):
                X=X1
                EC=ec1
                p=p1
            elif X2 and (not X1 or (X1 and p1>p2)):
                X=X2
                EC=ec2
                p=p2
            else:
                continue
            C+=1
            for i in range(len(EC)):
                ec[i]+=EC[i]
            for i in X:
                if bcr[i][1] and X[i] not in BC[i]:
                    for n in BC[i]:
                        if dbl.diff((X[i],n))==1:
                            X[i]=n
                            ec[3]+=1
                            p+=1
                            break
                if X[i] in BC[i]:
                    X[i]=BC[i][X[i]]
            if p:
                ec[4]+=1
            Y=''
            if len(rfiles)>1:
                Y=rfile[0]
            for n in dr:
                x=[]
                for m in n:
                    x.append(X[m])
                x=tuple(x)
                if Y:
                    Y+=','
                if x in DEF[n]:
                    Y+=DEF[n][x]
                else:
                    Y+=','.join(x)
            counts[Y]+=1
            dbl.progress_check(c,show,t)
        dbl.progress_end()
        f.close()
    x=[DEF[k].values() for k in DEF]
    y=0
    for n in counts:
        for m in n.split(','):
            if m not in x:
                y+=1
                break
    x=sum([k[2] for k in rfiles])
    dbl.pr2(r,('\n  Reads processed:').ljust(45)+f'{x:,}'.rjust(15))
    dbl.pr2(r,('  Successful reads (all barcodes detected):').ljust(45)+f'{C:,}'.rjust(15)+' ('+f'{C/x*100:.2f}'.rjust(6)+'% of total reads)')
    dbl.pr2(r,('  Error-corrected reads:').ljust(45)+f'{ec[4]:,}'.rjust(15)+' ('+f'{ec[4]/C*100:.2f}'.rjust(6)+'% of successful reads)')
    dbl.pr2(r,('  Reads with all barcodes identified:').ljust(45)+f'{C-y:,}'.rjust(15)+' ('+f'{(C-y)/C*100:.2f}'.rjust(6)+'% of successful reads)')
    dbl.pr2(r,'\n  Error-corrected barcodes in successful reads:')
    dbl.pr2(r,'  Alternate position:'.ljust(40)+f'{ec[0]:,}'.rjust(15))
    dbl.pr2(r,'  Indel within homopolymer in probe:'.ljust(40)+f'{ec[1]:,}'.rjust(15))
    dbl.pr2(r,'  Alternate position + indel:'.ljust(40)+f'{ec[2]:,}'.rjust(15))
    dbl.pr2(r,'  Nucleotide substitution:'.ljust(40)+f'{ec[3]:,}'.rjust(15))
    dbl.csv_write(proj+'_count.csv',None,counts,None,'Barcode distribution',r)
    r.close()
    print('\n  Report was saved into file: '+rname+'\n')

def analyze(args):
    cf=args.configuration_file
    if args.new:
        dbl.rename(cf)
    if args.new or not dbl.check_file(cf,False):
        anaconf(cf,args)
        return
    format=dbl.check_plot_format(args.file_format)
    print('\n  Checking configuration... ',end='')
    f=open(cf,'r')
    fail=''
    read=''
    bcount=''
    V=[]
    mix=''
    thr=0
    S={}
    gunit=''
    gtit={}
    eunit=''
    etit={}
    comb=[]
    titles=[]
    xaxis=None
    showind=True
    showerr=None
    bcmap=[]
    hcmap=[]
    for line in f:
        ln=line.strip()
        if ln[:1]=='#':
            read=''
        if ln[:9]=='# BARCODE':
            read='bcount'
        if ln[:10]=='# VARIANTS':
            read='variants'
        if ln[:13]=='# VARIANT MIX':
            read='mix'
        if ln[:11]=='# THRESHOLD':
            read='thr'
        if ln[:9]=='# SAMPLES':
            read='samples'
        if ln[:15]=='# GLOBAL GENOME':
            read='gtit'
        if ln[:19]=='# GLOBAL EXPRESSION':
            read='etit'
        if ln[:9]=='# COMBINE':
            read='comb'
        if read=='comb' and not ln:
            comb.append({})
            titles.append([])
        if ln[:8]=='# X-AXIS':
            read='xaxis'
        if ln[:10]=='# SHOW IND':
            read='showind'
        if ln[:10]=='# SHOW ERR':
            read='showerr'
        if ln[:10]=='# BAR PLOT':
            read='bcmap'
        if ln[:10]=='# HEAT MAP':
            read='hcmap'
        if ln[:3]=='===' or ln[:2]=='# ' or ln[:13]=='Instructions:' or not ln:
            continue
        if read=='bcount':
            if bcount:
                fail+='\n  Only one file should be listed under BARCODE COUNT FILE!'
            elif dbl.check_file(ln,False):
                proj=ln[:ln.rfind('_count')]
                _,bcount=dbl.csv_read(ln,False,None)
            else:
                fail+='\n  Barcode count file '+ln+' not found!'
        if read=='variants':
            if not ln in V:
                V.append(ln)
        if read=='mix':
            if mix:
                fail+='\n  Only one item should be present under VARIANT MIX!'
            else:
                mix=ln
        if read=='thr':
            x=ln.split()
            if x[0].replace('.','').replace('%','').isdigit():
                x=float(x[0].replace('%',''))
                if not 0<x<100 and not thr:
                    fail+='\n  Threshold nust be a number (integer or decimal) larger than 0 and lower than 100!'
                    continue
                if thr:
                    fail+='\n  Only one value should be present under THRESHOLD!'
                    continue
                thr=x
        if read=='samples':
            if not ln in S:
                S[ln]=set(ln.split(','))
        if read in ('gtit','etit'):
            x=ln.split()
            if len(x)==1 or not x[1].replace('.','').isdigit():
                fail+='\n  Titer value missing: '+ln
            elif x[0] not in S:
                fail+='\n  Unknown sample name: '+x[0]
            elif read=='gtit' and gunit and gunit!=' '.join(x[2:]):
                fail+='\n  Genome titers must all have the same unit!'
            elif read=='etit' and eunit and eunit!=' '.join(x[2:]):
                fail+='\n  Expression titers must all have the same unit!'
            else:
                if read=='gtit':
                    if len(x)>2:
                        gunit=' '.join(x[2:])
                    gtit[x[0]]=float(x[1])
                if read=='etit':
                    if len(x)>2:
                        eunit=' '.join(x[2:])
                    etit[x[0]]=float(x[1])
        if read=='comb':
            x=ln.split()
            if not comb or (comb[-1] and not x[-1] in S):
                comb.append({})
                titles.append([])
            if not comb[-1] and not any([k in S for k in x[1:]]):
                if not titles[-1]:
                    titles[-1].extend(x)
                else:
                    fail+='\n  Each plot can only have one title: '+' '.join(titles[-1])+', '+ln
                continue
            if len(x)>1 and all([k in S for k in x[1:]]):
                comb[-1][x[0]]=x[1:]
            elif len(x)==1:
                fail+='\n  In each group under COMBINE DATA, each line must contain a label followed by sample name(s) separated by spaces or tabs.'
            else:
                fail+='\n  Unknown sample in this list: '+ln
        if read=='bcmap':
            bcmap.append(ln)
        if read=='hcmap':
            hcmap.append(ln)
        ln=ln.lower()
        if read=='xaxis':
            if xaxis:
                fail+='\n  Only one item should be present under X-AXIS!'
            elif ln[0]=='v':
                xaxis=True
            elif ln[0]=='s':
                xaxis=False
            else:
                fail+='\n  X-AXIS must be either Variants or Samples!'
        if read=='showind':
            if ln[0] in ('p','.'):
                showind='k.'
            elif ln[0] in ('c','o'):
                showind='ko'
            elif ln[0] in ('s','*'):
                showind='k*'
            elif ln[0] in ('f','n'):
                showind=None
            else:
                fail+='\n  SHOW INDIVIDUL DATA POINTS can only be Points, Circles, Stars or None!'
        if read=='showerr':
            if ln[0]=='r':
                showerr='r'
            elif 'se' in ln or 'er' in ln:
                showerr='se'
            elif 'sd' in ln or 'de' in ln:
                showerr='sd'
            elif ln[0]=='n':
                showerr=None
            else:
                fail+='\n  SHOW ERROR BARS can only be Range, Standard error, Standard deviation or None!'
    for n in comb:
        if not n:
            continue
        if len(set([len(k) for k in n.values()]))>1:
            fail+='\n  Within groups under COMBINE DATA, all lines must contain the same number of samples!'
        x=len(list(n.values())[0])
        y=len(titles[comb.index(n)])
        if x>1 and (x>y or x<y-1):
            fail+='\n  Mismatch between number of identifiers and number of samples per line under COMBINE DATA!'
        if x>1 and x==y:
            titles[comb.index(n)].insert(0,'')
    if xaxis==None:
        xaxis=True
    if not bcmap or (len(set(bcmap))==1 and bcmap[0].lower()=='default'):
        bcmap=['turbo_r','viridis','viridis','turbo_r']
    elif len(bcmap)!=4:
        fail+='\n  There should be exactly 4 color maps for bar plots!'
    if not hcmap or (len(set(hcmap))==1 and hcmap[0].lower()=='default'):
        hcmap=['seismic','hot_r']
    elif len(hcmap)!=2:
        fail+='\n  There should be exactly 2 color maps for heat maps!'
    if fail:
        print('Problems found!\n'+fail+'\n')
        sys.exit()
    else:
        print('OK')
    for i in range(4):
        if bcmap[i].lower()=='default':
            bcmap[i]=['turbo_r','viridis','viridis','turbo_r'][i]
    for i in range(2):
        if hcmap[i].lower()=='default':
            hcmap[i]=['seismic','hot_r'][i]
    pre=proj+'_'+str(thr)+'_'+str(int(xaxis))+'_'
    if format:
        mppdf=''
        print()
    else:
        fname=pre+'figs.pdf'
        dbl.rename(fname)
        mppdf=PdfPages(fname)
    V.sort()
    M={}
    if mix:
        M={n:0 for n in V}
        m=set(mix.split(','))
        for n in bcount:
            for q in n[:-1]:
                if q in V:
                    x=set(k for k in n[:-1] if k!=q)
                    break
            else:
                continue
            if x==m:
                M[q]=int(n[-1])
        x=sum(M.values())
        q={n:M[n]/x for n in M}
        z=1/len(M)
        y=[k/z*100 for k in q.values()]
        colors,fig=dbl.plot_start(bcmap[0],len(M),'Variant mix composition')
        plt.bar(M.keys(),y,color=colors.colors)
        plt.xticks(rotation=90)
        plt.ylabel('% of equimolar frequency')
        plt.margins(x=0.015)
        a=pre+'variant-mix-composition'
        dbl.plot_end(fig,a,format,mppdf)
        dbl.csv_write(a+'.csv',list(M.keys()),y,None,'Variant mix composition',None)
        V=[V[i] for i in range(len(V)) if y[i]>=thr]
        if len(V)<len(M):
            print('\n  The following variants will be excluded from further analysis because their mix frequencies are below the threshold:')
            print('  '+'\n  '.join([k for k in M if not k in V])+'\n')
        z=sum([M[n] for n in V])
        M={n:M[n]/z for n in V}
    BC={k:{n:1 for n in V} for k in S}
    for n in bcount:
        for q in n[:-1]:
            if q in V:
                x=set(k for k in n[:-1] if k!=q)
                break
        else:
            continue
        if x in S.values():
            y=[k for k in S if S[k]==x][0]
            BC[y][q]=n[-1]
    del bcount
    x=[]
    y=[]
    if mix:
        x.append(mix)
        y.append(z)
    x.extend(list(S.keys()))
    y.extend([sum(BC[k].values()) for k in BC])
    colors,fig=dbl.plot_start(bcmap[1],len(x),'Global read count per sample')
    plt.bar(x,y,color=colors.colors)
    plt.xticks(rotation=90)
    plt.ylabel('Read count')
    plt.margins(x=0.015)
    a=pre+'read-count-per-sample'
    dbl.plot_end(fig,a,format,mppdf)
    dbl.csv_write(a+'.csv',x,y,None,'Global read cont per sample',None)
    for n in BC:
        for m in BC[n]:
            z=BC[n][m]/y[x.index(n)]
            if mix:
                BC[n][m]=z/M[m]
            else:
                BC[n][m]=z*len(V)
    colors,fig=dbl.plot_start(None,None,'Global variant enrichment')
    x=[[math.log10(n) for n in list(BC[k].values())] for k in BC]
    a=(list(BC.keys()),V)
    if not xaxis:
        x=list(zip(*x))
    plt.imshow(x,aspect='auto',cmap=hcmap[0])
    lim=max(max([n for m in x for n in m]),abs(min([n for m in x for n in m])))
    plt.clim(-lim, lim)
    plt.colorbar(shrink=0.7,pad=0.015,label='Log(Enrichment factor)')
    plt.xticks(range(len(a[xaxis])),a[xaxis],rotation=90)
    plt.yticks(range(len(a[not xaxis])),a[not xaxis])
    b=pre+'global_enrichment'
    dbl.plot_end(fig,b,format,mppdf)
    dbl.csv_write(b+'.csv',a[not xaxis],x,a[xaxis],'Global enrichment',None)
    a=[]
    if gtit:
        a.append(gtit)
    if etit:
        a.append(etit)
    for n in BC:
        for q in a:
            if n in q:
                x=sum(BC[n].values())
                for m in BC[n]:
                    BC[n][m]=BC[n][m]/x*q[n]
                break
    for i in range(len(comb)):
        if not comb[i]:
            continue
        if all([k in gtit for n in comb[i].values() for k in n]):
            unit=gunit
            a=gtit
        elif all([k in etit for n in comb[i].values() for k in n]):
            unit=eunit
            a=etit
        else:
            continue
        locs=list(range(len(titles[i][1:])))
        labels=list(comb[i].keys())
        colors,fig=dbl.plot_start(bcmap[2],len(labels),'Global '+titles[i][0].lower()+' biodistribution')
        z=0.8/len(labels)
        w=z*0.95
        for k in range(len(labels)):
            plt.bar([n+z*k for n in locs],[a[comb[i][labels[k]][j]] for j in range(len(comb[i][labels[k]]))],width=w,label=labels[k],color=colors.colors[k])
        plt.xticks([n+z*k/2 for n in locs],titles[i][1:])
        plt.ylabel(unit)
        plt.legend()
        dbl.plot_end(fig,pre+'global_'+titles[i][0].lower()+'_biodistribution',format,mppdf)
    for i in range(len(comb)):
        if not comb[i]:
            continue
        if all([k in gtit for n in comb[i].values() for k in n]):
            unit=gunit
            a=gtit
        elif all([k in etit for n in comb[i].values() for k in n]):
            unit=eunit
            a=etit
        else:
            unit='Enrichment factor'
        colors,fig=dbl.plot_start(None,None,titles[i][0]+' biodistribution')
        X=[[dbl.mean([BC[n][m] for n in comb[i][k]]) for m in V] for k in comb[i]]
        if showerr=='r':
            yerr0=[[min([BC[n][m] for n in comb[i][k]]) for m in V] for k in comb[i]]
            yerr1=[[max([BC[n][m] for n in comb[i][k]]) for m in V] for k in comb[i]]
        elif showerr=='se':
            yerr0=[[np.std([BC[n][m] for n in comb[i][k]],ddof=1)/np.sqrt(np.size([BC[n][m] for n in comb[i][k]])) for m in V] for k in comb[i]]
        elif showerr=='sd':
            yerr0=[[np.std([BC[n][m] for n in comb[i][k]],ddof=1) for m in V] for k in comb[i]]
        if showind:
            show=[[[BC[n][m] for n in comb[i][k]] for m in V] for k in comb[i]]
        a=(list(comb[i].keys()),V)
        if not xaxis:
            X=list(zip(*X))
        plt.imshow(X,aspect='auto',cmap=hcmap[1])
        plt.colorbar(shrink=0.7,pad=0.015,label=unit)
        plt.xticks(range(len(a[xaxis])),a[xaxis],rotation=90)
        plt.yticks(range(len(a[not xaxis])),a[not xaxis])
        dbl.plot_end(fig,pre+titles[i][0],format,mppdf)
        locs=list(range(len(a[not xaxis])))
        labels=a[xaxis]
        colors,fig=dbl.plot_start(bcmap[3],len(a[xaxis]),titles[i][0]+' biodistribution')
        z=0.8/len(labels)
        w=z*0.95 
        for j in range(len(labels)):
            x=[n+z*j for n in locs]
            h=[X[k][j] for k in range(len(a[not xaxis]))]
            plt.bar(x,[X[k][j] for k in range(len(a[not xaxis]))],width=w,label=labels[j],color=colors.colors[j])
            if showerr:
                if showerr=='r':
                    e0=[X[k][j]-yerr0[k][j] for k in range(len(a[not xaxis]))]
                    e1=[yerr1[k][j]-X[k][j] for k in range(len(a[not xaxis]))]
                else:
                    e0=[yerr0[k][j] for k in range(len(a[not xaxis]))]
                    e1=e0
                plt.errorbar(x,h,yerr=[e0,e1],fmt='none',ecolor='black',capsize=2)
            if showind:
                y=[show[k][j] for k in range(len(a[not xaxis]))]
                for x1, y1 in zip(x, y):
                    plt.plot([x1] * len(y1), y1,showind)
        plt.xticks([n+z*j/2 for n in locs],a[not xaxis])
        plt.ylabel(unit)
        plt.legend()
        dbl.plot_end(fig,pre+titles[i][0]+'_biodistribution',format,mppdf)
        dbl.csv_write(pre+titles[i][0]+'_biodistribution.csv',a[not xaxis],X,a[xaxis],titles[i][0]+' biodistribution',None)
    if not format:
        mppdf.close()
        print('\n  All figures were saved into single multipage file: '+fname+'\n')
    f.close()

def anaconf(fname,args):
    x=glob('*_count_report.txt')
    if not x:
        print('\n  Count report file not found! It must be present in the working directory!\n')
        sys.exit()
    report=max(x,key=os.path.getctime)
    f=open(report,'r')
    tlen=0
    defn=[]
    rfp=[]
    cfile=''
    a=b=c=False
    for l in f:
        l=l.strip()
        if not l:
            continue
        if 'Read file prefix' in l:
            a=True
            continue
        if 'Template sequence:' in l:
            a=False
            b=True
            continue
        if 'Barcodes:' in l:
            b=False
            continue
        if all(k in l for k in ('Position','Barcode','Definition')):
            c=True
            continue
        if 'Reads processed:' in l:
            c=False
            continue
        if 'Barcode distribution' in l:
            cfile=l.strip().split()[-1]
        if a:
            ln=l.strip().split()
            if ln:
                rfp.append(ln[0])
        if b:
            ln=l.strip().split()
            if ln and not ln[0].isdigit():
                tlen+=len(ln[0])
            continue
        if c:
            ln=l.strip().replace('(','').replace(')','').replace(',','').split()
            d=[(int(ln[0]),),ln[-1]]
            for i in range(1,len(ln)-1):
                if not ln[i].isdigit():
                    break
                d[0]+=(int(ln[i]),)
            if not d in defn:
                defn.append(d)
    f.close()
    if not defn:
        print('\n  No variant or sample definition found!\n')
        sys.exit()
    x={}
    for n in set([k[0] for k in defn]):
        x[n]=max([min(k,abs(k-tlen)) for k in n])
    v=sorted([k[1] for k in defn if x[k[0]]==max(x.values())])
    s=[]
    if len(rfp)>1:
        s.append(sorted(rfp))
    for n in set([k[0] for k in defn if k[1] not in v]):
        s.append(sorted([k[1] for k in defn if k[0]==n]))
    s=[','.join(k) for k in itertools.product(*s)]
    mix=''
    m=[k for k in s if 'mix' in k.lower()]
    if len(m)!=1:
        m=[]
        for q in s:
            a=q.replace('_','-').replace(',','-').split('-')
            b=[k.replace('_','-').replace(',','-').split('-') for k in s if k!=q]
            for k in b:
                for m in a:
                    if m in k:
                        break
                else:
                    continue
                break
            else:
                m.append(q)
    if len(m)==1:
        mix=m[0]
        s=[k for k in s if k!=mix]
    f=open(fname,'w')
    f.write('=== '+script.upper()+' '+args.command.upper()+' CONFIGURATION FILE ===\n\n')
    f.write('# BARCODE COUNT FILE:\n\n'+cfile+'\n\n')
    f.write('# VARIANTS:\n\n'+'\n'.join(v)+'\n\n')
    f.write('# VARIANT MIX:\n\n'+mix+'\n\n')
    f.write('# THRESHOLD:\nInstructions: Threshold in % of equimolar distribution in the variant mix. Variants with a frequency lower than the threshold in the mix will be excluded from analysis.\n\n10\n\n')
    f.write('# SAMPLES:\n\n'+'\n'.join(s)+'\n\n')
    dna=[k for k in s if any([n in [m.lower() for m in k.replace('_','-').replace(',','-').split('-')] for n in ('dna','gdna','genome')])]
    rna=[k for k in s if any([n in [m.lower() for m in k.replace('_','-').replace(',','-').split('-')] for n in ('rna','cdna','expression')])]
    other=[k for k in s if (k not in dna and k not in rna)]
    f.write('# GLOBAL GENOME TITERS:\nInstructions: Enter the names of all samples for which you have a global variant genome relative to host genome titer, one per line, with the name followed by the titer value and the unit, separated by tabs or blank spaces. The unit must be the same for all samples, and can be entered once (will be applied to all other samples) or more.\n\n'+'\t\tvg/cell\n'.join(dna)+'\t\tvg/cell\n\n')
    f.write('# GLOBAL EXPRESSION TITERS:\nInstructions: Enter the names of all samples for which you have a global variant expression relative to host housekeeping gene expression level titer, one per line, with the name followed by the titer value and the unit, separated by tabs or blank spaces. The unit must be the same for all samples, and can be entered once (will be applied to all other samples) or more.\n\n'+'\t\t% control\n'.join(rna)+'\t\t% control\n\n')
    f.write('# COMBINE BIOLOGICAL REPLICATES:\nInstructions: Each group of lines separated by empty lines will be a separate plot. Multiple samples within lines will be averaged. First line in each group is the figure title followed by sample identifiers (separated by space or tab). Each subsequent line: label followed by sample name(s) separated by tabs or spaces.\n\n')
    for n in (dna,rna,other):
        if not n:
            continue
        p=[k.replace('_','-').replace(',','-').split('-') for k in n]
        x=[[k[i] for i in range(len(k)) if len(set([m[i] for m in p]))>1] for k in p]
        w=[[k[i] for i in range(len(k)) if all([m[i].isdigit() for m in p])] for k in p]
        q=min([len(k) for k in w])
        if not q:
            w=[[k[i] for i in range(len(k)) if all([m[i][-1].isdigit() for m in x])] for k in x]
            q=min([len(k) for k in w])
        if q:
            w=[k[-1] for k in w]
        x=[[x[i][j] for j in range(len(x[i])) if x[i][j]!=w[i]] for i in range(len(x))]
        y=''
        if n==dna:
            y='Genome'
        elif n==rna:
            y='Expression'
        z=sorted(set(['-'.join(k) for k in x]))
        if w[0]:
            w=sorted(set(w))
        else:
            w=[]
        if len(w)==x.count(z[0].split('-')):
            y+='\t'+'\t'.join(w)
        f.write(y+'\n')
        for m in z:
            f.write(m+'\t'+'\t'.join(sorted([n[i] for i in range(len(n)) if all([q in p[i] for q in m.split('-')])]))+'\n')
        f.write('\n')
    f.write('# X-AXIS:\nInstructions: either variants or samples can be on the x-axis.\n\nVariants\n\n')
    f.write('# SHOW INDIVIDUAL DATA POINTS:\nInstructions: on bar plots of combined data, superimpose individual data points (points / circles / stars / none).\n\npoints\n\n')
    f.write('# SHOW ERROR BARS:\nInstructions: display error bars on bar plots of combined data (Range / Standard error / Standard deviation / None).\n\nNone\n\n')
    f.write('# BAR PLOT COLOR MAPS:\nInstructions: list color maps to be used in bar plots, one per line for variant mix, global read count, global biodistribution, detailed biodistribution (default or leave empty for default color maps, or any compatible matplotlib color map such as turbo, viridis, plasma, inferno, magma, cividis, twilight, tab20 or their reversed versions with _r appended).\n\ndefault\ndefault\ndefault\ndefault\n\n')
    f.write('# HEAT MAP COLOR MAPS:\nInstructions: list color maps to be used in heat maps, one per line for global variant enrichment and biodistribution (default or leave empty for default color maps, or any compatible color map such as seismic, bwr, hot, RdBu or their reversed versions with _r appended).\n\ndefault\ndefault\n\n')
    f.write('=== END OF CONFIGURATION FILE ===')
    f.close()
    print('\n  Edit the file '+fname+' before running '+script+' '+args.command+' again !\n\n')

def countconf(fname,args):
    z=script+' '+args.command
    f,dirname,rfiles=dbl.conf_start(fname,z)
    f.write('# PROJECT NAME\nInstructions: Project name to be used as prefix in output file name.\n\n')
    f.write(dirname+'\n\n')
    f.write('# READ FILE(S)\nInstructions: List prefix (to be used as sample name in output files if multiple read files are present, will be created automatically if missing) and read file names, one file (single-end or merged paired-end data) per line, separated by space or tab.\n\n')
    y=[k for k in rfiles if k.count(' ')==2]
    x=[k for k in rfiles if k not in y]
    if not x and y:
        print('\n  Only unmerged paired-end read files found in the current directory! Please merge before using '+script+'!\n')
        sys.exit()
    f.write('\n'.join(x)+'\n\n')
    f.write('# TEMPLATE SEQUENCE\nInstructions: Name of the file, in FASTA format, containing the PCR template sequence, with the barcode shown as a stretch of N. Alternatively, paste the sequence itself. Primer binding regions must be included.\n\n')
    x=glob('*.f*a')
    y='\n'
    if len(x)==1:
        y=x[0]+2*y
    f.write(y)
    f.write("# PRIMERS/BARCODES\nInstructions: List all primers (name followed by sequence), whether barcoded or not, one per line (forward and reverse PCR primers used to amplify the template, as well as primers used to generate the template). Barcodes only (without primer sequence) can be entered instead, if no ambiguity. If no primer information is entered, the exact consensus read sequence must be entered in the template sequence section above, with all barcodes indicated by stretches of N, and barcodes should be grouped (each grpup surrounded by empty lines) by their corresponding region in the order the barcoded regions appear in the read sequence. Barcodes from the reverse strand need to be followed by a - sign (their reverse-complement will then be used).\n\n\n")
    f.write("# DEFINITIONS\nInstructions: Each line contains a sample/variant name followed by one or more barcode name(s), separated by a space or tab. Barcodes from different regions indicate that the sample is defined by a combination of barcodes from those regions. Multiple barcodes from the same region are allowed and indicate that a sample/variant is represented by more than one barcoode. All primers/barcodes must have a different name. In case of more than one sample/conditions is defined by a barcode or barcode combination, use a separator (- or _ or ,) between them.\n\n\n")
    f.write('# PROBE LENGTH\nInstructions: Minimum length in nt of sequences used as probes to locate barcodes (integer between 1 and 50.\n\n5\n\n')
    dbl.conf_end(f,fname,z)

def find_bc(l,templ,bcr,cl,ctempl,cbcr):
    def fb(l,templ,i,a,b,c,bcr):
        r=''
        x=[k.start()+i-a for k in regex.finditer(templ[a:i],l,overlapped=True)]
        if x:
            y=[k.start() for k in regex.finditer(templ[b:c],l,overlapped=True)]
            z=[(n,m) for n in x for m in y if m>n and m-n==bcr[i][0]]
            if z:
                if len(z)==1:
                    z=z[0]
                else:
                    q=min([abs(k[0]-i) for k in z])
                    w=[k for k in z if abs(k[0]-i)==q]
                    if len(w)==1:
                        z=w[0]
                    else:
                        z=None
                if z:
                    r=l[z[0]:z[1]]
        return r
    X={}
    p=0
    ec=[0,0,0]
    for i in bcr:
        a=bcr[i][-3]
        b=bcr[i][-2]
        c=bcr[i][-1]
        if l[a:i]==templ[a:i] and l[b:min(c,len(l))]==templ[b:c]:
            X[i]=l[i:b]
            continue
        x=fb(l,templ,i,a,b,c,bcr)
        if x:
            X[i]=x
            p+=1
            ec[0]+=1
            continue
        if bcr[i][2]:
            j=list(cbcr.keys())[list(bcr.keys()).index(i)]
            a=cbcr[j][-3]
            b=cbcr[j][-2]
            c=cbcr[j][-1]
            if cl[a:j]==ctempl[a:j] and cl[b:min(c,len(cl))]==ctempl[b:c]:
                X[i]=cl[j:b]
                p+=1
                ec[1]+=1
                continue
            x=fb(cl,ctempl,j,a,b,c,cbcr)
            if x:
                X[i]=x
                p+=2
                ec[2]+=1
                continue
        X={}
        break
    return X,p,ec

def maxmatch(sample,target,probe):
    a=b=x=y=w=z=0
    for i in range(probe,len(sample)):
        if sample[-i:] in target:
            a=i
        else:
            w=1
        if sample[:i] in target:
            b=i
        else:
            z=1
        if w and z:
            break
    if a:
        x=target.find(sample[-a:])
    if b:
        y=target.find(sample[:b])
    return a,x,b,y

def override(func):
    class OverrideAction(argparse.Action):
        def __call__(self,parser,namespace,values,option_string):
            func()
            parser.exit()
    return OverrideAction

def version():
    print('\n  Project: '+script+'\n  Version: '+__version__+'\n  Latest update: '+last_update+'\n  Author: '+author+'\n  License: '+license+'\n')

if __name__ == '__main__':
    main()
