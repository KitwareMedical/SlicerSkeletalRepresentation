# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 13:57:11 2014

@author: jvicory
"""

import srep
import datetime

def readSrepFromM3D(filename):
    f = open(filename,'r')
    
    lines = [line.strip() for line in f.readlines()]
    
    figidx = lines.index('figure[0] {')
    coloridx = lines.index('color {')
    endidx = lines.index('}',coloridx)

    figparams = getSectionDict(lines,figidx+1,coloridx)
    numrows = int(figparams['numRows'])
    numcols = int(figparams['numColumns'])

    fig = srep.figure(numrows,numcols)

    primidx = endidx + 1;
    
    for row in range(numrows):
        for col in range(numcols):
            endprimidx = lines.index('}',primidx)
            primsection = getSectionDict(lines,primidx+1,endprimidx)
            fig.addAtomFromDict(row,col,primsection)
            primidx = endprimidx + 1;
    
    return fig
    
def writeSrepToM3D(filename, fig):    
    f = open(filename,'w')
    
    f.write('pabloVersion = 9974 2009/07/24 19:36:23;\n'.expandtabs(4));
    f.write('coordSystem {\n'.expandtabs(4));
    f.write('\tyDirection = 1;\n'.expandtabs(4));
    f.write('}\n'.expandtabs(4));
    
    f.write('model {\n'.expandtabs(4));
    f.write('\tfigureCount = 1;\n'.expandtabs(4));
    f.write('\tname = test;\n'.expandtabs(4));
    f.write('\tfigureTrees {\n'.expandtabs(4));
    f.write('\t\tcount=1;\n'.expandtabs(4));
    f.write('\t\ttree[0] {\n'.expandtabs(4));
    f.write('\t\t\tattachmentMode = 0;\n'.expandtabs(4));
    f.write('\t\t\tblendAmount = 0;\n'.expandtabs(4));    
    f.write('\t\t\tblendExtent = 0;\n'.expandtabs(4));
    f.write('\t\t\tchildCount = 0;\n'.expandtabs(4));
    f.write('\t\t\tfigureId = 0;\n'.expandtabs(4));
    f.write('\t\t\tlinkCount = 0;\n'.expandtabs(4));
    f.write('\t\t}\n'.expandtabs(4));
    f.write('\t}\n'.expandtabs(4));
    
    f.write('\tfigure[0] {\n'.expandtabs(4));
    f.write('\t\tname = test;\n'.expandtabs(4));
    
    f.write('\t\tnumColumns = %s;\n'.expandtabs(4) % fig.numCols);
    
    f.write('\t\tnumLandmarks = 0;\n'.expandtabs(4));
    
    f.write('\t\tnumRows = %s;\n'.expandtabs(4) % fig.numRows);
    
    f.write('\t\tpositivePolarity = 1;\n'.expandtabs(4));
    f.write('\t\tpositiveSpace = 1;\n'.expandtabs(4));
    f.write('\t\tsmoothness = 50;\n'.expandtabs(4));
    f.write('\t\ttype = QuadFigure;\n'.expandtabs(4));
    
    f.write('\t\tcolor {\n'.expandtabs(4));
    f.write('\t\t\tblue = 0;\n'.expandtabs(4));
    f.write('\t\t\tgreen = 1;\n'.expandtabs(4));
    f.write('\t\t\tred = 0;\n'.expandtabs(4));
    f.write('\t\t}\n'.expandtabs(4));
    
    for row in range(fig.numRows):
        for col in range(fig.numCols):
            
            isCrest = fig.atoms[row,col].isCrest();
            f.write('\t\tprimitive[{0}][{1}] {{\n'.format(row,col).expandtabs(4));

            f.write('\t\t\tr[0] = {0};\n'.format(fig.atoms[row,col].topSpoke.r).expandtabs(4));
            f.write('\t\t\tr[1] = {0};\n'.format(fig.atoms[row,col].botSpoke.r).expandtabs(4));
            if isCrest:
                f.write('\t\t\tr[2] = {0};\n'.format(fig.atoms[row,col].crestSpoke.r).expandtabs(4));
            
            f.write('\t\t\tselected = 1;\n'.format(fig.atoms[row,col].selected).expandtabs(4));
            
            if isCrest:
                f.write('\t\t\ttype = EndPrimitive;\n'.expandtabs(4));
            else:
                f.write('\t\t\ttype = StandardPrimitive;\n'.expandtabs(4));
            
            f.write('\t\t\tux[0] = {0};\n'.format(fig.atoms[row,col].topSpoke.U[0]).expandtabs(4));
            f.write('\t\t\tux[1] = {0};\n'.format(fig.atoms[row,col].botSpoke.U[0]).expandtabs(4));
            
            if isCrest:
                f.write('\t\t\tux[2] = {0};\n'.format(fig.atoms[row,col].crestSpoke.U[0]).expandtabs(4));
            else:
                f.write('\t\t\tux[2] = 1;\n'.expandtabs(4));
                
            f.write('\t\t\tuy[0] = {0};\n'.format(fig.atoms[row,col].topSpoke.U[1]).expandtabs(4));
            f.write('\t\t\tuy[1] = {0};\n'.format(fig.atoms[row,col].botSpoke.U[1]).expandtabs(4));
            
            if isCrest:
                f.write('\t\t\tuy[2] = {0};\n'.format(fig.atoms[row,col].crestSpoke.U[1]).expandtabs(4));
            else:
                f.write('\t\t\tuy[2] = 0;\n'.expandtabs(4));
                
            f.write('\t\t\tuz[0] = {0};\n'.format(fig.atoms[row,col].topSpoke.U[2]).expandtabs(4));
            f.write('\t\t\tuz[1] = {0};\n'.format(fig.atoms[row,col].botSpoke.U[2]).expandtabs(4));
            
            if isCrest:
                f.write('\t\t\tuz[2] = {0};\n'.format(fig.atoms[row,col].crestSpoke.U[2]).expandtabs(4));
            else:
                f.write('\t\t\tuz[2] = 0;\n'.expandtabs(4));
                
            f.write('\t\t\tx = {0};\n'.format(fig.atoms[row,col].hub.P[0]).expandtabs(4));
            f.write('\t\t\ty = {0};\n'.format(fig.atoms[row,col].hub.P[1]).expandtabs(4));
            f.write('\t\t\tz = {0};\n'.format(fig.atoms[row,col].hub.P[2]).expandtabs(4));
            
            f.write('\t\t}\n'.expandtabs(4));
            
    f.write('\t}\n'.expandtabs(4));
    f.write('}');
            
            
    
    f.close()
    
    
    
def getSectionDict(lines, start, stop):
    section = [line.strip(';') for line in lines[start:stop]]
    sectiondict = dict(line.split(' = ') for line in section)
    return sectiondict
    
    