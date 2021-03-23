
from general_tools import pload



data=pload('93J.p')
out_folder='template_data/sn1993J/'

phases=open(out_folder+'Spec.JD.phases','w')
phases.write('#\n#\n')
for fname in data.keys():
    p=float(fname.split('_')[2])
    specname=fname[:-13]
    phases.write('%s\tfill\t%f\n' %(specname,p))
    var=open(out_folder+fname,'w')    
    spec=open(out_folder+specname,'w')
    for line in data[fname]:
        a,b,c=line.split()
        var.write(line)
        spec.write('%s\t%g\n' %(a,float(c)-float(b)))
    var.close()
    spec.close()
phases.close()