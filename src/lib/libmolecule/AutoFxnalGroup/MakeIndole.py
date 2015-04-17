import itertools
import io
Groups=["Amine","Alkenyl","Alkenyl","Alkenyl","Alkenyl","Alkenyl","Alkenyl","Alkenyl","Alkenyl"]
Values=[2,2,2,3,2,2,2,2,3]
CanChange=[0,1,2,4,5,6,7]
f1=open('Indole.hh','w')
f2=open('IndoleTypes.hh','w')
f3=open('PrintIndole.hh','w')
def PrintCarbon(number,Label1,NewValues,index,f1):
   f1.write("F_t::C"+str(number)+Label1.upper()+",")
   if NewValues[number]==2:
      f1.write("F_t::H"+str(index)+Label1.upper()+",")
      return index+1
   return index

def PrintTypes(Label1,f2):
    f2.write(Label1.upper()+",\n")
    f2.write("N"+Label1.upper()+",")
    for i in range(1,9):
       f2.write("C"+str(i)+Label1.upper()+",")
    for i in range(1,8):
       f2.write("H"+str(i)+Label1.upper()+",")
    f2.write("\n")

def TypeCases(Label1,f3):
    BaseLabel="\""
    Subs=Label1.split('_')
    if len(Subs)>=2:
       for i in Subs[1]:
          BaseLabel+=str(i)
          if not i==Subs[1][len(Subs[1])-1]:
             BaseLabel+=","
       BaseLabel+=" Substituted Indolyl"
    elif Label1=="Indolyl0":
       BaseLabel+="Indole"
    else:
       BaseLabel+="Heptasubstituted Indole"
    f3.write("CASE("+Label1.upper()+","+BaseLabel+"\")\n")
    f3.write("CASE(N"+Label1.upper()+","+BaseLabel+" Nitrogen\")\n")
    for i in range(1,9):
       f3.write("CASE(C"+str(i)+Label1.upper()+","+BaseLabel+" Carbon "+str(i)+"\")\n")
    for i in range(1,8):
       f3.write("CASE(H"+str(i)+Label1.upper()+","+BaseLabel+" Hydrogen "+str(i)+"\")\n")
 


for i in range(0,8):
   c=itertools.combinations(CanChange,i)
   for j in c:
       Label1="Indolyl"+str(i)
       if not i==0 and not i==7:
          Label1+="_"
       NewValues=Values[:]
       for k in j:
          if not i==0 and not i==7:
             Label1+=str(k+1)
          NewValues[k]+=1
       #print Label1+",",
       TypeCases(Label1,f3)
       PrintTypes(Label1,f2)
       f1.write("class "+Label1+" : public RingFinder<")
       Label=""
       for k in range(0,len(NewValues)):
           Label+=Groups[k]+str(NewValues[k])
           if not k==len(NewValues)-1:
              Label+=","
       f1.write(Label+">{\n")
       f1.write("   private:\n")
       f1.write("      typedef RingFinder<"+Label+"> Base_t;\n")
       f1.write("   public:\n")
       f1.write("   "+Label1+"():Base_t(")
       f1.write("F_t::"+Label1.upper()+",")
       f1.write("F_t::N"+Label1.upper()+","),
       index=1
       if NewValues[0]==2:
          f1.write("F_t::H"+str(index)+Label1.upper()+",")
          index+=1
       for C1 in range(1,8):
           index=PrintCarbon(C1,Label1,NewValues,index,f1)
       f1.write("F_t::C8"+Label1.upper()+"){}\n};\n")
